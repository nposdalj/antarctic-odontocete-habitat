library(tidyverse)
library(raster)
library(sf)
library(ggspatial)
library(ncdf4)
library(maps)
library(ggquiver)
library(metR)
library(magick)

# Code to visualize FSLE magnitude + orientation features for a date range,
# animated into a two-panel GIF (magnitude on the left, orientation on the
# right).

# =====================================================================
# USER SETTINGS -- edit everything in this block for your own machine,
# date range, and sites. Nothing below "END USER SETTINGS" should need
# to change for normal use.
# =====================================================================

# ---- Date range to animate ----
# options: 03-05-2014 to 07-17-2014, 02-10-2015 to 1-29-2016, and 02-04-2016 to 12-02-2016
start_date <- as.Date("2014-03-05")
end_date   <- as.Date("2014-07-17")

# ---- Map bounding box (covers all sites below) ----
# NOTE: longitude here uses the 0-360 convention (matches the FSLE NetCDF
# files), not the usual -180/180 convention.
lat_min <- -62
lat_max <- -60
lon_min <- 301
lon_max <- 308

# ---- Site coordinates to plot on the map ----
# input as standard -180/180 longitude; +360 is applied automatically below
sites <- data.frame(
  name    = c("EI",       "KGI",        "CI"),
  lon_180 = c(-55.95400,  -57.941917,   -53.483433),
  lat     = c(-60.8869,   -61.457817,   -61.251867)
)

# ---- Input data paths ----
# folder of daily FSLE NetCDF files (fsle_max / theta_max variables)
fsle_folder <- "F:/Trisha's Old Stuff/FSLE/AVISO_Global_FSLE"

# GMRT bathymetry raster used for depth contour lines (see siteMap.R)
bathymetry_path <- paste0(
  "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/",
  "Environmental Data/GMRT/GMRTbathymetry_lowerres_widerView.grd"
)

# ---- Output paths ----
# daily frame PNGs are written here, then assembled into the GIF below
image_folder <- paste0(
  "F:/FSLEs/FSLE_images/"
)

gif_path <- file.path(image_folder, "FSLE_2014.gif")
gif_fps  <- 4

# ---- Output resolution ----
# pixel size of EACH panel (magnitude and orientation are each rendered at
# this size, then placed side by side -- so the combined frame/GIF is
# (2 * frame_width) x frame_height). Bump these up for a sharper GIF (bigger
# file size, slower to render); the source data resolution is unaffected --
# this only controls how many pixels the plot is rasterized to.
frame_width  <- 480
frame_height <- 480

# ---- Bathymetry contour style ----
bathy_contour_color     <- "black"
bathy_contour_linewidth <- 0.3

# depth contour lines are drawn every bathy_contour_interval meters, from
# -bathy_contour_interval down to -bathy_contour_max_depth
bathy_contour_interval  <- 500
bathy_contour_max_depth <- 4500

# which isobaths (in meters, positive) get a "<depth> m" text label directly
# on the line -- must also be one of the levels generated above (i.e. a
# multiple of bathy_contour_interval, at or shallower than bathy_contour_max_depth)
bathy_label_depths <- c(1000)
bathy_label_size   <- 3

# =====================================================================
# END USER SETTINGS
# =====================================================================


# ---------------- Derived Setup -------------------

dates <- seq(from = start_date, to = end_date, by = "day")

if (!dir.exists(image_folder)) {
  dir.create(image_folder, recursive = TRUE)
}

# get land data (to gray out land)
land_map <- map_data("world", region = 'Antarctica')
land_map <- land_map %>% filter(long >= (lon_min-360), long <= (lon_max-360),
                                lat >= lat_min, lat <= lat_max)
land_map$lon <- land_map$long + 360
land_map <- subset(land_map, select=-long)

# site coordinates, shifted to the 0-360 longitude convention used by FSLE data
sites$lon <- sites$lon_180 + 360
sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

# ---------------- Bathymetry Contours (see siteMap.R) -------------------
# Vector depth-contour lines drawn straight into the ggplot, same
# raster -> rasterToContour -> geom_sf pattern used in siteMap.R,
# replacing the old magick::image_composite PNG overlay approach.

if (!file.exists(bathymetry_path)) {
  stop("Bathymetry file not found at: ", bathymetry_path)
}

bathymetry <- raster(bathymetry_path)

# crop to the FSLE bounding box (bathymetry raster uses -180/180 longitude,
# FSLE data uses 0-360, so shift the crop box down by 360 to match)
bathy_bbox <- extent(lon_min - 360, lon_max - 360, lat_min, lat_max)
bathymetry <- crop(bathymetry, bathy_bbox)

# depth contour lines, spaced every bathy_contour_interval meters
bathy_contour_levels <- -seq(
  bathy_contour_interval,
  bathy_contour_max_depth,
  by = bathy_contour_interval
)

bathy_contours <- rasterToContour(bathymetry, levels = bathy_contour_levels)
bathy_contours_sf <- st_as_sf(bathy_contours)
bathy_contours_sf$level_num <- as.numeric(as.character(bathy_contours_sf$level))

# shift contour longitudes back to the FSLE plots' 0-360 convention
bathy_contours_sf <- st_shift_longitude(bathy_contours_sf)

# midpoint of each labeled isobath's line, used to place a "<depth> m" label
bathy_labels_sf <- bathy_contours_sf[bathy_contours_sf$level_num %in% -bathy_label_depths, ]
if (nrow(bathy_labels_sf) > 0) {
  bathy_labels_sf <- st_set_geometry(
    bathy_labels_sf,
    st_point_on_surface(st_geometry(bathy_labels_sf))
  )
}


# ---------------- Step 2: Load File ------------------

if (!dir.exists(fsle_folder)) {
  stop("fsle_folder does not exist: ", fsle_folder)
}

fsle_files <- list.files(
  path = fsle_folder,
  pattern = "\\.nc$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(fsle_files) == 0) {
  stop(
    "No .nc files found under fsle_folder: ", fsle_folder,
    "\nDouble-check the path and that the files are actually .nc."
  )
}

file_dates <- as.Date(
  sub(
    ".*fsle_([0-9]{8})_.*",
    "\\1",
    basename(fsle_files)
  ),
  format = "%Y%m%d"
)

if (any(is.na(file_dates))) {
  stop(
    "Could not parse a date from these filenames (expected an 8-digit date ",
    "after 'fsle_'):\n",
    paste(basename(fsle_files)[is.na(file_dates)], collapse = "\n")
  )
}

# Check whether every desired date has a corresponding FSLE file
missing_dates <- dates[!dates %in% file_dates]

if (length(missing_dates) > 0) {
  stop(
    "Missing FSLE files for:\n",
    paste(missing_dates, collapse = ", ")
  )
}

loadFSLE <- function(date) {

  # Find the file corresponding to this date
  file_index <- which(file_dates == date)

  if (length(file_index) == 0) {
    stop("No FSLE file found for date: ", date)
  }

  if (length(file_index) > 1) {
    stop("Multiple FSLE files found for date: ", date)
  }

  # Open the matching NetCDF file
  FSLE_nc <- nc_open(fsle_files[file_index])

  # Extract longitude, latitude, and time
  lon <- ncvar_get(FSLE_nc, "lon")
  lat <- ncvar_get(FSLE_nc, "lat")
  time <- ncvar_get(FSLE_nc, "time")

  # Find index ranges for bounding box
  lat_idx <- which(lat >= lat_min & lat <= lat_max)
  lon_idx <- which(lon >= lon_min & lon <= lon_max)

  # Extract FSLE magnitude
  fsle <- ncvar_get(
    FSLE_nc,
    "fsle_max",
    start = c(min(lon_idx), min(lat_idx), 1),
    count = c(length(lon_idx), length(lat_idx), 1)
  )

  # Extract FSLE orientation
  fsle_orient <- ncvar_get(
    FSLE_nc,
    "theta_max",
    start = c(min(lon_idx), min(lat_idx), 1),
    count = c(length(lon_idx), length(lat_idx), 1)
  )

  # Create data frame
  fsle_df <- expand.grid(
    lon = lon[lon_idx],
    lat = lat[lat_idx],
    date = time
  )

  fsle_df$fsle <- as.vector(fsle)
  fsle_df$fsle_orient <- as.vector(fsle_orient)

  # Convert NetCDF time to Date
  fsle_df$date <- as.Date(
    fsle_df$date,
    origin = "1950-01-01"
  )

  # Close NetCDF file
  nc_close(FSLE_nc)

  return(fsle_df)
}

# ---------------- Step 3: Plot Data ------------------
# Magnitude (fsle_max) plotted with the inferno scale, orientation (theta_max)
# plotted separately below with a cyclic rainbow colormap, since orientation
# is axial data: -90 deg and 90 deg represent the same line orientation.
plotFSLE <- function(date) {
  fsle_df <- loadFSLE(date)
  date_str <- format(date, "%b %d %Y")

  plot <- ggplot() +
    # fsle values
    geom_tile(data = fsle_df, aes(x = lon, y = lat, fill = fsle)) +

    # gray out land
    geom_polygon(data = land_map, aes(x = lon, y = lat, group = group), fill = "#22192d") +

    # bathymetry depth contour lines (see siteMap.R)
    geom_sf(data = bathy_contours_sf, color = bathy_contour_color, linewidth = bathy_contour_linewidth) +

    # isobath depth labels (bathy_label_depths)
    geom_sf_text(data = bathy_labels_sf, aes(label = paste0(abs(level_num), " m")),
                 color = bathy_contour_color, size = bathy_label_size, fontface = "bold") +

    # add in sites and site labels
    geom_sf(data = sites_sf, shape = 21, fill = "darkmagenta", size = 5, color = "#22192d") +
    geom_text(data = sites, aes(x = lon, y = lat, label = name),
              nudge_y = -0.12, nudge_x = 0.2, color = "#22192d", size = 7) +

    # style
    coord_sf(ylim=c(lat_min, lat_max), xlim=c(lon_min,lon_max),expand=FALSE) +
    scale_fill_viridis_c(option = "inferno", name ='FSLE Strength', limits = c(-0.6, 0)) +
    labs(y=NULL,x=NULL,title=date_str) +
    theme(axis.text = element_text(size=10, vjust=-0.2))

  return(plot)
}

# Cyclic rainbow palette for orientation: -90 and 90 both map to red since
# orientation is axial (a line at -90 deg and one at 90 deg point the same way)
orientation_colors <- c("red", "yellow", "green", "cyan", "blue", "magenta", "red")

# Orientation (theta_max) map, R port of Theta_Map_lol_2.m
plotTheta <- function(date) {
  fsle_df <- loadFSLE(date)
  date_str <- format(date, "%b %d %Y")

  plot <- ggplot() +
    # fsle orientation values
    geom_tile(data = fsle_df, aes(x = lon, y = lat, fill = fsle_orient)) +

    # gray out land
    geom_polygon(data = land_map, aes(x = lon, y = lat, group = group), fill = "#22192d") +

    # bathymetry depth contour lines (see siteMap.R)
    geom_sf(data = bathy_contours_sf, color = bathy_contour_color, linewidth = bathy_contour_linewidth) +

    # isobath depth labels (bathy_label_depths)
    geom_sf_text(data = bathy_labels_sf, aes(label = paste0(abs(level_num), " m")),
                 color = bathy_contour_color, size = bathy_label_size, fontface = "bold") +

    # add in sites and site labels
    geom_sf(data = sites_sf, shape = 21, fill = "darkmagenta", size = 5, color = "#22192d") +
    geom_text(data = sites, aes(x = lon, y = lat, label = name),
              nudge_y = -0.12, nudge_x = 0.2, color = "#22192d", size = 7) +

    # style
    coord_sf(ylim=c(lat_min, lat_max), xlim=c(lon_min,lon_max),expand=FALSE) +
    scale_fill_gradientn(
      colors = orientation_colors,
      name = "FSLE Orientation (°)",
      limits = c(-90, 90),
      breaks = c(-90, -45, 0, 45, 90),
      labels = c("-90 E <-> W", "-45 SW <-> NE", "0 N <-> S", "45 NW <-> SE", "90 E <-> W")
    ) +
    labs(y=NULL,x=NULL,title=date_str) +
    theme(axis.text = element_text(size=10, vjust=-0.2))

  return(plot)
}

# ---------------- Step 4: Plot and Save Frames ------------------
# Create each daily frame: magnitude panel (left) + orientation panel (right)
for (d in seq_along(dates)) {

  date_str <- format(dates[d], "%m-%d-%Y")

  frame_path       <- file.path(image_folder, paste0(date_str, ".png"))
  fsle_frame_path  <- file.path(image_folder, paste0(date_str, "_fsle.png"))
  theta_frame_path <- file.path(image_folder, paste0(date_str, "_theta.png"))

  # ---- Magnitude panel ----
  png(
    filename = fsle_frame_path,
    width = frame_width,
    height = frame_height,
    units = "px"
  )

  print(plotFSLE(dates[d]))

  dev.off()

  fsle_frame <- image_read(fsle_frame_path)

  # ---- Orientation panel ----
  png(
    filename = theta_frame_path,
    width = frame_width,
    height = frame_height,
    units = "px"
  )

  print(plotTheta(dates[d]))

  dev.off()

  theta_frame <- image_read(theta_frame_path)

  # ---- Combine side-by-side (magnitude | orientation) ----
  combined_frame <- image_append(c(fsle_frame, theta_frame))

  image_write(
    combined_frame,
    path = frame_path
  )

  # clean up the intermediate single-panel files
  file.remove(fsle_frame_path, theta_frame_path)

  print(paste("Combined magnitude/orientation plot for", dates[d], "done."))
}

# --------------- Step 5: Make GIF --------------------

frames <- file.path(
  image_folder,
  paste0(format(dates, "%m-%d-%Y"), ".png")
)

# Check that all expected PNG frames exist
missing_frames <- frames[!file.exists(frames)]

if (length(missing_frames) > 0) {
  stop(
    "The following PNG frames are missing:\n",
    paste(missing_frames, collapse = "\n")
  )
}

# Read images
gif <- image_read(frames)

# Animate GIF
gif <- image_animate(
  gif,
  fps = gif_fps
)

# Save GIF
image_write(
  gif,
  path = gif_path
)
