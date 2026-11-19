# =============================================================================
# Script  : currents_map.R
# Section : 2.1 Study area
# Purpose : Data-based map of the regional circulation around the three HARP
#           sites: long-term mean surface currents from the Copernicus GLORYS12
#           reanalysis (monthly means averaged over the climatology window)
#           drawn as vectors on the GMRT bathymetry base map, with the
#           climatological Antarctic Circumpolar Current (ACC) fronts of
#           Orsi et al. (1995). Reproducible replacement for the manual
#           PowerPoint currents/water-mass overlay.
# Inputs  : Copernicus Marine GLOBAL_MULTIYEAR_PHY_001_030, layer
#             cmems_mod_glo_phy_my_0.083deg_P1M-m (monthly), variables uo, vo
#             at the surface level (0.494 m); downloaded on first run (needs a
#             Copernicus Marine account; credentials as environment variables
#             COPERNICUSMARINE_SERVICE_USERNAME / _PASSWORD) and cached in
#             DOWNLOAD_DIR/copernicus_regional/
#           Environmental Data/Copernicus/Physics Reanalysis/EI40km_*.nc
#             (original daily portal download; used only to validate the
#             regional download)
#           Environmental Data/GMRT/GMRTbathymetry_highres_widerView.grd
#           R package orsifronts (0.2.0; front positions of Orsi et al. 1995)
#           final_scripts_odontocetes/config.R, 01_study_area/map_helpers.R,
#           05_environmental_data/copernicus_arco_helpers.R
# Outputs : DOWNLOAD_DIR/copernicus_regional/phy_surface_currents_P1M_<years>.nc
#           outputs/01_study_area/currents_map.png (600 dpi)
#           outputs/01_study_area/currents_map.pdf
#           outputs/01_study_area/currents_download_validation.csv
# Author  : Natalie Posdaljian
# Date    : 2026-11-17
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Averaging window = climatology window of the project: earliest
#             deployment year - CLIM_HALF_WINDOW_YEARS to latest deployment
#             year + CLIM_HALF_WINDOW_YEARS (2009-2021), all months weighted
#             equally.
#           - Front positions are climatological (hydrographic sections up to
#             ~1990) and are shown as published, not recomputed from the model.
#           - Water-mass and current names from the earlier manual overlay
#             (mCDW, TWW/WSW, BC/TBW/BSW, Antarctic Coastal Current, "PF") are
#             NOT drawn until their positions are taken from a cited source;
#             see section 7.
#           - Data citations: E.U. Copernicus Marine Service Information,
#             GLOBAL_MULTIYEAR_PHY_001_030, doi:10.48670/moi-00021
#             [CONFIRM product version (layer suffix 202311) and access date];
#             Orsi et al. (1995) Deep-Sea Res. I 42:641-673.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "01_study_area", "map_helpers.R"))
source(file.path(FINAL_DIR, "05_environmental_data", "copernicus_arco_helpers.R"))

suppressPackageStartupMessages({
  library(orsifronts)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
currents_layer <- "cmems_mod_glo_phy_my_0.083deg_P1M-m"   # monthly means
current_vars   <- c("uo", "vo")
grid_res_deg   <- 1 / 12

# Download region (xmin, ymin, xmax, ymax): one degree wider than the map in
# longitude and half a degree in latitude, so vectors reach the map edges.
download_region <- c(MAP_XLIM[1] - 1, MAP_YLIM[1] - 0.5,
                     MAP_XLIM[2] + 1, MAP_YLIM[2] + 0.5)

clim_years <- c(min(year(SITES$rec_start)) - CLIM_HALF_WINDOW_YEARS,
                max(year(SITES$rec_end)) + CLIM_HALF_WINDOW_YEARS)

regional_dir  <- file.path(DOWNLOAD_DIR, "copernicus_regional")
currents_file <- file.path(regional_dir,
                           sprintf("phy_surface_currents_P1M_%d-%d.nc",
                                   clim_years[1], clim_years[2]))

# Vector display: one arrow every `arrow_stride` model cells (1/12 deg), and
# arrow length = `arrow_deg_per_ms` degrees of latitude per 1 m/s (the
# longitude component is stretched by 1/cos(latitude) so arrows keep their
# true direction and length on the map).
arrow_stride      <- 4
arrow_deg_per_ms  <- 1.5
reference_speed   <- 0.2   # m/s, legend arrow

# -----------------------------------------------------------------------------
# 2. Download monthly surface currents (first run only)
# -----------------------------------------------------------------------------
# The cached NetCDF is re-used on later runs. The download goes through the
# workarounds in copernicus_arco_helpers.R (packed int16 scaling, longitude
# offset); vertical range 0 to -1 m returns only the 0.494 m surface level.
download_currents <- function(file) {
  check_copernicus_setup()
  product <- COPERNICUS$phy$product
  layer_years <- year(layer_time_extent(product, currents_layer))
  if (clim_years[1] < layer_years[1] || clim_years[2] > layer_years[2]) {
    stop("Averaging window ", paste(clim_years, collapse = "-"),
         " is not fully covered by layer ", currents_layer)
  }
  attrs <- layer_zarr_attrs(product, currents_layer, current_vars)
  months <- seq(as.Date(sprintf("%d-01-01", clim_years[1])),
                as.Date(sprintf("%d-12-01", clim_years[2])), by = "month")

  # One request per month: with CopernicusMarine 0.3.2, a request spanning
  # several time chunks of this layer (one month per chunk) comes back with
  # the cells scrambled across time and latitude, whereas single-time
  # requests are correct (checked against the land mask and against the
  # portal download in section 4).
  fields <- lapply(setNames(current_vars, current_vars), function(v) {
    message("Downloading monthly ", v, " ", paste(clim_years, collapse = "-"))
    per_month <- lapply(months, function(m) {
      t_m <- as.POSIXct(rep(format(m), 2), tz = "UTC")
      x <- download_subset(product, currents_layer, v, download_region, t_m,
                           verticalrange = c(0, -1))
      subset_to_array(x, v, attrs[[v]], grid_res_deg,
                      plausible_range = c(-3, 3))
    })
    first <- per_month[[1]]
    same_grid <- vapply(per_month, function(a) {
      identical(a$lon, first$lon) && identical(a$lat, first$lat)
    }, logical(1))
    if (!all(same_grid)) stop(v, ": monthly grids differ")
    list(values = simplify2array(lapply(per_month, function(a) a$values[, , 1, 1])),
         lon = first$lon, lat = first$lat, depth = first$depth,
         time = do.call(c, lapply(per_month, `[[`, "time")))
  })
  ref <- fields[[1]]

  # Land cells (NA) must be identical in every month; a changing mask would
  # indicate mis-assembled chunks.
  na_mask <- is.na(ref$values)
  if (any(na_mask != as.vector(na_mask[, , 1]))) {
    stop("Land mask differs between months: download is mis-assembled")
  }

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  nc_dims <- list(
    ncdim_def("longitude", "degrees_east", ref$lon),
    ncdim_def("latitude", "degrees_north", ref$lat),
    ncdim_def("time", "seconds since 1970-01-01 00:00:00", as.numeric(ref$time))
  )
  nc_vars <- lapply(current_vars, function(v) {
    ncvar_def(v, units = attrs[[v]]$zattrs$units, dim = nc_dims,
              missval = NaN, prec = "float",
              longname = paste(attrs[[v]]$zattrs$standard_name,
                               "at", ref$depth, "m"))
  })
  tmp <- paste0(file, ".part")
  nc <- nc_create(tmp, nc_vars)
  for (i in seq_along(current_vars)) {
    ncvar_put(nc, nc_vars[[i]], fields[[current_vars[i]]]$values)
  }
  ncatt_put(nc, 0, "source", paste(COPERNICUS$phy$product, currents_layer))
  ncatt_put(nc, 0, "depth_m", ref$depth)
  nc_close(nc)
  file.rename(tmp, file)
}

if (!file.exists(currents_file)) download_currents(currents_file)

# -----------------------------------------------------------------------------
# 3. Read the regional monthly fields
# -----------------------------------------------------------------------------
read_currents <- function(file) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  lon  <- ncvar_get(nc, "longitude")
  lat  <- ncvar_get(nc, "latitude")
  time <- as.POSIXct(ncvar_get(nc, "time"), origin = "1970-01-01", tz = "UTC")
  grid <- expand_grid(time = time, lat = lat, lon = lon)
  # ncdf4 returns (lon, lat, time); as.vector runs lon fastest, then lat,
  # then time, matching the expand_grid order above (last column fastest).
  grid |>
    mutate(month = as.Date(time),
           uo = as.vector(ncvar_get(nc, "uo")),
           vo = as.vector(ncvar_get(nc, "vo"))) |>
    select(month, lon, lat, uo, vo)
}

currents_monthly <- read_currents(currents_file)

# -----------------------------------------------------------------------------
# 4. Validate the download against the original portal file
# -----------------------------------------------------------------------------
# The EI +/- 40 km box was downloaded manually from the Copernicus portal as
# DAILY fields (2014-03-05 to 2014-07-17). For every complete calendar month
# in that file, the mean of the daily surface currents in each cell should
# equal the monthly layer in the same cell up to the int16 packing precision
# (~0.0006 m/s). This checks scaling, grid alignment (longitude offset), and
# variable order in one step.
validate_against_portal <- function(currents_monthly) {
  portal_file <- list.files(file.path(ENV_DIR, "Copernicus",
                                      "Physics Reanalysis"),
                            pattern = "^EI40km_cmems_mod_glo_phy",
                            full.names = TRUE)
  nc <- nc_open(portal_file)
  on.exit(nc_close(nc))
  lon  <- round(ncvar_get(nc, "longitude"), 4)
  lat  <- round(ncvar_get(nc, "latitude"), 4)
  date <- as.Date(as.POSIXct(ncvar_get(nc, "time"), origin = "1970-01-01",
                             tz = "UTC"))
  surface <- function(v) {
    ncvar_get(nc, v, start = c(1, 1, 1, 1),
              count = c(-1, -1, 1, -1), collapse_degen = FALSE)
  }
  daily <- expand_grid(date = date, lat = lat, lon = lon) |>
    mutate(uo_daily = as.vector(surface("uo")),
           vo_daily = as.vector(surface("vo")))

  complete_months <- tibble(date = date) |>
    count(month = floor_date(date, "month")) |>
    filter(n == days_in_month(month)) |>
    pull(month)

  daily |>
    mutate(month = floor_date(date, "month")) |>
    filter(month %in% complete_months) |>
    group_by(month, lon, lat) |>
    summarise(across(c(uo_daily, vo_daily), mean), .groups = "drop") |>
    inner_join(currents_monthly |>
                 mutate(lon = round(lon, 4), lat = round(lat, 4)),
               by = c("month", "lon", "lat")) |>
    filter(!is.na(uo_daily), !is.na(uo)) |>
    group_by(month) |>
    summarise(n_cells     = n(),
              max_abs_du  = max(abs(uo - uo_daily)),
              max_abs_dv  = max(abs(vo - vo_daily)),
              cor_u       = cor(uo, uo_daily),
              cor_v       = cor(vo, vo_daily),
              .groups = "drop")
}

validation <- validate_against_portal(currents_monthly)
print(validation)
write_csv(validation, out_path("01_study_area",
                               "currents_download_validation.csv"))
if (nrow(validation) == 0 ||
      any(validation$max_abs_du > 0.005 | validation$max_abs_dv > 0.005)) {
  stop("Regional monthly currents do not match the original portal download")
}

# -----------------------------------------------------------------------------
# 5. Long-term mean surface currents and arrow geometry
# -----------------------------------------------------------------------------
# Mean over all months of the window; cells that are land (NA) in any month
# are dropped. Arrows are thinned to every `arrow_stride`-th cell in both
# directions.
mean_currents <- currents_monthly |>
  group_by(lon, lat) |>
  summarise(n_months = sum(!is.na(uo)),
            uo = mean(uo), vo = mean(vo), .groups = "drop") |>
  filter(n_months == max(n_months)) |>
  mutate(speed = sqrt(uo^2 + vo^2))

lon_keep <- sort(unique(mean_currents$lon))[c(TRUE, rep(FALSE, arrow_stride - 1))]
lat_keep <- sort(unique(mean_currents$lat))[c(TRUE, rep(FALSE, arrow_stride - 1))]

arrows_df <- mean_currents |>
  filter(lon %in% lon_keep, lat %in% lat_keep) |>
  mutate(lon_end = lon + uo * arrow_deg_per_ms / cos(lat * pi / 180),
         lat_end = lat + vo * arrow_deg_per_ms)

message(sprintf("Mean surface speed in map extent: median %.3f m/s, max %.3f m/s",
                median(mean_currents$speed), max(mean_currents$speed)))

# Reference-arrow key in a white box in the top-right corner of the map:
# label on top, arrow of `reference_speed` centred below it.
key_box <- c(xmin = -52.35, xmax = -51.05, ymin = -60.42, ymax = -60.04)
ref_len_lon <- reference_speed * arrow_deg_per_ms / cos(-60.3 * pi / 180)
ref_arrow <- tibble(lon = mean(key_box[c("xmin", "xmax")]) - ref_len_lon / 2,
                    lat = -60.32,
                    lon_end = lon + ref_len_lon,
                    lat_end = lat)
ref_label <- paste0(reference_speed, " m s", intToUtf8(c(8315, 185)))

# -----------------------------------------------------------------------------
# 6. ACC fronts (Orsi et al. 1995)
# -----------------------------------------------------------------------------
# Only the fronts that cross the map are kept: the Southern ACC Front (SACCF)
# and the Southern Boundary of the ACC (SB). The Polar Front lies north of
# the map (~58-59 S in Drake Passage).
front_names <- c(saccf = "Southern ACC Front (SACCF)",
                 sbdy  = "Southern Boundary of the ACC (SB)")

fronts <- st_as_sf(orsifronts) |>
  st_make_valid() |>
  filter(front %in% names(front_names)) |>
  st_set_agr("constant") |>
  st_crop(xmin = download_region[1], ymin = download_region[2],
          xmax = download_region[3], ymax = download_region[4]) |>
  mutate(front_label = factor(front_names[as.character(front)],
                              levels = front_names))

# -----------------------------------------------------------------------------
# 7. Optional water-mass / current labels (from a cited source only)
# -----------------------------------------------------------------------------
# The previous manual overlay labelled mCDW, TWW/WSW, BC/TBW/BSW, the
# Antarctic Coastal Current, and "PF" (probably the Peninsula Front in
# Bransfield Strait, not the Polar Front). Labels are drawn only for rows
# whose position has been filled in from a cited source.
# [INSERT positions (lon, lat) and citations for each label, e.g., from
#  Sangra et al. (2011, 2017) for Bransfield Strait or Thompson et al. (2009)
#  for the Antarctic Slope Front / Coastal Current]
water_mass_labels <- tribble(
  ~label,                       ~lon,         ~lat,         ~citation,
  "mCDW",                       NA_real_,     NA_real_,     "[INSERT]",
  "TWW/WSW",                    NA_real_,     NA_real_,     "[INSERT]",
  "BC/TBW/BSW",                 NA_real_,     NA_real_,     "[INSERT]",
  "Antarctic Coastal Current",  NA_real_,     NA_real_,     "[INSERT]"
) |>
  filter(!is.na(lon), !is.na(lat))

# -----------------------------------------------------------------------------
# 8. Map
# -----------------------------------------------------------------------------
# Bathymetry is slightly faded so the vectors stand out; each arrow is drawn twice
# (white halo, then black) for contrast over both light and dark water.
bathy <- load_bathymetry()
arrow_head <- arrow(length = unit(0.07, "cm"), type = "closed", angle = 25)
front_colors <- setNames(c("#D81B60", "#FF8F00"), front_names)
front_legend_title <- "ACC fronts\n(Orsi et al. 1995)"

currents_map <- ggplot() +
  bathymetry_layers(bathy, ocean_alpha = 0.8) +
  # ACC fronts
  geom_sf(data = fronts, aes(color = front_label, linetype = front_label),
          linewidth = 0.9) +
  scale_color_manual(values = c(SITE_COLORS, front_colors),
                     breaks = unname(front_names),
                     name = front_legend_title,
                     guide = guide_legend(order = 2)) +
  scale_linetype_manual(values = c("longdash", "dotdash"),
                        breaks = unname(front_names),
                        name = front_legend_title,
                        guide = guide_legend(order = 2)) +
  # Mean surface current vectors
  geom_segment(data = arrows_df,
               aes(x = lon, y = lat, xend = lon_end, yend = lat_end),
               color = "white", linewidth = 0.75, arrow = arrow_head) +
  geom_segment(data = arrows_df,
               aes(x = lon, y = lat, xend = lon_end, yend = lat_end),
               color = "black", linewidth = 0.3, arrow = arrow_head) +
  # Reference-arrow key
  annotate("rect", xmin = key_box[["xmin"]], xmax = key_box[["xmax"]],
           ymin = key_box[["ymin"]], ymax = key_box[["ymax"]],
           fill = alpha("white", 0.9), color = "grey40", linewidth = 0.2) +
  annotate("text", x = mean(key_box[c("xmin", "xmax")]), y = -60.14,
           label = ref_label, size = 3, color = INK_COLOR) +
  geom_segment(data = ref_arrow,
               aes(x = lon, y = lat, xend = lon_end, yend = lat_end),
               color = "black", linewidth = 0.4, arrow = arrow_head) +
  # Optional cited water-mass labels (none until positions are filled in)
  geom_text(data = water_mass_labels, aes(x = lon, y = lat, label = label),
            size = 3.2, fontface = "italic", color = INK_COLOR) +
  # Sites (drawn without their own colour scale entry: the shared colour
  # scale above carries both site and front colours)
  site_layers()[-2] +
  map_frame() +
  theme(legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.5, "cm"))

save_map(currents_map, "currents_map", width_in = 9.5, height_in = 5.6)
