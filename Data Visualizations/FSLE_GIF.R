library(tidyverse)
library(raster)
library(sf)
library(ggspatial)
library(ncdf4)
library(maps)
library(ggquiver)
library(metR)
library(magick)

# Code to visualize FSLE features for a particular day
# Will run this for the first of the month from 02/15 to 12/16
# Also running for month of july 2016 to visualize granular movements

# ---------------- Step 1: Set Dates and Site Info -------------------
# list desired dates 
# options: 03-05-2014 to 07-17-2014, 02-10-2015 to 1-29-2016, and 02-04-2016 to 12-02-2016
# monthly dates:
# dates <- c('02-01-2015','03-01-2015','04-01-2015','05-01-2015','06-01-2015','07-01-2015',
#            '08-01-2015','09-01-2015','10-01-2015','11-01-2015','12-01-2015','01-01-2016',
#            '02-01-2016','03-01-2016','04-01-2016','05-01-2016','06-01-2016','07-01-2016',
#            '08-01-2016','09-01-2016','10-01-2016','11-01-2016','12-01-2016')

# sequence of dates for July 2016:
start_date <- as.Date("07-01-2016",format = '%m-%d-%Y')
end_date <- as.Date("07-31-2016",format = '%m-%d-%Y')
dates <- seq(from = start_date, to = end_date, by = "day")

dates <- as.Date(dates, format = '%m-%d-%Y')

# get land data (to gray out land)
land_map <- map_data("world", region = 'Antarctica')
land_map <- land_map %>% filter(long >= (lon_min-360), long <= (lon_max-360),
                                lat >= lat_min, lat <= lat_max)
land_map$lon <- land_map$long + 360
land_map <- subset(land_map, select=-long)

# site coordinates
sites <- data.frame(
  name = c("EI", "KGI", "CI"),
  lon = c((-55.95400 + 360), (-57.941917 + 360), (-53.483433+ 360)),
  lat = c(-60.8869, -61.457817, -61.251867))
sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

# Define bounding box to cover all sites
lat_min <- -62
lat_max <- -60
lon_min <- 301
lon_max <- 308

# ---------------- Step 2: Load File ------------------
loadFSLE <- function(date) {
  # loading file for specified date
  date_str <- format(date, "%Y%m%d")
  FSLE_nc <- nc_open(paste0("D:/FSLE/AVISO_Global_FSLE/dt_global_allsat_madt_fsle_",
                            date_str,"_20180704.nc"))
  
  # Extract lat long and constrain to bounding box
  lon <- ncvar_get(FSLE_nc, 'lon')
  lat <- ncvar_get(FSLE_nc,'lat')
  date <- ncvar_get(FSLE_nc, 'time')
  
  # Find index ranges for the bounding box
  lat_idx <- which(lat >= lat_min & lat <= lat_max)
  lon_idx <- which(lon >= lon_min & lon <= lon_max)
  
  # extract FSLE magnitude and orientation
  fsle <- ncvar_get(FSLE_nc, "fsle_max", start = c(min(lon_idx), min(lat_idx),1),
                    count = c(length(lon_idx), length(lat_idx),1))
  fsle_orient <- ncvar_get(FSLE_nc, "theta_max", start = c(min(lon_idx), min(lat_idx),1),
                           count = c(length(lon_idx), length(lat_idx),1))
  
  # Create dataframe with FSLE values
  fsle_df <- expand.grid(lon = lon[lon_idx], lat = lat[lat_idx], date = date)
  fsle_df$fsle <- as.vector(fsle)
  fsle_df$fsle_orient <- as.vector(fsle_orient)
  fsle_df$date <- as.Date(fsle_df$date, origin = '1950-01-01')
  
  return(fsle_df)
}

# ---------------- Step 3: Plot Data ------------------
# Note: NOT going to plot FSLE orientations (only plotting magnitude)
#    Plots dare more difficult to visually interpret mesoscale features from
plotFSLE <- function(date) {
  fsle_df <- loadFSLE(date)
  date_str <- format(date, "%b %d %Y")
  
  plot <- ggplot() + 
    # fsle values
    geom_tile(data = fsle_df, aes(x = lon, y = lat, fill = fsle)) +
    
    # gray out land
    geom_polygon(data = land_map, aes(x = lon, y = lat, group = group), fill = "#22192d") + 
    
    # add in sites and site labels
    geom_sf(data = sites_sf, shape = 21, fill = "darkmagenta", size = 5, color = "#22192d") +
    geom_text(data = sites, aes(x = lon, y = lat, label = name),
              nudge_y = -0.12, nudge_x = 0.2, color = "#22192d", size = 7) +
    
    # style
    coord_sf(ylim=c(-62,-60), xlim=c(301,308),expand=FALSE) + 
    scale_fill_viridis_c(option = "inferno", name ='FSLE Strength', limits = c(-0.6, 0)) +
    labs(y=NULL,x=NULL,title=date_str) + 
    theme(axis.text = element_text(size=10, vjust=-0.2))
  
  return(plot)
}

# call and save plots
for(d in 1:length(dates)) {
  date_str <- format(dates[d], '%m-%d-%Y')
  png(paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Presentation Plots/FSLE Images/",
             date_str,'.png'))
  print(plotFSLE(dates[d]))
  dev.off()
  print(paste("Plot for ",dates[d],"done."))
}

# now, manually save plots from R studio viewer
# file title: month-year.png (e.g. 2-2015 for february 2015)
#      or month-day-year.png (e.g. 2-01-2015 for February 1 2015)

# --------------- Step 4: Make GIF --------------------
# adding images to list for GIF
# IMPORTANT: Manually edit 'frames' filepaths for the dates you are interested in animating

# monthly image files:
# frames <- paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Presentation Plots/FSLE Images/", 
#                  2:12, "-2015.png")
# frames <- append(frames, paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Presentation Plots/FSLE Images/", 
#                                      1:12, "-2016.png"))

# july 2016 front image files:
frames <- paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Presentation Plots/FSLE Images/07-0",
                 1:9, "-2016.png")
frames <- append(frames, paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Presentation Plots/FSLE Images/07-",
                                10:31, "-2016.png"))

# making GIF
gif <- image_read(frames)
gif <- image_animate(gif, fps=4)