library(tidyverse)
library(raster)
library(sf)
library(ggspatial)
library(ncdf4)

# Code to visualize FSLE features for a particular day

# ---------------- Step 1: Load & Prep Data -------------------
# loading file for some arbitrary day
FSLE_nc <- nc_open("D:/FSLE/AVISO_Global_FSLE/dt_global_allsat_madt_fsle_20151027_20180704.nc")

# Extract lat long and constrain to bounding box
lon <- ncvar_get(FSLE_nc, 'lon')
lat <- ncvar_get(FSLE_nc,'lat')
date <- ncvar_get(FSLE_nc, 'time')

# Define bounding box to cover all sites
lat_min <- -62
lat_max <- -60
lon_min <- 301
lon_max <- 308
# Find index ranges for the bounding box
lat_idx <- which(lat >= lat_min & lat <= lat_max)
lon_idx <- which(lon >= lon_min & lon <= lon_max)

fsle <- ncvar_get(FSLE_nc, "fsle_max", start = c(min(lon_idx), min(lat_idx),1),
                  count = c(length(lon_idx), length(lat_idx),1))
fsle_orient <- ncvar_get(FSLE_nc, "theta_max", start = c(min(lon_idx), min(lat_idx),1),
                         count = c(length(lon_idx), length(lat_idx),1))

fsle_df <- expand.grid(lon = lon[lon_idx], lat = lat[lat_idx], date = date)
fsle_df$fsle <- as.vector(fsle)
fsle_df$fsle_orient <- as.vector(fsle_orient)
fsle_df$date <- as.Date(grid$date, origin = '1950-01-01')


# ------------------------ Step 2: Spatial Plot ---------------------------
# coordinate setup
coord <- st_as_sf(fsle_df, coords = c('lon', 'lat'))
#coord <- st_crs('EPSG:4326')$coord

#Plot it:

ggplot(coord) + 
  #geom_sf(aes(color=fsle)) +
  geom_sf(aes(ifelse(fsle==0, color = 'gray40', color=fsle)))

ggplot(coord) +
  geom_sf(aes(color = fsle_orient))