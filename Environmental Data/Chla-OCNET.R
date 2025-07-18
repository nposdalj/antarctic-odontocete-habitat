# Load required libraries
library(ncdf4)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(fields)

# ----------- Step 1: Extract Data ------------
EI_chla <- nc_open("I:/Chlorophyll - OCNET/OCNET_chla_2014.nc")
KGI_chla <- nc_open("I:/Chlorophyll - OCNET/OCNET_chla_2015.nc")
CI_chla <- nc_open("I:/Chlorophyll - OCNET/OCNET_chla_2016.nc")

dfFromNC <- function(data,site) {
  if (site == 'EI'){
    # site lat-long: -60.8869, -55.95400
    latitude = c(-59.3869,-61.3869)
    longitude = c(-56.594, -55.594)
    years = c(2014)
    start <- as.Date('20140305',format='%Y%m%d')
    end <- as.Date('20140717',format='%Y%m%d') 
    ## Test you're plotting the right place
    plot(1, type = "n", xlim = c(-60, -50), ylim = c(-65, -55),
         xlab = "Longitude", ylab = "Latitude", main = "Zoomed Bounding Box View")
    maps::map("world", add = TRUE)
    rect(xleft = min(longitude), xright = max(longitude),
         ybottom = min(latitude), ytop = max(latitude),
         border = "red", lwd = 2)}
  if (site == 'KGI'){
    # site lat-long: -61.457817, -57.941917
    latitude = c(-60.957817,-61.957817)
    longitude = c(-57.44917, -56.441917)
    years = c(2015,2016)
    start <- as.Date('20150210',format='%Y%m%d')
    end <- as.Date('20160129',format='%Y%m%d')
    plot(1, type = "n", xlim = c(-60, -50), ylim = c(-65, -55),
         xlab = "Longitude", ylab = "Latitude", main = "Zoomed Bounding Box View")
    maps::map("world", add = TRUE)
    rect(xleft = min(longitude), xright = max(longitude),
         ybottom = min(latitude), ytop = max(latitude),
         border = "red", lwd = 2)}
  if (site == 'CI'){
    # site lat-long: -61.251867, -53.483433
    latitude = c(-60.751867,-61.751867)
    longitude = c(-53.983433, -52.983433)
    years = c(2016)
    start <- as.Date('20160204',format='%Y%m%d')
    end <- as.Date('20161202',format='%Y%m%d') 
    plot(1, type = "n", xlim = c(-60, -50), ylim = c(-65, -55),
         xlab = "Longitude", ylab = "Latitude", main = "Zoomed Bounding Box View")
    maps::map("world", add = TRUE)
    rect(xleft = min(longitude), xright = max(longitude),
         ybottom = min(latitude), ytop = max(latitude),
         border = "red", lwd = 2)}
  
  lat <- ncvar_get(data, 'latitude')
  lat_vec <- lat[,1]  # Latitude is constant along columns
  lat_order <- order(lat_vec)
  lat_vec_sorted <- lat_vec[lat_order]
  lat_min <- min(latitude)
  lat_max <- max(latitude)
  lat_idx_sorted <- which(lat_vec_sorted >= lat_min & lat_vec_sorted <= lat_max)
  lat_idx <- lat_order[lat_idx_sorted]
  
  lon <- ncvar_get(data, 'longitude')
  lon_vec <- lon[1,]  # Longitude is constant along rows
  lon_min <- min(longitude)
  lon_max <- max(longitude)
  lon_idx <- which(lon_vec >= lon_min & lon_vec <= lon_max)
  
  time <- ncvar_get(data, 'time')
  time <- as.Date(time,origin='1900-01-01')
  time_min <- as.Date(start)
  time_max <- as.Date(end)
  time_idx <- which(time >= time_min & time <= time_max)

    if (length(lat_idx) == 0 || length(lon_idx) == 0 || length(time_idx) == 0) {
    stop("latitude, longitude, or time index is empty - fix")
  }
  chla <- ncvar_get(data, 'Chla', start = c(min(lat_idx), min(lon_idx),min(time_idx)),
                    count = c(length(lat_idx), length(lon_idx),length(time_idx)))
  
  dates <- as.Date(seq(start, end, by = "day"))
  df <- data.frame()
  
  grid <- expand.grid(lon = lon[lon_idx], lat = lat[lat_idx], time=time[time_idx])
  grid$chla <- as.vector(chla)
  grid$date <- as.Date(current)
  
  df <- grid %>%
    group_by(time) %>%
    summarise(chla = mean(chla, na.rm = TRUE), .groups = "drop") %>%
    rename(date = time)
  
  return(df)
}
EI_df <- dfFromNC(EI_chla,'EI')
KGI_df <- dfFromNC(KGI_chla,'KGI')
CI_df <- dfFromNC(CI_chla,'CI')

# ------------ Step 2: Make Timeseries -----------
# Compare to remotely sensed data

# ------------ Step 3: Save Dataframe ------------
EI_df$Site <- 'EI'
KGI_df$Site <- 'KGI'
CI_df$Site <- 'CI'

final <- rbind(EI_df,KGI_df,CI_df)
write.csv(final, "FILEPATH")