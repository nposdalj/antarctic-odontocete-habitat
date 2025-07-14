library(rerddap)
library(rerddapXtracto)
library(ncdf4)
library(parsedate)
library(sp)
library(gganimate)
library(ggplot2)
library(plotdap)
library(abind)
library(dplyr)

# Only relevant AVISO variables are FSLE
# Ignore functions for other variables (data already obtained from Copernicus)

# --------------------Step 1: Subset and Save Relevant AVISO data--------------------------
subsetAVISO <- function(site) {
  directory <- "D:/FSLE/AVISO_Global_FSLE/dt_global_allsat_madt_fsle_"
  
  if (site == 'EI'){
    # site lat-long: -60.8869, -55.95400
    latitude = c(-59.3869,-61.3869)
    longitude = c(303.406, 304.406)
    years = c(2014)
    start <- as.Date('20140305',format='%Y%m%d')
    end <- as.Date('20140717',format='%Y%m%d') }
  if (site == 'KGI'){
    # site lat-long: -61.457817, -57.941917
    latitude = c(-60.957817,-61.957817)
    longitude = c(302.558083, 303.558083)
    years = c(2015,2016)
    start <- as.Date('20150210',format='%Y%m%d')
    end <- as.Date('20160129',format='%Y%m%d') }
  if (site == 'CI'){
    # site lat-long: -61.251867, -53.483433
    latitude = c(-60.751867,-61.751867)
    longitude = c(306.016567, 307.016567)
    years = c(2016)
    start <- as.Date('20160204',format='%Y%m%d')
    end <- as.Date('20161202',format='%Y%m%d') }
  
  dates <- seq(start, end, by = "day")
  df <- data.frame()
  for(current in dates) {
    date_str <- format(as.Date(current), "%Y%m%d")
    file_path <- paste(directory,date_str,'_20180704.nc',sep="")
    
    if (!file.exists(file_path)) {
      warning(paste("File does not exist:", file_path))
      next
    }
    nc <- nc_open(file_path)
    
    # Extract lat long and constrain to bounding box
    lon <- ncvar_get(nc, 'lon')
    lat <- ncvar_get(nc,'lat')
    time <- ncvar_get(nc, 'time')
    # Define bounding box
    lat_min <- latitude[2]
    lat_max <- latitude[1]
    lon_min <- longitude[1]
    lon_max <- longitude[2]
    # Find index ranges for the bounding box
    lat_idx <- which(lat >= lat_min & lat <= lat_max)
    lon_idx <- which(lon >= lon_min & lon <= lon_max)
    
    fsle <- ncvar_get(nc, "fsle_max", start = c(min(lon_idx), min(lat_idx),1),
                     count = c(length(lon_idx), length(lat_idx),1))
    
    grid <- expand.grid(lon = lon[lon_idx], lat = lat[lat_idx], time = time)
    grid$temp <- as.vector(fsle)
    grid$date <- as.Date(current)
    
    # Average over space for each depth
    daily_avg <- grid %>%
      group_by(date) %>%
      summarise(fsle = mean(fsle, na.rm = TRUE),
                .groups = "drop")
    
    df <- bind_rows(df, daily_avg)
    print(paste("File for ", current, "done."))
  }
  return(df)
}
EI_fsle <- subsetAVISO('EI')
write.csv(EI_fsle, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/AVISO/EI_fsle")
KGI_fsle <- subsetAVISO('KGI')
write.csv(KGI_fsle, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/AVISO/KGI_fsle")
CI_fsle <- subsetAVISO('CI')
write.csv(CI_fsle, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/AVISO/CI_fsle")

# ------------------ Step 2: FSLE Timeseries -------------------------
EI <- ggplot(data = EI_fsle, mapping = aes(x = date, y = fsle)) + 
  geom_line(color = "mediumvioletred", linewidth = 0.7) + 
  labs(x = NULL, y = "FSLE", title="Elephant Island") +
  scale_x_date(date_labels = "%b %Y")
KGI <- ggplot(data = KGI_fsle, mapping = aes(x = date, y = fsle)) + 
  geom_line(color = "mediumvioletred", linewidth = 0.7) + 
  labs(x = NULL, y = "FSLE", title="King George Island") +
  scale_x_date(date_labels = "%b %Y")
CI <- ggplot(data = CI_fsle, mapping = aes(x = date, y = fsle)) + 
  geom_line(color = "mediumvioletred", linewidth = 0.7) + 
  labs(x = NULL, y = "FSLE", title="Clarence Island") +
  scale_x_date(date_labels = "%b %Y")


