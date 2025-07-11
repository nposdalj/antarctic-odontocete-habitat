#This script queries the ERDDAP database for data
#Learn more about the database here: https://coastwatch.pfeg.noaa.gov/erddap/index.html
#TS 06/30/2025
#Uses data from this link - https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chla1day_Lon0360.html
#Load packages
library(rerddap)
library(rerddapXtracto)
library(ncdf4)
library(parsedate)
library(sp)
library(gganimate)
library(ggplot2)
library(plotdap)

# ---------------Step 1: Access ERDDAP data----------------
saveDir = paste("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Salinity_SMOS/",sep="")
SMOS_info <- info('coastwatchSMOSv662SSS1day', url="https://coastwatch.pfeg.noaa.gov/erddap/")


GetSSS <- function(site) {
  # Specify Lat and Long and time
  # Doing +/- 0.5 degree for now, can change if something else is better
  if (site == 'EI'){
    # site lat-long: -60.8869, -55.95400
    latitude = c(-59.3869,-61.3869)
    longitude = c(-56.594, -55.594)
    years = c(2014)
    time = c(paste(years,"-03-05",sep=""),paste(years,"-07-17",sep="")) }
  if (site == 'KGI'){
    # site lat-long: -61.457817, -57.941917
    latitude = c(-60.957817,-61.957817)
    longitude = c(-57.441917, -56.441917)
    years = c(2015,2016)
    time = c(paste(years[1],"-02-10",sep=""),paste(years[2],"-01-29",sep="")) }
  
  if (site == 'CI'){
    # site lat-long: -61.251867, -53.483433
    latitude = c(-60.751867,-61.751867)
    longitude = c(-53.983433, -52.983433)
    years = c(2016)
    time = c(paste(years,"-02-04",sep=""),paste(years,"-12-02",sep="")) }
  year_str <- paste(years, collapse = "-")
  SMOS <- griddap(SMOS_info, latitude = latitude, longitude = longitude,
                     time = time, fields = 'sss',
                     store = disk(path = paste(saveDir,site,"_SSS_.nc",sep =""), overwrite = TRUE))
  
}
GetSSS('EI')
GetSSS('KGI')
GetSSS('CI')

# ------------------Step 2: Extract Data-----------------
# Load ERDDAP SMOS salinity data
EI_SSS <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Salinity_SMOS/EI_SSS_.nc/c0bd29c5fba6884e6028ecfb2164b8d1.nc")
KGI_SSS <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Salinity_SMOS/KGI_SSS_.nc/noaacwSMOSsssDaily_d213_6ee2_4c0a_U1752249399382.nc")
CI_SSS <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Salinity_SMOS/CI_SSS_.nc/9a0327e1abdf2a14c0907b1ed6ff19b3.nc")

dfFromNC <- function(data) {
  # Extracting relevant data
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  salinity <- ncvar_get(data, 'sss')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  
  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   salinity = as.vector(salinity))
  df <- na.omit(df)
  
  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(salinity = mean(salinity))
  return(avg_df)
}
# Constructing site dataframes
EI_df <- dfFromNC(EI_SSS)
KGI_df <- dfFromNC(KGI_SSS)
CI_df <- dfFromNC(CI_SSS)
write.csv(EI_df, paste(saveDir,"EI_salinity.csv"))
write.csv(KGI_df, paste(saveDir,"KGI_salinity.csv"))
write.csv(CI_df, paste(saveDir,"CI_salinity.csv"))
