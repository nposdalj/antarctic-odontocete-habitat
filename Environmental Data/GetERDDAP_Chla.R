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

saveDir = paste("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data",sep="")
AQUA_sstInfo <- info('erdMH1chla1day_Lon0360')


GetChla <- function(site) {
  # Specify Lat and Long and time
  # Doing +/- 0.5 degree for now, can change if something else is better
  if (site == 'EI'){
    # site lat-long: -60.8869, -55.95400
    latitude = c(-59.3869,-61.3869)
    longitude = c(303.406, 304.406)
    years = c(2014)
    time = c(paste(years,"-03-05",sep=""),paste(years,"-07-17",sep="")) }
  if (site == 'KGI'){
    # site lat-long: -61.457817, -57.941917
    latitude = c(-60.957817,-61.957817)
    longitude = c(302.558083, 303.558083)
    years = c(2015,2016)
    time = c(paste(years[1],"-02-10",sep=""),paste(years[2],"-01-29",sep="")) }

  if (site == 'CI'){
    # site lat-long: -61.251867, -53.483433
    latitude = c(-60.751867,-61.751867)
    longitude = c(306.016567, 307.016567)
    years = c(2016)
    time = c(paste(years,"-02-04",sep=""),paste(years,"-12-02",sep="")) }
  year_str <- paste(years, collapse = "-")
  AQUASST <- griddap(AQUA_sstInfo, latitude = latitude, longitude = longitude,
                     time = time, fields = 'chlorophyll',
                     store = disk(path = paste(saveDir, "/",site,"_ChlA_",year_str,".nc",sep =""), overwrite = TRUE))

}
GetChla('EI')
GetChla('KGI')
GetChla('CI')
