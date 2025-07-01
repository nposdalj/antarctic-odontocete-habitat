#This script queries the ERDDAP database for data
#Learn more about the database here: https://coastwatch.pfeg.noaa.gov/erddap/index.html
#TS 06/30/2025
#Uses data from this link - https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1sstd1dayR20190SQ_Lon0360.html

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
# # Regions = c('ENP','CCE','CentralPac')
# Regions = c('CentralPac')
AQUA_sstInfo <- info('erdMH1sstd1dayR20190SQ_Lon0360')


GetSST <- function(site) {
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
                     time = time, fields = 'sstMasked',
                     store = disk(path = paste(saveDir, "/",site,"_SST_",year_str,".nc",sep =""), overwrite = TRUE))
  
}
GetSST('EI')
GetSST('KGI')
GetSST('CI')

# Natalie code:
# for (i in 1:length(Regions)){
#   Region = Regions[i]
# 
#   #Specify Lat and Long and time
#   if (Region == 'CCE'){
#     latitude = c(48.,28.)
#     longitude = c(230., 247)
#     years = c(2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022)}
#   
#   if (Region == 'ENP'){
#     latitude = c(60.,50.)
#     longitude = c(170., 220.)
#     years = c(2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)}
#   
#   if (Region == 'CentralPac'){
#     latitude = c(30.,0.)
#     longitude = c(145., 216)
#     years = c(2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)}
#   
# for (yy in 1:length(years)){
#   year = as.character(years[yy])
#   if (Region == 'CCE' & year == "2022"){
#     time = c(paste(year,"-01-01",sep=""),paste(year,"-06-29",sep=""))  
#   }else{
#     time = c(paste(year,"-01-01",sep=""),paste(year,"-12-31",sep=""))}
#   AQUASST <- griddap(AQUA_sstInfo, latitude = latitude, longitude = longitude, 
#                      time = time, fields = 'sstMasked',
#                      store = disk(path = paste(saveDir,Region,"_SST_",year,".nc",sep =""), overwrite = TRUE))
#   
# }
# }