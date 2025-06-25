#This script queries the ERDDAP database for data
#Learn more about the database here: https://coastwatch.pfeg.noaa.gov/erddap/index.html
#NP 07/09/2021
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

saveDir = paste("L:/My Drive/FourLoko/HabitatVariables/",sep="")
# Regions = c('ENP','CCE','CentralPac')
Regions = c('CentralPac')
AQUA_sstInfo <- info('erdMH1chla1day_Lon0360')

for (i in 1:length(Regions)){
  Region = Regions[i]

  #Specify Lat and Long and time
  if (Region == 'CCE'){
    latitude = c(48.,28.)
    longitude = c(230., 247)
    years = c(2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022)}
  
  if (Region == 'ENP'){
    latitude = c(60.,50.)
    longitude = c(170., 220.)
    years = c(2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)}
  
  if (Region == 'CentralPac'){
    latitude = c(30.,0.)
    longitude = c(145., 216)
    years = c(2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)}
  
for (yy in 1:length(years)){
  year = as.character(years[yy])
  if (Region == 'CCE' & year == "2022"){
    time = c(paste(year,"-01-01",sep=""),paste(year,"-07-27",sep=""))  
  }else{
    time = c(paste(year,"-01-01",sep=""),paste(year,"-12-31",sep=""))}
  AQUASST <- griddap(AQUA_sstInfo, latitude = latitude, longitude = longitude, 
                     time = time, fields = 'chlorophyll',
                     store = disk(path = paste(saveDir,Region,"_ChlA_",year,".nc",sep =""), overwrite = TRUE))
  
}
}