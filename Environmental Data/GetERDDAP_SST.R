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

# --------------Step 1: Access ERRDAP data-------------------
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

# ------------------Step 2: Extract Data-----------------
# Load ERDDAP SST data
EI_sst <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/EI_SST_2014.nc/683aefaef697b5b3566df97963a5795b.nc")
KGI_sst <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/KGI_SST_2015-2016.nc/d845de2ea7f37f907cc6bbcfc07dd9a3.nc")
CI_sst <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/CI_SST_2016.nc/eab0cca18b09862e5e68457478cf25d0.nc")

dfFromNC <- function(data) {
  # Extracting relevant data
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  sst <- ncvar_get(data, 'sstMasked')
  # fillvalue <- ncatt_get(data, "sea_surface_temperature", "_FillValue")
  # sst[sst==fillvalue$value] <- NA
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  
  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   sst = as.vector(sst))
  df <- na.omit(df) # removing missing data
  
  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(sst = mean(sst))
  return(avg_df)
}
# Constructing site dataframes
EI_df <- dfFromNC(EI_sst)
KGI_df <- dfFromNC(KGI_sst)
CI_df <- dfFromNC(CI_sst)

# Function to create timeseries by site
SSTtimeseries <- function(data, site) {
  data$date <- as.Date(data$date)
  
  windows(16,4)
  sst <- ggplot(data = data, mapping = aes(x = date, y = sst)) + geom_line(color = "tan", size = 2) + 
    geom_point(color = "saddlebrown", size = 1) + 
    labs(x = "Date", y = "Sea Surface Temperature (C)", title = site) + scale_x_date(date_labels = "%b %Y") +
    theme(axis.title = element_text(size = 13))
  return(sst)
}
# Construct timeseries for all sites
SSTtimeseries(EI_df, "Elephant Island")
SSTtimeseries(KGI_df, "King George Island")
SSTtimeseries(CI_df, "Clarence Island")