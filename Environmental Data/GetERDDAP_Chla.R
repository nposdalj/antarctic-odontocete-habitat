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

# ------------------Step 2: Extract Data-----------------
# Load ERDDAP chlorophyll-a data
EI_Chla <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/EI_ChlA_2014.nc/740f1cfbd25414d9582871d5d75f46e7.nc")
KGI_Chla <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/KGI_ChlA_2015-2016.nc/3eea4e994f89226515fd0f0dc7b7e210.nc")
CI_Chla <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/CI_ChlA_2016.nc/b18777cf8092a4ed9e97382edc7ba6b8.nc")

dfFromNC <- function(data) {
  # Extracting relevant data
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  chla <- ncvar_get(data, 'chlorophyll')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  
  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   chla = as.vector(chla))
  df <- na.omit(df)
  
  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(chla = mean(chla))
  return(avg_df)
}
# Constructing site dataframes
EI_df <- dfFromNC(EI_Chla)
KGI_df <- dfFromNC(KGI_Chla)
CI_df <- dfFromNC(CI_Chla)

# Function to create timeseries by site
ChlTimeseries <- function(data, site) {
  data$date <- as.Date(data$date)
  
 windows(16,4)
  chla <- ggplot(data = data, mapping = aes(x = date, y = chla)) + geom_line(color = "yellowgreen", size = 2) + 
    geom_point(color = "darkolivegreen", size = 1) +
    labs(x = "Date", y = "Chlorophyll-a (mg/m^3)", title = site) + scale_x_date(date_labels = "%b %Y") +
    theme(axis.title = element_text(size = 13))
  return(chla)
}
# Construct timeseries for all sites
ChlTimeseries(EI_df, "Elephant Island")
ChlTimeseries(KGI_df, "King George Island")
ChlTimeseries(CI_df, "Clarence Island")