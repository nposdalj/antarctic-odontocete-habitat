library(lubridate)
library(CopernicusMarine)
library(stars)
library(ncdf4) # package for netcdf manipulation
#library(RNetCDF)
library(raster) # package for raster manipulation
library(tidyverse)

# -----------Step 1: Access Copernicus data-----------
# Was getting 404 error when trying to download data (even though I changed function to match the updated CopernicusMarine library)
# Downloaded Copernicus data directly from the website
# Variables (daily resolution): east water velocity, north water velocity, SSH, mixed layer thickness, salinity

# YEARS = c(2014, 2015, 2016)
# 
# for (i in 1:length(YEARS)){
#   year = YEARS[i]
#   out_dir = paste("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data",as.character(year),".nc",sep = "")
#   timerange1 = paste(as.character(year),"-01-01",sep="")
#   timerange2 = paste(as.character(year),"-12-31",sep="")
#   cms_download_subset(
#     username = getOption("CopernicusMarine_uid", "tsharan"),
#     password = getOption("CopernicusMarine_pwd", "H!ghfr3qu3ncy"),
#     destination = out_dir,
#     product = "GLOBAL_MULTIYEAR_PHY_001_030", # changed from GLOBAL_REANALYSIS_PHY_001_031
#     layer = "cmems_mod_glo_phy_my_0.083deg_P1D-m_202311/zos",
#     variable = "zos",
#     # output = "netcdf",
#     region = c(-59, -62, -52, -60), # bounding box 52 t0 59 W & 60 to 62 S
#     timerange = c(timerange1,timerange2),
#   # sub_variables = c("mlotst_mean","so_mean","thetao_mean","uo_mean",
#   # "vo_mean","zos_mean"),
#   # sub_variables = c("mlotst_mean","mlotst_std","so_mean","so_std","thetao_mean","thetao_std","uo_mean",
#   #                   "uo_std","vo_mean","vo_std","zos_mean","zos_std"),
# #   #verticalrange = c(0.5057600140571594,1.56)
#   #verticalrange = c(0,-2)
# )
# }

# -----------------Step 2: Construct Timeseries-------------------
# Load Copernicus data
EI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/EI_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751320063762.nc")
KGI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/KGI_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751320609977.nc")
CI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/CI_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751320549137.nc")

dfFromNC <- function(data) {
  # Extracting relevant data
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  ssh <- as.vector(ncvar_get(data, 'zos'))
  n_velocity <- as.vector(ncvar_get(data, 'vo'))
  e_velocity <- as.vector(ncvar_get(data, 'uo'))
  salinity <- as.vector(ncvar_get(data, 'so'))
  mixed_layer <- as.vector(ncvar_get(data, 'mlotst'))
  
  # Cleaning up latitude, longitude, and time variables
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  lonlattime <- as.matrix(expand.grid(lon,lat,time_obs))
  
  # Creating a dataframe with all variables
  lonlat <- expand.grid(lon = lon, lat = lat, time = time_obs)
  
  df <- data.frame(
    long = lonlat$lon,
    lat = lonlat$lat,
    date = lonlat$time,
    ssh = as.vector(ssh),
    n_velocity = as.vector(n_velocity),
    e_velocity = as.vector(e_velocity),
    salinity = as.vector(salinity),
    mixed_layer = as.vector(mixed_layer)
  )
  # df <- data.frame(cbind(lonlattime, ssh, n_velocity, e_velocity, salinity, mixed_layer))
  # colnames(df) <- c("long","lat","date","ssh", "n_velocity",
  #                         "e_velocity", "salinity", "mixed_layer")

  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(ssh = mean(ssh), n_velocity = mean(n_velocity),
                                                e_velocity = mean(e_velocity), salinity = mean(salinity),
                                                mixed_layer = mean(mixed_layer))

  return(avg_df)
  }
dfFromNC(EI_cop)
KGI_df <- dfFromNC(KGI_cop)
CI_df <- dfFromNC(CI_cop)

CopTimeseries <- function(data) {
  
}