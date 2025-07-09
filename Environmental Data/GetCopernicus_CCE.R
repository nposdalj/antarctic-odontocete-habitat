library(lubridate)
library(CopernicusMarine)
library(stars)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(tidyverse)
library(gridExtra) # for grid.arrange

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

# -----------------Step 2: Extract Data-------------------
# Load Copernicus data
# surface data, all vars 
EI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/EI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990714294.nc")
KGI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/KGI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990804735.nc")
CI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/CI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990301834.nc")

dfFromNC <- function(data) {
  # Extracting relevant data
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  ssh <- ncvar_get(data, 'zos')
  n_velocity <- ncvar_get(data, 'vo')
  e_velocity <- ncvar_get(data, 'uo')
  salinity <- ncvar_get(data, 'so')
  mixed_layer <- ncvar_get(data, 'mlotst')
  # Ice variables
  sice_conc <- ncvar_get(data,'siconc')
  sice_thick <- ncvar_get(data, 'sithick')
  sice_e_veloc <- ncvar_get(data, 'usi')
  sice_n_veloc <- ncvar_get(data, 'vsi')

  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   ssh = as.vector(ssh), n_velocity = as.vector(n_velocity),
                   e_velocity = as.vector(e_velocity), salinity = as.vector(salinity),
                   mixed_layer = as.vector(mixed_layer), sice_conc = as.vector(sice_conc), sice_thick = as.vector(sice_thick),
                   sice_e_veloc = as.vector(sice_e_veloc), sice_n_veloc = as.vector(sice_n_veloc))

  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(ssh = mean(ssh), n_velocity = mean(n_velocity),
                                                e_velocity = mean(e_velocity), salinity = mean(salinity),
                                                mixed_layer = mean(mixed_layer), sice_conc = mean(sice_conc),
                                                sice_thick = mean(sice_thick), sice_e_veloc = mean(sice_e_veloc),
                                                sice_n_veloc = mean(sice_n_veloc))
  return(avg_df)
}
# Constructing site dataframes
EI_df <- dfFromNC(EI_cop)
KGI_df <- dfFromNC(KGI_cop)
CI_df <- dfFromNC(CI_cop)

# --------------------------Step 3: Compute EKE--------------------------
getEKE <- function(data) {
  # get u' and v' for EKE equation
  u_mean <- mean(data$e_velocity) 
  v_mean <- mean(data$n_velocity)
  data$u_prime <- data$e_velocity - u_mean
  data$v_prime <- data$n_velocity - v_mean
  
  # EKE = 1/2((u')^2 + (v')^2)
  # units: cm^2/s^2
  data$EKE <- (0.5 * ((data$u_prime^2) + (data$v_prime^2))) * 10000
  return(data)
}
EI_df <- getEKE(EI_df)
KGI_df <- getEKE(KGI_df)
CI_df <- getEKE(CI_df)

# -------------------------Step 4: Make Oceanographic Data Timeseries------------------------
# Function to create timeseries by site
CopTimeseries <- function(data, site) {
  data$date <- as.Date(data$date)
  # sea surface height
  ssh <- ggplot(data = data, mapping = aes(x = date, y = ssh)) + geom_line(color = "tomato", size = 1) + 
    labs(x = "Sea Surface Height", y = "meters") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # eddy kinetic energy
    EKE <- ggplot(data = data, mapping = aes(x = date, y = EKE)) + 
    geom_line(color = "mediumpurple", size = 1) + 
    labs(x = "Eddy Kinetic Energy", y = "cm^2/s^2") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # salinity
  salinity <- ggplot(data = data, mapping = aes(x = date, y = salinity)) + geom_line(color = "gold", size = 1) + 
    labs(x = "Salinity", y = "psu") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # mixed layer depth
  mixed <- ggplot(data = data, mapping = aes(x = date, y = mixed_layer)) + 
    geom_line(color = "mediumseagreen", size = 1) + 
    labs(x = "Mixed Layer Thickness", y = "meters") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  # Arrange plots to one figure
  windows()
  final_plot <- grid.arrange(ssh, EKE, salinity, mixed, nrow = 4, top = site)
  
}
# Construct timeseries for all sites
CopTimeseries(EI_df, "Elephant Island")
CopTimeseries(KGI_df, "King George Island")
CopTimeseries(CI_df, "Clarence Island")

# -------------------------Step 4: Make Sea Ice Data Timeseries------------------------
# Function to create timeseries by site
IceTimeseries <- function(data, site) {
  data$date <- as.Date(data$date)
  # sea ice thickness
  thickness <- ggplot(data = data, mapping = aes(x = date, y = sice_thick)) + geom_line(color = "tomato", size = 1) + 
    labs(x = "Sea Ice Thickness", y = "meters") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # 
  EKE <- ggplot(data = data, mapping = aes(x = date, y = EKE)) + 
    geom_line(color = "mediumpurple", size = 1) + 
    labs(x = "Eddy Kinetic Energy", y = "cm^2/s^2") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # salinity
  salinity <- ggplot(data = data, mapping = aes(x = date, y = salinity)) + geom_line(color = "gold", size = 1) + 
    labs(x = "Salinity", y = "psu") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # mixed layer depth
  mixed <- ggplot(data = data, mapping = aes(x = date, y = mixed_layer)) + 
    geom_line(color = "mediumseagreen", size = 1) + 
    labs(x = "Mixed Layer Thickness", y = "meters") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  # Arrange plots to one figure
  windows()
  final_plot <- grid.arrange(ssh, EKE, salinity, mixed, nrow = 4, top = site)
  
}
# Construct timeseries for all sites
IceTimeseries(EI_df, "Elephant Island")
IceTimeseries(KGI_df, "King George Island")
IceTimeseries(CI_df, "Clarence Island")