library(CopernicusMarine)
library(stars)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(tidyverse)
library(gridExtra) # for grid.arrange

# -----------------Step 2: Extract Oceanographic Data-------------------
# Load Copernicus data
# surface data, all vars
# EI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/EI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990714294.nc")
# KGI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/KGI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990804735.nc")
# CI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/CI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990301834.nc")
# 0.5 degree, all vars
# EI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/EIhfdeg_cmems_mod_glo_phy_my_0.083deg_P1D-m_1752518931575.nc")
# KGI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/KGIhfdeg_cmems_mod_glo_phy_my_0.083deg_P1D-m_1752529248292.nc")
# CI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/CIhfdeg_cmems_mod_glo_phy_my_0.083deg_P1D-m_1752529427864.nc")

# all physical vars 
EI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/EI40km_cmems_mod_glo_phy_my_0.083deg_P1D-m_1753134290538.nc")
KGI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/KGI40km_cmems_mod_glo_phy_my_0.083deg_P1D-m_1753134447092.nc")
CI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/CI40km_cmems_mod_glo_phy_my_0.083deg_P1D-m_1753134569824.nc")

# biogeochem vars (oxygen and chlorophyll)
EI_bio <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/EI40km_cmems_mod_glo_bgc_my_0.25deg_P1D-m_1753133719978.nc")
KGI_bio <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/KGI40km_cmems_mod_glo_bgc_my_0.25deg_P1D-m_1753133880987.nc")
CI_bio <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/CI40km_cmems_mod_glo_bgc_my_0.25deg_P1D-m_1753134029874.nc")

# Ensuring that variable dimensions index correctly
fix_grid <- function(var, four_dim = FALSE) {
  if (four_dim) {
    # From [lon, lat, depth, time] → [lon, lat, time, depth]
    var <- aperm(var, c(1, 2, 4, 3))  # moves depth to last
  } else {
    # 3D case [lon, lat, time] → [lon, lat, time], no change
    var <- aperm(var, c(1, 2, 3))
  }
  return(as.vector(var))
}

physFromNC <- function(data) {
  # Extracting relevant data
  # Dimensions
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  depth <- ncvar_get(data, 'depth')
  # Depths (based off dive depths of detected mammals)
  # surface: 0.494025 
  # Long-finned pilot: 15.810070, 643.566772
  # Killer: 11.405000, 453.937714 
  # Sperm: 380.213013, 1684.284058
  # Gray's and Strap-toothed: 65.807266, 902.339294
  # Southern bottlenose: 763.333130
  
  # Oceanographic variables
  temp <- fix_grid(ncvar_get(data, 'thetao'), four_dim=TRUE)
  ssh <- fix_grid(ncvar_get(data, 'zos'))
  n_velocity <- fix_grid(ncvar_get(data, 'vo'), four_dim=TRUE)
  e_velocity <- fix_grid(ncvar_get(data, 'uo'), four_dim=TRUE)
  salinity <- fix_grid(ncvar_get(data, 'so'), four_dim=TRUE)
  mixed_layer <- fix_grid(ncvar_get(data, 'mlotst'))
  # Ice variables
  sice_conc <- fix_grid(ncvar_get(data,'siconc'))
  sice_thick <- fix_grid(ncvar_get(data, 'sithick'))
  sice_e_veloc <- fix_grid(ncvar_get(data, 'usi'))
  sice_n_veloc <- fix_grid(ncvar_get(data, 'vsi'))

  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs, depth = depth)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   depth = lonlattime$depth,
                   ssh = as.vector(ssh), n_velocity = as.vector(n_velocity),
                   e_velocity = as.vector(e_velocity), salinity = as.vector(salinity),
                   mixed_layer = as.vector(mixed_layer), sice_conc = as.vector(sice_conc), sice_thick = as.vector(sice_thick),
                   sice_e_veloc = as.vector(sice_e_veloc), sice_n_veloc = as.vector(sice_n_veloc), temp=as.vector(temp))
  
  # Filtering by different depths & correcting for rounding errors
  df <- df %>% mutate(depth = round(depth, 3))
  filtered <- df %>% filter(depth %in% round(c(0.494025, 15.810070, 453.937714, 1684.284058, 380.213013, 643.566772, 11.405000,
                      65.807266, 763.333130, 902.339294),3))
  filtered <- filtered %>% mutate(date = as.Date(date))

  # Creating a dataframe with daily spatial averages
  # spatial standard deviations as well
  # Mean absolute deviation for north & east velocity to make variability calculations easier for EKE later on
  avg_df <- filtered %>% group_by(date, depth) %>% summarize(ssh_mean = mean(ssh, na.rm=TRUE), ssh_sd = sd(ssh, na.rm=TRUE),
                                                             n_velocity_mean = mean(n_velocity, na.rm=TRUE), n_velocity_mad = mad(n_velocity, na.rm=TRUE),
                                                             e_velocity_mean = mean(e_velocity, na.rm=TRUE), e_velocity_mad = mad(e_velocity,na.rm=TRUE),
                                                             salinity_mean = mean(salinity, na.rm=TRUE), salinity_sd = sd(salinity,na.rm=TRUE),
                                                             mixed_layer_mean = mean(mixed_layer, na.rm=TRUE), mixed_layer_sd = sd(mixed_layer, na.rm=TRUE),
                                                             sice_conc_mean = mean(sice_conc, na.rm=TRUE), ice_conc_sd = sd(sice_conc, na.rm=TRUE),
                                                             sice_thick_mean = mean(sice_thick, na.rm=TRUE), sice_e_veloc_mean = mean(sice_e_veloc, na.rm=TRUE),
                                                             sice_n_veloc_mean = mean(sice_n_veloc, na.rm=TRUE), 
                                                             temp_mean=mean(temp, na.rm=TRUE), temp_sd=sd(temp, na.rm=TRUE))
  avg_df <- ungroup(avg_df)
  return(avg_df)
}
# Constructing site dataframes
EI_phys <- physFromNC(EI_cop)
KGI_phys <- physFromNC(KGI_cop)
CI_phys <- physFromNC(CI_cop)

# Similar function to above, for biogeochem vars
bioFromNC <- function(data) {
  # Extracting relevant data
  # Dimensions
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  depth <- ncvar_get(data, 'depth')
  # Depths (based off dive depths of detected mammals)
  # surface: 0.506 
  # Long-finned pilot: 16.525, 628.026
  # Killer: 11.774, 457.626 
  # Sperm: 370.689, 1652.568
  # Gray's and Strap-toothed: 69.022, 947.448
  # Southern bottlenose: 773.368
  
  # biogeochemical variables
  chla <- fix_grid(ncvar_get(data, 'chl'), four_dim=TRUE)
  o2 <- fix_grid(ncvar_get(data, 'o2'), four_dim=TRUE)
  productivity <- fix_grid(ncvar_get(data, 'nppv'), four_dim=TRUE)

  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs, depth = depth)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   depth = lonlattime$depth,
                   chla = as.vector(chla), o2 = as.vector(o2), productivity = as.vector(productivity))
  
  # Filtering by different depths & correcting for rounding errors
  df <- df %>% mutate(depth = round(depth, 3))
  filtered <- df %>% filter(depth %in% round(c(0.506, 16.525, 457.626, 1652.568, 370.689, 628.026, 11.774,
                                               69.022, 773.368, 947.448),3))
  
  # Creating a dataframe with daily spatial averages
  avg_df <- filtered %>% group_by(date, depth) %>% summarize(chla_mean = mean(chla, na.rm=TRUE), chla_sd = sd(chla,na.rm=TRUE),
                                                             o2_mean = mean(o2, na.rm=TRUE), o2_sd = sd(o2, na.rm=TRUE),
                                                             productivity_mean = mean(productivity, na.rm=TRUE), productivity_sd = sd(productivity,na.rm=TRUE))
  avg_df <- ungroup(avg_df)
  return(avg_df)
}
# Constructing site dataframes
EI_biodf <- bioFromNC(EI_bio)
KGI_biodf <- bioFromNC(KGI_bio)
CI_biodf <- bioFromNC(CI_bio)

# combine biogeochemical and physical oceanographic data
# first, match depth classes
depthMatch <- function(data) {
  data$depth[data$depth %in% c(0.506, 0.494)] <- 0.5
  data$depth[data$depth %in% c(11.774, 11.405)] <- 11
  data$depth[data$depth %in% c(15.810, 16.525)] <- 16
  data$depth[data$depth %in% c(65.807, 69.022)] <- 67
  data$depth[data$depth %in% c(370.689, 380.213)] <- 375
  data$depth[data$depth %in% c(457.626, 453.938)] <- 455
  data$depth[data$depth %in% c(628.026, 643.567)] <- 635
  data$depth[data$depth %in% c(773.368, 763.333)] <- 768
  data$depth[data$depth %in% c(947.448, 902.339)] <- 920
  data$depth[data$depth %in% c(1652.568, 1684.284)] <- 1665
  return(data)
}
EI_final <- left_join(depthMatch(EI_biodf),depthMatch(EI_phys), by = c('date', 'depth'))
KGI_final <- left_join(depthMatch(KGI_biodf),depthMatch(KGI_phys), by = c('date', 'depth'))
CI_final <- left_join(depthMatch(CI_biodf),depthMatch(CI_phys), by = c('date', 'depth'))


# --------------------------Step 3: Compute EKE--------------------------
getEKE <- function(data) {
  # get u' and v' for EKE equation
  data <- data %>% group_by(depth) %>% mutate(u_mean = mean(e_velocity_mean,na.rm=TRUE),
                                                   v_mean = mean(n_velocity_mean, na.rm=TRUE))
  data <- data %>% group_by(depth) %>% mutate(u_prime = e_velocity_mean-u_mean, 
                                                 v_prime = n_velocity_mean-v_mean)
  
  # EKE = 1/2((u')^2 + (v')^2)
  # units: cm^2/s^2
  data$EKE_mean <- (0.5 * ((data$u_prime^2) + (data$v_prime^2))) * 10000
  data <- subset(data, select = -c(u_mean, v_mean, u_prime,v_prime))
  data <- data %>% ungroup()
  
  # similar process, but applied on velocity median absolute deviations to get daily MAD for EKE
  # get u' and v' for EKE equation
  data <- data %>% group_by(depth) %>% mutate(u_mean = mean(e_velocity_mad,na.rm=TRUE),
                                              v_mean = mean(n_velocity_mad, na.rm=TRUE))
  data <- data %>% group_by(depth) %>% mutate(u_prime = e_velocity_mad-u_mean, 
                                              v_prime = n_velocity_mad-v_mean)
  
  # EKE = 1/2((u')^2 + (v')^2)
  # units: cm^2/s^2
  data$EKE_mad <- (0.5 * ((data$u_prime^2) + (data$v_prime^2))) * 10000
  data <- subset(data, select = -c(u_mean, v_mean, u_prime,v_prime))
  data <- data %>% ungroup()
  return(data)
}
EI_final <- getEKE(EI_final)
KGI_final <- getEKE(KGI_final)
CI_final <- getEKE(CI_final)

# -------------------------Step 4: Make Oceanographic Data Timeseries------------------------
# Function to create timeseries by site (for surface only)
CopTimeseries <- function(data, site) {
  data$date <- as.Date(data$date)
  data <- filter(data, depth == 0.5)
  # sea surface height
  ssh <- ggplot(data = data, mapping = aes(x = date, y = ssh_mean)) + geom_line(color = "tomato", linewidth = 1) + 
    labs(x = "Sea Surface Height", y = "meters") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # eddy kinetic energy
    EKE <- ggplot(data = data, mapping = aes(x = date, y = EKE_mean)) + 
    geom_line(color = "mediumpurple", linewidth = 1) + 
    labs(x = "Eddy Kinetic Energy", y = "cm^2/s^2") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # salinity
  salinity <- ggplot(data = data, mapping = aes(x = date, y = salinity_mean)) + geom_line(color = "gold", linewidth = 1) + 
    labs(x = "Salinity", y = "psu") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # mixed layer depth
  mixed <- ggplot(data = data, mapping = aes(x = date, y = mixed_layer_mean)) + 
    geom_line(color = "mediumseagreen", linewidth = 1) + 
    labs(x = "Mixed Layer Thickness", y = "meters") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # temperature
  temp <- ggplot(data = data, mapping = aes(x = date, y = temp_mean)) + 
    geom_line(color = "darkred", linewidth = 1) + 
    labs(x = "Temperature", y = "C") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # chlorophyll
  chl <- ggplot(data = data, mapping = aes(x = date, y = chla_mean)) + 
    geom_line(color = "darkgreen", linewidth = 1) + 
    labs(x = "Chlorophyll-a", y = "mg/m3") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # oxygen
  o2 <- ggplot(data = data, mapping = aes(x = date, y = o2_mean)) + 
    geom_line(color = "navy", linewidth = 1) + 
    labs(x = "Oxygen Concentration", y = "mmol/m3") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # productivity
  prod <- ggplot(data = data, mapping = aes(x = date, y = productivity_mean)) + 
    geom_line(color = "deeppink", linewidth = 1) + 
    labs(x = "Primary Production", y = "mg/m3/day") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  
  # Arrange plots to one figure
  windows()
  final_plot <- grid.arrange(ssh, EKE, salinity, mixed, temp, chl, o2, prod, nrow = 8, top = site)
  
}
# Construct timeseries for all sites
CopTimeseries(EI_final, "Elephant Island")
CopTimeseries(KGI_final, "King George Island")
CopTimeseries(CI_final, "Clarence Island")

# -------------------------Step 5: Make Sea Ice Data Timeseries------------------------
# Function to create timeseries by site
IceTimeseries <- function(data, site) {
  data$sice_conc <- data$sice_conc_mean * 100
  data$date <- as.Date(data$date)
  data <- filter(data, depth == 0.5)
  # sea ice thickness
  thickness <- ggplot(data = data, mapping = aes(x = date, y = sice_thick_mean)) + geom_line(color = "navy", linewidth = 1) + 
    labs(x = "Sea Ice Thickness", y = "meters") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # sea ice concentration
  concentration <- ggplot(data = data, mapping = aes(x = date, y = sice_conc)) + 
    geom_line(color = "blue", linewidth = 1) + 
    labs(x = "Sea Ice Concentration", y = "%") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # northward velocity
  nvelocity <- ggplot(data = data, mapping = aes(x = date, y = sice_n_veloc_mean)) + geom_line(color = "mediumslateblue", linewidth = 1) + 
    labs(x = "Sea Ice North Velocity", y = "m/s") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # eastward velocity
  evelocity <- ggplot(data = data, mapping = aes(x = date, y = sice_e_veloc_mean)) + 
    geom_line(color = "dodgerblue", linewidth = 1) + 
    labs(x = "Sea Ice East Velocity", y = "m/s") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  # Arrange plots to one figure
  windows()
  final_plot <- grid.arrange(thickness, concentration, nvelocity, evelocity, nrow = 4, top = site)
  
}
# Construct timeseries for all sites
IceTimeseries(EI_final, "Elephant Island")
IceTimeseries(KGI_final, "King George Island")
IceTimeseries(CI_final, "Clarence Island")

# --------------------- Step 6: Save Data ---------------
EI_final$Site <- 'EI'
KGI_final$Site <- 'KGI'
CI_final$Site <- 'CI'

final <- rbind(EI_final,KGI_final,CI_final)
write.csv(final, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/copernicus.csv") # write destination