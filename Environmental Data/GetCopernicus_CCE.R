library(CopernicusMarine)
library(stars)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(tidyverse)
library(gridExtra) # for grid.arrange

# -----------------Step 2: Extract Data-------------------
# Load Copernicus data
# surface data, all vars
# EI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/EI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990714294.nc")
# KGI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/KGI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990804735.nc")
# CI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/CI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990301834.nc")

# all vars 
EI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/EI_cmems_mod_glo_phy_my_0.083deg_P1D-m_1752518931575.nc")
KGI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/KGI_cmems_mod_glo_phy_my_0.083deg_P1D-m_1752529248292.nc")
CI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/CI_cmems_mod_glo_phy_my_0.083deg_P1D-m_1752529427864.nc")

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

dfFromNC <- function(data) {
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

  # Creating a dataframe with daily spatial averages
  avg_df <- filtered %>% group_by(date, depth) %>% summarize(ssh = mean(ssh, na.rm=TRUE), n_velocity = mean(n_velocity, na.rm=TRUE),
                                                e_velocity = mean(e_velocity, na.rm=TRUE), salinity = mean(salinity, na.rm=TRUE),
                                                mixed_layer = mean(mixed_layer, na.rm=TRUE), sice_conc = mean(sice_conc, na.rm=TRUE),
                                                sice_thick = mean(sice_thick, na.rm=TRUE), sice_e_veloc = mean(sice_e_veloc, na.rm=TRUE),
                                                sice_n_veloc = mean(sice_n_veloc, na.rm=TRUE), temp=mean(temp, na.rm=TRUE))
  avg_df <- ungroup(avg_df)
  return(avg_df)
}
# Constructing site dataframes
EI_df <- dfFromNC(EI_cop)
KGI_df <- dfFromNC(KGI_cop)
CI_df <- dfFromNC(CI_cop)

# --------------------------Step 3: Compute EKE--------------------------
getEKE <- function(data) {
  # get u' and v' for EKE equation
  data <- data %>% group_by(depth) %>% mutate(u_mean = mean(e_velocity,na.rm=TRUE),
                                                   v_mean = mean(n_velocity, na.rm=TRUE))
  data <- data %>% group_by(depth) %>% mutate(u_prime = e_velocity-u_mean, 
                                                 v_prime = n_velocity-v_mean)
  
  # EKE = 1/2((u')^2 + (v')^2)
  # units: cm^2/s^2
  data$EKE <- (0.5 * ((data$u_prime^2) + (data$v_prime^2))) * 10000
  data <- subset(data, select = -c(u_mean, v_mean, u_prime,v_prime))
  data <- data %>% ungroup()
  return(data)
}
EI_df <- getEKE(EI_df)
KGI_df <- getEKE(KGI_df)
CI_df <- getEKE(CI_df)

# -------------------------Step 4: Make Oceanographic Data Timeseries------------------------
# Function to create timeseries by site (for surface only)
CopTimeseries <- function(data, site) {
  data$date <- as.Date(data$date)
  data <- filter(data, depth == 0.494)
  # sea surface height
  ssh <- ggplot(data = data, mapping = aes(x = date, y = ssh)) + geom_line(color = "tomato", size = 1) + 
    labs(x = "Sea Surface Height", y = "meters") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # eddy kinetic energy
    EKE <- ggplot(data = data, mapping = aes(x = date, y = EKE)) + 
    geom_line(color = "mediumpurple", size = 1) + 
    labs(x = "Eddy Kinetic Energy", y = "cm^2/s^2") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # salinity
  salinity <- ggplot(data = data, mapping = aes(x = date, y = salinity)) + geom_line(color = "gold", size = 1) + 
    labs(x = "Salinity", y = "psu") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # mixed layer depth
  mixed <- ggplot(data = data, mapping = aes(x = date, y = mixed_layer)) + 
    geom_line(color = "mediumseagreen", size = 1) + 
    labs(x = "Mixed Layer Thickness", y = "meters") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  # temperature
  temp <- ggplot(data = data, mapping = aes(x = date, y = temp)) + 
    geom_line(color = "darkred", size = 1) + 
    labs(x = "Temperature", y = "C") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(.5, 0.5, .5, 0.5), units = "line"))
  
  # Arrange plots to one figure
  windows()
  final_plot <- grid.arrange(ssh, EKE, salinity, mixed, temp, nrow = 5, top = site)
  
}
# Construct timeseries for all sites
CopTimeseries(EI_df, "Elephant Island")
CopTimeseries(KGI_df, "King George Island")
CopTimeseries(CI_df, "Clarence Island")

# -------------------------Step 5: Make Sea Ice Data Timeseries------------------------
# Function to create timeseries by site
IceTimeseries <- function(data, site) {
  data$sice_conc <- data$sice_conc * 100
  data$date <- as.Date(data$date)
  data <- filter(data, depth == 0.494)
  # sea ice thickness
  thickness <- ggplot(data = data, mapping = aes(x = date, y = sice_thick)) + geom_line(color = "navy", size = 1) + 
    labs(x = "Sea Ice Thickness", y = "meters") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # sea ice concentration
  concentration <- ggplot(data = data, mapping = aes(x = date, y = sice_conc)) + 
    geom_line(color = "blue", size = 1) + 
    labs(x = "Sea Ice Concentration", y = "%") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # northward velocity
  nvelocity <- ggplot(data = data, mapping = aes(x = date, y = sice_n_veloc)) + geom_line(color = "mediumslateblue", size = 1) + 
    labs(x = "Sea Ice North Velocity", y = "m/s") + scale_x_date(date_labels = "%b %Y") + 
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  # eastward velocity
  evelocity <- ggplot(data = data, mapping = aes(x = date, y = sice_e_veloc)) + 
    geom_line(color = "dodgerblue", size = 1) + 
    labs(x = "Sea Ice East Velocity", y = "m/s") + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  # Arrange plots to one figure
  windows()
  final_plot <- grid.arrange(thickness, concentration, nvelocity, evelocity, nrow = 4, top = site)
  
}
# Construct timeseries for all sites
IceTimeseries(EI_df, "Elephant Island")
IceTimeseries(KGI_df, "King George Island")
IceTimeseries(CI_df, "Clarence Island")

# --------------------- Step 6: Save Data ---------------
EI_df$Site <- 'EI'
KGI_df$Site <- 'KGI'
CI_df$Site <- 'CI'

final <- rbind(EI_df,KGI_df,CI_df)
write.csv(final, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/copernicus.csv") # write destination