library(CopernicusMarine)
library(stars)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(tidyverse)
library(gridExtra) # for grid.arrange

# -----------------Step 2: Extract Oceanographic Data-------------------
# Load Copernicus data
# all physical vars 
# downloaded from: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/download
APEI <- nc_open("L:/Shared drives/Antarctic Marine Mammals/Krill Data/cmems_mod_glo_phy_my_0.083deg_P1D-m_1763418042285.nc")
APDE <- nc_open("L:/Shared drives/Antarctic Marine Mammals/Krill Data/cmems_mod_glo_phy_my_0.083deg_P1D-m_1763440909642.nc")

# Ensuring that variable dimensions index correctly
fix_grid <- function(var) {
  d <- length(dim(var))
  
  if (is.null(d)) {
    # 1D vector
    return(as.vector(var))
  }
  
  if (d == 3) {
    # [lon, lat, time] -> keep as is
    var <- aperm(var, c(1, 2, 3))
  } else if (d == 2) {
    # [lon, lat] -> keep as is
    var <- aperm(var, c(1, 2))
  }
  
  return(as.vector(var))
}

physFromNC <- function(data) {
  # Dimensions
  time <- ncvar_get(data, "time")
  lat  <- ncvar_get(data, "latitude")
  lon  <- ncvar_get(data, "longitude")
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")
  
  # Oceanographic variables (all 2D or 3D: lon x lat x time)
  temp        <- fix_grid(ncvar_get(data, "thetao"))
  ssh         <- fix_grid(ncvar_get(data, "zos"))
  n_velocity  <- fix_grid(ncvar_get(data, "vo"))
  e_velocity  <- fix_grid(ncvar_get(data, "uo"))
  salinity    <- fix_grid(ncvar_get(data, "so"))
  mixed_layer <- fix_grid(ncvar_get(data, "mlotst"))
  
  # Sea ice variables
  sice_conc    <- fix_grid(ncvar_get(data, "siconc"))
  sice_thick   <- fix_grid(ncvar_get(data, "sithick"))
  sice_e_veloc <- fix_grid(ncvar_get(data, "usi"))
  sice_n_veloc <- fix_grid(ncvar_get(data, "vsi"))
  
  # Build lon-lat-time grid (no depth)
  lonlattime <- expand.grid(lon = lon, lat = lat, date = time_obs)
  
  df <- data.frame(
    lon   = lonlattime$lon,
    lat   = lonlattime$lat,
    date  = lonlattime$date,
    ssh   = ssh,
    n_velocity  = n_velocity,
    e_velocity  = e_velocity,
    salinity    = salinity,
    mixed_layer = mixed_layer,
    sice_conc    = sice_conc,
    sice_thick   = sice_thick,
    sice_e_veloc = sice_e_veloc,
    sice_n_veloc = sice_n_veloc,
    temp        = temp
  )
  
  # Daily *spatial* averages (no depth term)
  avg_df <- df %>%
    group_by(date) %>%
    summarize(
      ssh_mean         = mean(ssh, na.rm = TRUE),
      ssh_sd           = sd(ssh, na.rm = TRUE),
      n_velocity_mean  = mean(n_velocity, na.rm = TRUE),
      n_velocity_mad   = mad(n_velocity, na.rm = TRUE),
      e_velocity_mean  = mean(e_velocity, na.rm = TRUE),
      e_velocity_mad   = mad(e_velocity, na.rm = TRUE),
      salinity_mean    = mean(salinity, na.rm = TRUE),
      salinity_sd      = sd(salinity, na.rm = TRUE),
      mixed_layer_mean = mean(mixed_layer, na.rm = TRUE),
      mixed_layer_sd   = sd(mixed_layer, na.rm = TRUE),
      sice_conc_mean   = mean(sice_conc, na.rm = TRUE),
      ice_conc_sd      = sd(sice_conc, na.rm = TRUE),
      sice_thick_mean  = mean(sice_thick, na.rm = TRUE),
      sice_e_veloc_mean = mean(sice_e_veloc, na.rm = TRUE),
      sice_n_veloc_mean = mean(sice_n_veloc, na.rm = TRUE),
      temp_mean        = mean(temp, na.rm = TRUE),
      temp_sd          = sd(temp, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(avg_df)
}

# Constructing site dataframes
EI_phys <- physFromNC(APEI)
KGI_phys <- physFromNC(APDE)

# -------------------- Step 3: Make all depth dataframe --------------
depthMatchALL <- function(data, type) {
  # matching depths from biogeochem and ocean physics data
  data$depth[data$depth %in% c(0.506, 0.494)] <- 0.5
  data$depth[data$depth %in% c(1.556, 1.541)] <- 1.5
  data$depth[data$depth %in% c(2.668, 2.646)] <- 2.6
  data$depth[data$depth %in% c(3.856, 3.819)] <- 3.8
  data$depth[data$depth %in% c(5.140, 5.078)] <- 5
  data$depth[data$depth %in% c(6.543, 6.441)] <- 6.5
  data$depth[data$depth %in% c(8.093, 7.930)] <- 8
  data$depth[data$depth %in% c(9.823, 9.573)] <- 10
  data$depth[data$depth %in% c(11.774, 11.405)] <- 11
  data$depth[data$depth %in% c(13.991, 13.467)] <- 14
  data$depth[data$depth %in% c(16.525, 15.810)] <- 16
  data$depth[data$depth %in% c(19.430, 18.496)] <- 19
  data$depth[data$depth %in% c(22.758, 21.599)] <- 22
  data$depth[data$depth %in% c(26.558, 25.211)] <- 26
  data$depth[data$depth %in% c(30.875, 29.445)] <- 30
  data$depth[data$depth %in% c(35.740, 34.434)] <- 35
  data$depth[data$depth %in% c(41.180, 40.344)] <- 41
  data$depth[data$depth %in% c(47.212, 47.374)] <- 47
  data$depth[data$depth %in% c(53.851, 55.764)] <- 54
  data$depth[data$depth %in% c(69.022, 65.807)] <- 67 
  data$depth[data$depth %in% c(77.611, 77.854)] <- 77
  data$depth[data$depth %in% c(97.041, 92.326)] <- 95 
  data$depth[data$depth %in% c(108.030, 109.729)] <- 109
  data$depth[data$depth %in% c(133.076, 130.666)] <- 130
  data$depth[data$depth %in% c(180.550, 186.126)] <- 183
  data$depth[data$depth %in% c(221.141, 222.475)] <- 222
  data$depth[data$depth %in% c(271.356, 266.040)] <- 270
  data$depth[data$depth %in% c(333.863, 318.127)] <- 325
  data$depth[data$depth %in% c(370.689, 380.213)] <- 375
  data$depth[data$depth %in% c(457.626, 453.938)] <- 455
  data$depth[data$depth %in% c(565.292, 541.089)] <- 555
  data$depth[data$depth %in% c(628.026, 643.567)] <- 635
  data$depth[data$depth %in% c(773.368, 763.333)] <- 768
  data$depth[data$depth %in% c(947.448, 902.339)] <- 920
  data$depth[data$depth %in% c(1045.854, 1062.440)] <- 1055
  data$depth[data$depth %in% c(1265.861, 1245.291)] <- 1255
  data$depth[data$depth %in% c(1516.364, 1452.251)] <- 1585
  data$depth[data$depth %in% c(1652.568, 1684.284)] <- 1665
  # Depths in biogeochem that are not in physics:
  #    61.113, 86.929 120.000 147.406 163.165 199.790 244.891 300.888 
  #    411.794 508.640 697.259 856.679 1151.991 1387.377
  return(data)
}
EI_depths <- left_join(depthMatchALL(EI_biodf),depthMatchALL(EI_phys), by = c('date', 'depth'))
KGI_depths <- left_join(depthMatchALL(KGI_biodf),depthMatchALL(KGI_phys), by = c('date', 'depth'))
CI_depths <- left_join(depthMatchALL(CI_biodf),depthMatchALL(CI_phys), by = c('date', 'depth'))

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
EI_depths <- getEKE(EI_depths)
KGI_depths <- getEKE(KGI_depths)
CI_depths <- getEKE(CI_depths)

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
EI_depths$Site <- 'EI'
KGI_depths$Site <- 'KGI'
CI_depths$Site <- 'CI'

final <- rbind(EI_final,KGI_final,CI_final)
final_depths <- rbind(EI_depths,KGI_depths,CI_depths)
write.csv(final, "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/copernicus.csv") # write destination
write.csv(final_depths, "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/copernicus_depths.csv") # write destination