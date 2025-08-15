library(CopernicusMarine)
library(stars)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(tidyverse)
library(gridExtra) # for grid.arrange

# ---------------- Step 0: Subsetting and Loading Data -------------------
# Accessing Copernicus surface primary proudction, chlorophyll, and EKE at time lags of 1-6 months

# Date and grid range for EI site:
# site date, coordinate: 3/5/14 to 7/17/14, -60.8869 lat 55.95396667 lon
# date bounds (6 months before start to 1 month before end), month length = 30 days
#   9/6/13 to 6/17/14
# 40 km bounding box:
#    latitude: -60.52602,-61.24778
#    longitude: -56.69255, -55.21545

# Date and grid range for KGI site:
# site date, coordinate: 2/10/15 to 1/29/16, -61.457817 lat 57.94191667 lon
# date bounds (6 months before start to 1 month before end), month length = 30 days
#   8/14/14 to 12/30/15
# 40 km bounding box:
#    latitude: -61.09694, -61.8187
#    longitude: -58.69396, -57.18988

# Date and grid range for CI site:
# site date, coordinate: 2/4/16 to 12/2/16, -61.25186667 lat 53.48343333 lon
# date bounds (6 months before start to 1 month before end), month length = 30 days
#   8/8/15 to 11/2/16
# 40 km bounding box:
#    latitude: -60.89099, -61.61275
#    longitude: -54.23054, -52.73632

# temperature, salinity, velocities (to derive EKE)
# downloaded from: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/download
EI_phys <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/lagged/EI-lag_cmems_mod_glo_phy_my_0.083deg_P1D-m_1755194554475.nc')
KGI_phys <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/lagged/KGI-lag_cmems_mod_glo_phy_my_0.083deg_P1D-m_1755130224784.nc')
CI_phys <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/lagged/CI-lag_cmems_mod_glo_phy_my_0.083deg_P1D-m_1755130273890.nc')

# primary production, chlorophyll
# downloaded from: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/download 
EI_bio <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/lagged/EI-lag_cmems_mod_glo_bgc_my_0.25deg_P1D-m_1755197529356.nc')
KGI_bio <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/lagged/KGI-lag_cmems_mod_glo_bgc_my_0.25deg_P1D-m_1755129104796.nc')
CI_bio <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/lagged/CI-lag_cmems_mod_glo_bgc_my_0.25deg_P1D-m_1755129164010.nc')

# ------------------- Step 1: Prepping Data ---------------
physFromNC <- function(data) {
  # Extracting relevant data
  # Dimensions
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time

  # Oceanographic variables
  temp <- ncvar_get(data,'thetao')
  n_velocity <- ncvar_get(data, 'vo')
  e_velocity <- ncvar_get(data,'uo')
  salinity <- ncvar_get(data,'so')

  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   n_velocity = as.vector(n_velocity),
                   e_velocity = as.vector(e_velocity), salinity = as.vector(salinity),
                   temp=as.vector(temp))
  
  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(n_velocity = mean(n_velocity, na.rm=TRUE),
                                                       e_velocity = mean(e_velocity, na.rm=TRUE), 
                                                       salinity = mean(salinity, na.rm=TRUE), 
                                                       temp =mean(temp, na.rm=TRUE))
  avg_df <- ungroup(avg_df)
  return(avg_df)
}
# Constructing site dataframes
EI_physdf <- physFromNC(EI_phys)
KGI_physdf <- physFromNC(KGI_phys)
CI_physdf <- physFromNC(CI_phys)

# Similar function to above, for biogeochem vars
bioFromNC <- function(data) {
  # Extracting relevant data
  # Dimensions
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  
  # biogeochemical variables
  chla <- ncvar_get(data, 'chl')
  productivity <- ncvar_get(data, 'nppv')
  
  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   chla = as.vector(chla), productivity = as.vector(productivity))
  
  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(chla = mean(chla, na.rm=TRUE), 
                                                       productivity = mean(productivity, na.rm=TRUE))
  avg_df <- ungroup(avg_df)
  return(avg_df)
}
# Constructing site dataframes
EI_biodf <- bioFromNC(EI_bio)
KGI_biodf <- bioFromNC(KGI_bio)
CI_biodf <- bioFromNC(CI_bio)

# joining dataframes
EI_final <- left_join(EI_biodf,EI_physdf, by = c('date'))
EI_final$site <- 'EI'
KGI_final <- left_join(KGI_biodf,KGI_physdf, by = c('date'))
KGI_final$site <- 'KGI'
CI_final <- left_join(CI_biodf,CI_physdf, by = c('date'))
CI_final$site <- 'CI'

# joining all sites to one dataframe
unlagged <- rbind(EI_final, KGI_final, CI_final)
unlagged$date <- as.Date(unlagged$date)


# --------------- Step 2: Match with site dates -------------
# setting dates and site ID for each site
date <- seq(from = as.Date('03-05-2014', format = '%m-%d-%Y'),
             to = as.Date('07-17-2014', format = '%m-%d-%Y'), by = 1)
date <- append(date, seq(from = as.Date('02-10-2015', format = '%m-%d-%Y'),
                           to = as.Date('01-29-2016', format = '%m-%d-%Y'), by = 1))
date <- append(date, seq(from = as.Date('02-04-2016', format = '%m-%d-%Y'),
                           to = as.Date('12-02-2016', format = '%m-%d-%Y'), by = 1))
site <- c(rep("EI", 135), rep("KGI", 354), rep("CI", 303))

# making basis data frame for each site
lagged <- data.frame(date, site)

# function for adding lagged variables
lagVar <- function(var) {
  # running through for 1 to 6 month lags
  for(lag in 1:6) {
    if(lag == 1) {
      # 30 day date lags
      lag_dates_EI <- seq(from = as.Date('02-03-2014', format = '%m-%d-%Y'),
                          to = as.Date('06-17-2014', format = '%m-%d-%Y'), by = 1)
      lag_dates_KGI <-seq(from = as.Date('01-11-2015', format = '%m-%d-%Y'),
                          to = as.Date('12-30-2015', format = '%m-%d-%Y'), by = 1)
      lag_dates_CI <- seq(from = as.Date('01-05-2016', format = '%m-%d-%Y'),
                          to = as.Date('11-02-2016', format = '%m-%d-%Y'), by = 1)
      
    } else if (lag == 2) {
      # 60 day date lags
      lag_dates_EI <- seq(from = as.Date('01-04-2014', format = '%m-%d-%Y'),
                          to = as.Date('05-18-2014', format = '%m-%d-%Y'), by = 1)
      lag_dates_KGI <-seq(from = as.Date('12-12-2014', format = '%m-%d-%Y'),
                          to = as.Date('11-30-2015', format = '%m-%d-%Y'), by = 1)
      lag_dates_CI <- seq(from = as.Date('12-06-2015', format = '%m-%d-%Y'),
                          to = as.Date('10-03-2016', format = '%m-%d-%Y'), by = 1)
      
    } else if (lag == 3) {
      # 90 day date lags
      lag_dates_EI <- seq(from = as.Date('12-05-2013', format = '%m-%d-%Y'),
                          to = as.Date('04-18-2014', format = '%m-%d-%Y'), by = 1)
      lag_dates_KGI <-seq(from = as.Date('11-12-2014', format = '%m-%d-%Y'),
                          to = as.Date('10-31-2015', format = '%m-%d-%Y'), by = 1)
      lag_dates_CI <- seq(from = as.Date('11-06-2015', format = '%m-%d-%Y'),
                          to = as.Date('09-03-2016', format = '%m-%d-%Y'), by = 1)
      
    } else if (lag == 4) {
      # 120 day date lags
      lag_dates_EI <- seq(from = as.Date('11-05-2013', format = '%m-%d-%Y'),
                          to = as.Date('03-19-2014', format = '%m-%d-%Y'), by = 1)
      lag_dates_KGI <-seq(from = as.Date('10-13-2014', format = '%m-%d-%Y'),
                          to = as.Date('10-01-2015', format = '%m-%d-%Y'), by = 1)
      lag_dates_CI <- seq(from = as.Date('10-07-2015', format = '%m-%d-%Y'),
                          to = as.Date('08-04-2016', format = '%m-%d-%Y'), by = 1)
      
    } else if (lag == 5) {
      # 150 day date lags
      lag_dates_EI <- seq(from = as.Date('10-06-2013', format = '%m-%d-%Y'),
                          to = as.Date('02-17-2014', format = '%m-%d-%Y'), by = 1)
      lag_dates_KGI <-seq(from = as.Date('09-13-2014', format = '%m-%d-%Y'),
                          to = as.Date('09-01-2015', format = '%m-%d-%Y'), by = 1)
      lag_dates_CI <- seq(from = as.Date('09-07-2015', format = '%m-%d-%Y'),
                          to = as.Date('07-05-2016', format = '%m-%d-%Y'), by = 1)
      
    } else if (lag == 6) {
      # 180 day date lags
      lag_dates_EI <- seq(from = as.Date('09-06-2013', format = '%m-%d-%Y'),
                          to = as.Date('01-18-2014', format = '%m-%d-%Y'), by = 1)
      lag_dates_KGI <-seq(from = as.Date('08-14-2014', format = '%m-%d-%Y'),
                          to = as.Date('08-02-2015', format = '%m-%d-%Y'), by = 1)
      lag_dates_CI <- seq(from = as.Date('08-08-2015', format = '%m-%d-%Y'),
                          to = as.Date('06-05-2016', format = '%m-%d-%Y'), by = 1)
      
    }
    
    # get lagged data
    lagged_data <- doLags(lag_dates_EI, lag_dates_KGI, lag_dates_CI, var, lag)
    colName <- paste0(var,'_',lag,'mon')    # name for lagged column
    
    # adding specific month lag to existing dataframe
    lagged <- cbind(lagged, lagged_data)
    
    # renaming column to reflect variable name and lag
    lagged <- lagged %>% mutate(!!colName := lagged_data)
    lagged <- lagged %>% subset(select = -lagged_data)
    
    
  }
  
  return(lagged)
}

# function to create column for specified variable and lag
doLags <- function(EI_date, KGI_date, CI_date, var, lag) {
  lag_dates <- c(EI_date, KGI_date, CI_date)
  
  # filter out repeated dates across sites and add dataframes together
  EI_filtered <- unlagged %>% filter(date %in% EI_date, site == 'EI')
  KGI_filtered <- unlagged %>% filter(date %in% KGI_date, site == 'KGI')
  CI_filtered <- unlagged %>% filter(date %in% CI_date, site == 'CI')
  filtered <- rbind(EI_filtered, KGI_filtered, CI_filtered)
  
  # creating list with values of the lagged variables
  varVal <- c(NULL)
  for(d in lag_dates) {
    as.Date(d)
    rowIdx <- which(filtered$date == d)
    val <- pull(filtered[rowIdx, var])
    varVal <- append(varVal,val)
  }
  
  return(varVal)
}


lagged <- lagVar('chla')
lagged <- lagVar('productivity')
lagged <- lagVar('n_velocity')
lagged <- lagVar('e_velocity')
lagged <- lagVar('salinity')
lagged <- lagVar('temp')

getEKE <- function(north, east) {
  u_mean <- mean(east,na.rm=TRUE)
  v_mean <- mean(north, na.rm=TRUE)
  
  u_prime <- east - u_mean
  v_prime <- north - v_mean
  
  # EKE = 1/2((u')^2 + (v')^2)
  # units: cm^2/s^2
  EKE <- (0.5 * ((u_prime^2) + (v_prime^2))) * 10000
  return(EKE)
}

for(lag in 1:6) {
  north <- paste0('n_velocity_',lag,'mon')
  east <- paste0('e_velocity_',lag,'mon')
  EKE <- paste0('EKE_',lag,'mon')
  lagged[,EKE] <- getEKE(lagged[,north], lagged[,east])
  lagged <- lagged %>% select(-all_of(c(north, east)))
}

write.csv(lagged, '/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/lagged/copernicusLagged.csv')

# --------------- Step 3: Make timeseries ------------------
# To generate lag timeseries, call makePlots(var,site_name)
#   var = variable name as string, as it appears in lagged dataframe
#   site_name = site name as string
#   example call: makePlots('productivity','EI')

makePlots <- function(var,site_name) {
  if(var == 'chla') {
    name <- 'Chlorophyll'
    unit <- 'mg/m3'
  } else if(var == 'productivity') {
    name <- 'Net Primary Production'
    unit <- 'mg/m3 carbon'
  } else if(var == 'temp') {
    name <- 'Surface Temperature'
    unit <- '°C'
  } else if(var == 'salinity') {
    name <- 'Surface Salinity'
    unit <- 'psu'
  } else if(var == 'EKE') {
    name <- 'Surface EKE'
    unit <- 'cm2/s2'
  }
  
  site_df <- lagged %>% filter(site == site_name)
  
  p1 <- ggplot(data = site_df, mapping = aes(x = date, y = .data[[paste0(var,'_1mon')]])) + 
    geom_line(, linewidth = 1) + 
    labs(x = paste0(name," 1 Month Lag"), y = unit) + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  p2 <- ggplot(data = site_df, mapping = aes(x = date, y = .data[[paste0(var,'_2mon')]])) + 
    geom_line(, linewidth = 1) + 
    labs(x = paste0(name," 2 Month Lag"), y = unit) + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  p3 <- ggplot(data = site_df, mapping = aes(x = date, y = .data[[paste0(var,'_3mon')]])) + 
    geom_line(, linewidth = 1) + 
    labs(x = paste0(name," 3 Month Lag"), y = unit) + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  p4 <- ggplot(data = site_df, mapping = aes(x = date, y = .data[[paste0(var,'_4mon')]])) + 
    geom_line(, linewidth = 1) + 
    labs(x = paste0(name," 4 Month Lag"), y = unit) + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  p5 <- ggplot(data = site_df, mapping = aes(x = date, y = .data[[paste0(var,'_5mon')]])) + 
    geom_line(, linewidth = 1) + 
    labs(x = paste0(name," 5 Month Lag"), y = unit) + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  p6 <- ggplot(data = site_df, mapping = aes(x = date, y = .data[[paste0(var,'_6mon')]])) + 
    geom_line(, linewidth = 1) + 
    labs(x = paste0(name," 6 Month Lag"), y = unit) + scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(1, 0.5, 1, 0.5), units = "line"))
  
  final_plot <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 6, top = site_name)
  return(final_plot)
}