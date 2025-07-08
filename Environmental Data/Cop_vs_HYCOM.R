library(R.matlab)
library(lubridate)
library(tidyverse)
library(gridExtra)
library(CopernicusMarine)
library(stars)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(gridExtra) # for grid.arrange

# --------------------------Step 1: Read HYCOM ----------------------------
readHYCOM <- function(site, start, end) {
  start <- as.Date(start)
  end <- as.Date(end)
  dates <- seq(start, end, by = "day")
  
  df <- data.frame()
  
  for (current in dates) {
    # Format date into yyyymmdd
    date_str <- format(as.Date(current), "%Y%m%d")
    
    # Construct filename (you can make site-specific filename logic here if needed)
    if(site == 'EI'){
      file_path <- paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/HYCOM/Data/W303E304S61N59_",
                          date_str, "T0000Z.mat")
    } else if(site == 'KGI') {
      file_path <- paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/HYCOM/Data/W303E304S62N61_",
                          date_str, "T0000Z.mat")
    } else {
      file_path <- paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/HYCOM/Data/W306E307S62N61_",
                          date_str, "T0000Z.mat")
    }
    
    if (!file.exists(file_path)) {
      warning(paste("File does not exist:", file_path))
      next
    }
    
    # Read the .mat file
    mat <- readMat(file_path)
    
    # Extract relevant variables (adapt these to match the actual structure of your .mat files)
    lon <- as.vector(mat$lon)
    lat <- as.vector(mat$lat)
    depth <- as.vector(mat$depth)
    temp <- mat$temp
    salt <- mat$salt
    
    # Assume temp and salt are 3D arrays: lon x lat x depth
    # Flatten arrays to data frame
    nlon <- length(lon)
    nlat <- length(lat)
    ndepth <- length(depth)
    
    grid <- expand.grid(lon = lon, lat = lat, depth = depth)
    grid$temp <- as.vector(temp)
    grid$salt <- as.vector(salt)
    grid$date <- as.Date(current)
    
    # Average over space for each depth
    daily_avg <- grid %>%
      group_by(depth, date) %>%
      summarise(temp = mean(temp, na.rm = TRUE),
                salt = mean(salt, na.rm = TRUE),
                .groups = "drop")

    df <- bind_rows(df, daily_avg)
  }
  df$date <- as.POSIXct(df$date, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  return(df)
  
}
EI_HYCOM <- readHYCOM('EI','2014-03-05','2014-07-17')
KGI_HYCOM <- readHYCOM('KGI','2015-02-10','2016-01-29') 
CI_HYCOM <- readHYCOM('CI','2016-02-04','2016-12-02')


# -----------------Step 2: Load Copernicus Data-------------------
# Load Copernicus data
EI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/EI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990714294.nc")
KGI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/KGI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990804735.nc")
CI_cop <- nc_open("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/CI_ALLsurf_cmems_mod_glo_phy_my_0.083deg_P1D-m_1751990301834.nc")

dfFromNC <- function(site, data) {
  # Extracting relevant data
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  temp <- ncvar_get(data, 'thetao')
  salinity <- ncvar_get(data, 'so')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  
  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   sst = as.vector(temp), salinity = as.vector(salinity))
  
  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(sst = mean(sst, na.rm=TRUE), salinity = mean(salinity, na.rm=TRUE))
  
  return(avg_df)
}
# Constructing site dataframes
EI_copDF <- dfFromNC("EI",EI_cop)
KGI_copDF <- dfFromNC("KGI",KGI_cop)
CI_copDF <- dfFromNC("CI",CI_cop)


# -------------------- Step 3: Salinity/Temp Timeseries--------------
cop_HYCOM_ts <- function(cop, hycom, site) {
  # Extracting surface HYCOM values and seeing NAs
  hycom_surf <- filter(hycom, depth == 0)
  print(paste(sum(is.na(cop$sst)), "NA temp values and ", sum(is.na(cop$salinity)), 
              "NA salinity values for Copernicus at site ",site))
  print(paste(sum(is.na(hycom_surf$temp)), "NA temp values and ", sum(is.na(hycom_surf$salt)), 
              "NA salinity values for Copernicus at site ",site))
  
  # Combining two sources into one dataframe
  combined <- data.frame(cbind(date = cop$date, cop_temp = cop$sst, cop_sal = cop$salinity,
                    hycom_temp = hycom_surf$temp, hycom_sal = hycom_surf$salt))
  combined$date <- as.POSIXct(combined$date, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  combined$date <- as.Date(combined$date, "%Y-%m-%d", tx = "UTC")
  
  # Making combined timeseries plot
  temp <- ggplot(data = combined, aes(date)) + geom_line(aes(y=cop_temp,color = 'Copernicus')) +
    geom_line(aes(y=hycom_temp, color='HYCOM')) + labs(y = "Temperature (C)", color = "Source", x = NULL) +
    scale_x_date(date_labels = "%b %Y") + scale_color_manual(values = c("Copernicus"='red', "HYCOM"='blue'))
  salinity <- ggplot(data = combined, aes(date)) + geom_line(aes(y=cop_sal,color = 'Copernicus')) +
    geom_line(aes(y=hycom_sal, color='HYCOM')) + labs(y = "Salinity (psu)", color = "Source", x = NULL) +
    scale_x_date(date_labels = "%b %Y") + scale_color_manual(values = c("Copernicus"='red', "HYCOM"='blue'))
  final <- grid.arrange(temp, salinity, nrow = 2, top = site)
}
# Making timeseries for each site
cop_HYCOM_ts(EI_copDF, EI_HYCOM, 'Elephant Island')
cop_HYCOM_ts(KGI_copDF, KGI_HYCOM, 'King George Island')
cop_HYCOM_ts(CI_copDF, CI_HYCOM, 'Clarence Island')