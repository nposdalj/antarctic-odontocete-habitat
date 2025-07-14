library(R.matlab)
library(lubridate)
library(tidyverse)
library(gridExtra)
library(CopernicusMarine)
library(stars)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(gridExtra) # for grid.arrange
library(multDM) # For Diebold-Mariono test
library(ggpubr)

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
# Loading remotely sensed data
# commenting out salinity because values seemed unreasonable
# EI_sal <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Salinity_aqua/ EI_salinity.csv")
# KGI_sal <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Salinity_aqua/ KGI_salinity.csv")
# CI_sal <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Salinity_aqua/ CI_salinity.csv")
EI_aqua <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/SST_ERDDAP/EI_SST.csv")
KGI_aqua <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/SST_ERDDAP/KGI_SST.csv")
CI_aqua <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/SST_ERDDAP/CI_SST.csv")
# Merging into one
# EI_aqua <- merge(EI_temp, EI_sal, by=intersect(names(EI_temp), names(EI_sal)))
# KGI_aqua <- merge(KGI_temp, KGI_sal, by=intersect(names(KGI_temp), names(KGI_sal)))
# CI_aqua <- merge(CI_temp, CI_sal, by=intersect(names(CI_temp), names(CI_sal)))

cop_HYCOM_ts <- function(cop, hycom, aqua, site) {
  # Extracting surface HYCOM values and seeing NAs
  hycom_surf <- filter(hycom, depth == 0)
  print(paste(sum(is.na(cop$sst)), "NA temp values and ", sum(is.na(cop$salinity)), 
              "NA salinity values for Copernicus at site ",site))
  print(paste(sum(is.na(hycom_surf$temp)), "NA temp values and ", sum(is.na(hycom_surf$salt)), 
              "NA salinity values for HYCOM at site ",site))
  print(paste(sum(is.na(aqua$sst)), "NA sst values for Aqua MODIS at site ",site))
  
  # Combining two sources into one dataframe
  combined <- data.frame(cbind(date = cop$date, cop_temp = cop$sst, cop_sal = cop$salinity,
                    hycom_temp = hycom_surf$temp, hycom_sal = hycom_surf$salt, aqua_temp = aqua$sst))
  combined$date <- as.POSIXct(combined$date, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  combined$date <- as.Date(combined$date, "%Y-%m-%d", tx = "UTC")
  
  # Making combined timeseries plot
  temp <- ggplot(data = combined, aes(date)) + geom_line(aes(y=cop_temp,color = 'Copernicus')) +
    geom_line(aes(y=hycom_temp, color='HYCOM')) + 
    geom_point(aes(y=aqua_temp, color='Aqua MODIS'),size=1.5) +
    labs(y = "Temperature (C)", color = "Source", x = NULL) +
    scale_x_date(date_labels = "%b %Y") + ylim(-2,4) +
    scale_color_manual(values = c("Copernicus"='red', "HYCOM"='blue', 'Aqua MODIS' = 'darkgreen'))
  
  salinity <- ggplot(data = combined, aes(date)) + geom_line(aes(y=cop_sal,color = 'Copernicus')) +
    geom_line(aes(y=hycom_sal, color='HYCOM')) +  ylim(33,35) +
    labs(y = "Salinity (psu)", color = "Source", x = NULL) +
    scale_x_date(date_labels = "%b %Y") + scale_color_manual(values = c("Copernicus"='red', "HYCOM"='blue', 'aqua'='darkgreen'))
  final <- grid.arrange(temp, salinity, nrow = 2, top = site)
}
# geom_line(aes(y=aqua_sal, color = 'aqua')) +
# Making timeseries for each site
cop_HYCOM_ts(EI_copDF, EI_HYCOM, EI_aqua, 'Elephant Island')
cop_HYCOM_ts(KGI_copDF, KGI_HYCOM, KGI_aqua, 'King George Island')
cop_HYCOM_ts(CI_copDF, CI_HYCOM, CI_aqua, 'Clarence Island')



# ----------------Step 4: Statistical Testing------------------
statTest <- function(cop, hycom, aqua, site) {
  # First, matching dataframes to each other
  hycom_surf <- filter(hycom, depth == 0)
  cop <- rename(cop, cop_temp=sst, cop_sal=salinity)
  hycom_surf <- rename(hycom_surf, hycom_temp=temp, hycom_sal=salt)
  aqua <- rename(aqua,aqua_temp=sst)
  aqua$date <- as.Date(aqua$date, format="%Y-%m-%d")
  combined <- left_join(cop, hycom_surf, by='date')
  combined <- left_join(combined, aqua, by='date')
  cleaned <- combined %>% filter(!is.na(cop_temp) & !is.na(hycom_temp) & !is.na(aqua_temp))
  
  # Checking normality
  # QQ plots to visualize
  qqplots <- list(ggqqplot(combined,x='hycom_temp', title='HYCOM temp'), 
                  ggqqplot(combined,x='cop_temp', title='Copernicus temp'),
                  ggqqplot(combined,x='hycom_sal', title='HYCOM salinity'), 
                  ggqqplot(combined,x='cop_sal', title='Copernicus salinity'))
  do.call(grid.arrange, c(qqplots, nrow = 2, top = site))
  # Shapiro-Wilk Test to formally check for normality
  htemp_norm <- shapiro.test(combined$hycom_temp)
  ctemp_norm <- shapiro.test(combined$cop_temp)
  hsal_norm <- shapiro.test(combined$hycom_sal)
  csal_norm <- shapiro.test(combined$cop_sal)
  variable <- c(htemp_norm[4], ctemp_norm[4], hsal_norm[4], csal_norm[4])
  W_stat <- c(htemp_norm[1], ctemp_norm[1], hsal_norm[1], csal_norm[1])
  p_value <- c(htemp_norm[2], ctemp_norm[2], hsal_norm[2], csal_norm[2])
  norm_check <- cbind(variable, W_stat, p_value)
  print("Checking for normality of each variable: ")
  print(norm_check)
  
  # Testing to see if there is a difference between models
  # t-test (parametric)
  t_temp <- t.test(combined$hycom_temp, combined$cop_temp, paired=TRUE)
  t_sal <- t.test(combined$hycom_sal, combined$cop_sal, paired=TRUE)
  # Wilcoxon signed-rank test (non-parametric)
  wilc_temp <- wilcox.test(combined$hycom_temp, combined$cop_temp, paired=TRUE, alternative='two.sided')
  wilc_sal <- wilcox.test(combined$hycom_sal, combined$cop_sal, paired=TRUE, alternative='two.sided')
  # Making table
  test <- c('t-test','t-test','wilcoxon','wilcoxon')
  test_stat <- c(t_temp$statistic,t_sal$statistic,wilc_temp$statistic,wilc_sal$statistic)
  p_value <- c(t_temp$p.value,t_sal$p.value,wilc_temp$p.value,wilc_sal$p.value)
  diff_test <- cbind(test,test_stat,p_value)
  print("Checking to see if HYCOM and Copernicus significantly differ: ")
  print(diff_test)
  
  # Testing which matches Aqua MODIS better
  # Diebold-Mariano test to:
  #     (1) Check if both models match remotely sensed data equally well
  #     (2) Check which model matches better, and if this is statistically significant
  # For temperature
  same <- DM.test(cleaned$hycom_temp, cleaned$cop_temp, cleaned$aqua_temp, loss.type='AE',c=FALSE, H1='same')
  more <- DM.test(cleaned$hycom_temp, cleaned$cop_temp, cleaned$aqua_temp, loss.type='AE',c=FALSE, H1='more')
  less <- DM.test(cleaned$hycom_temp, cleaned$cop_temp, cleaned$aqua_temp, loss.type='AE',c=FALSE, H1='less')
  # Making table
  alternative_hypothesis <- c("different accuracy", "HYCOM more accurate", "HYCOM less accurate")
  test_stat <- c(same$statistic,more$statistic,less$statistic)
  p_value <- c(same$p.value,more$p.value,less$p.value)
  dm_test <- cbind(alternative_hypothesis, test_stat, p_value)
  print("Diebold-Mariano test to determine more accurate forecast (for temperature only): ")
  print(dm_test)

}
statTest(EI_copDF, EI_HYCOM, EI_aqua,"Elephant Island")
statTest(KGI_copDF, KGI_HYCOM, KGI_aqua,"King George Island")
statTest(CI_copDF, CI_HYCOM,CI_aqua, "Clarence Island")