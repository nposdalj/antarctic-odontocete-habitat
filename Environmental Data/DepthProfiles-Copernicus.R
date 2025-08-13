library(tidyverse)
library(zoo)
library(oce)

# ------------- Step 0: Data prepping --------------
# Load Copernicus Data
copernicus <- read.csv("/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/copernicus_depths.csv")
copernicus$date <- as.Date(copernicus$date,format='%Y-%m-%d')

# Extend data to make each measurement a depth range
copernicus_extended <- copernicus %>% arrange(date, depth) %>% group_by(date) %>%
  mutate(depth_min = depth,
         depth_max = lead(depth, default = max(depth) + 100)) %>% ungroup()
         
# ------------- Step 1: Depth profile for temperature ----------------
# filtering out depths that temperature has NA values for, changing depth ranges based on this
temp_filtered <- copernicus_extended %>% filter(!is.na(temp_mean)) %>% 
  arrange(date, depth) %>% group_by(date) %>%
  mutate(depth_min = depth, depth_max = lead(depth, default = max(depth) + 100)) %>% ungroup()

tempProfile <- function(site) {
  # setting date bounds
  if(site == 'EI') {
    start <- as.Date('2014-03-05') 
    end <- as.Date('2014-07-17')
  } else if(site == 'KGI') {
    start <- as.Date('2015-02-10') 
    end <- as.Date('2016-01-29')
  } else {
    start <- as.Date('2016-02-04') 
    end <- as.Date('2016-12-02')
  }
  
  # filtering for site date
  filtered <- temp_filtered %>% filter(date >= start & date <= end)
  
  profile <- ggplot(data = filtered, aes(x=date,y=-depth,fill=temp_mean)) + 
    # colored rectangle for each depth range
    geom_rect(aes(ymin = -depth_min, ymax = -depth_max,
                  xmin=date-0.5, xmax=date+0.5)) +
    scale_fill_viridis_c(option = "inferno", name ='Temperature (°C)') +
    
    # axis labels, plot scales, and theme
    labs(y = 'Depth Below Sea Surface (m)', x = 'Date', title = paste0('Temperature Profile at ',site)) +
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_date(breaks = seq(start,end, by = '2 month'),
                 date_labels = "%b %Y", expand = c(0, 0)) + theme_bw()
  
  return(profile)
}

EI_temp <- tempProfile('EI')
KGI_temp <- tempProfile('KGI')
CI_temp <- tempProfile('CI')


# ------------- Step 2: Depth profile for salinity ----------------
# filtering out depths that temperature has NA values for, changing depth ranges based on this
sal_filtered <- copernicus_extended %>% filter(!is.na(salinity_mean)) %>% 
  arrange(date, depth) %>% group_by(date) %>%
  mutate(depth_min = depth, depth_max = lead(depth, default = max(depth) + 100)) %>% ungroup()


salProfile <- function(site) {
  # settting date bounds
  if(site == 'EI') {
    start <- as.Date('2014-03-05') 
    end <- as.Date('2014-07-17')
  } else if(site == 'KGI') {
    start <- as.Date('2015-02-10') 
    end <- as.Date('2016-01-29')
  } else {
    start <- as.Date('2016-02-04') 
    end <- as.Date('2016-12-02')
  }
  
  # filtering for site date
  filtered <- temp_filtered %>% filter(date >= start & date <= end)
  
  profile <- ggplot(data = filtered, aes(x=date,y=-depth,fill=salinity_mean)) + 
    # colored rectangle for each depth range
    geom_rect(aes(ymin = -depth_min, ymax = -depth_max,
                  xmin=date-0.5, xmax=date+0.5)) +
    scale_fill_viridis_c(option = "viridis", name ='Salinity (psu)') +
    
    # axis labels, plot scales, and theme
    labs(y = 'Depth Below Sea Surface (m)', x = 'Date', title = paste0('Salinity Profile at ',site)) +
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_date(breaks = seq(start,end, by = '2 month'),
                 date_labels = "%b %Y", expand = c(0, 0)) + theme_bw()
  
  return(profile)
}

EI_sal <- salProfile('EI')
KGI_sal <- salProfile('KGI')
CI_sal <- salProfile('CI')


# ------------- Step 3: Depth profile for oxygen  ----------------
o2Profile <- function(site) {
  # setting date bounds
  if(site == 'EI') {
    start <- as.Date('2014-03-05') 
    end <- as.Date('2014-07-17')
  } else if(site == 'KGI') {
    start <- as.Date('2015-02-10') 
    end <- as.Date('2016-01-29')
  } else {
    start <- as.Date('2016-02-04') 
    end <- as.Date('2016-12-02')
  }
  
  # filtering for site dates
  filtered <- copernicus_extended %>% filter(date >= start & date <= end)
  
  
  profile <- ggplot(data = filtered, aes(x=date,y=-depth,fill=o2_mean)) + 
    # colored rectangle for each depth range
    geom_rect(aes(ymin = -depth_min, ymax = -depth_max,
                  xmin=date-0.5, xmax=date+0.5)) +
    scale_fill_viridis_c(option = "viridis", name ='O2 Concentration') +
    
    # axis labels, plot scales, and theme
    labs(y = 'Depth Below Sea Surface (m)', x = 'Date', title = paste0('Oxygen Profile at ',site)) +
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_date(breaks = seq(start,end, by = '2 month'),
                 date_labels = "%b %Y", expand = c(0, 0)) + theme_bw()
  
  return(profile)
}

EI_o2 <- o2Profile('EI')
KGI_o2 <- o2Profile('KGI')
CI_o2 <- o2Profile('CI')

# ---------------- Step 4: Rough density profile --------------
# calculate density based on temperature, salinity, and depth to get a rough gauge of water masses

# first, filtering dataframe to depth bins that are populated with temperature and salinity
# setting rectangle bins based on new depth bins
density_df <- copernicus_extended %>% filter(!is.na(salinity_mean)) %>%
  arrange(date, depth) %>% group_by(date) %>%
  mutate(depth_min = depth, depth_max = lead(depth, default = max(depth) + 100)) %>% ungroup()


# to derive pressure, using formula P = ρgh where ρ = 1023.6 kg/m3, g = 9.80665 m/s2, h = depth
#   values received from here: https://bluerobotics.com/learn/pressure-depth-calculator/
# also expressing pressure in decibars, dividng value from above formula by 10^
density_df <- density_df %>% mutate(pressure = ((1023.6*9.80665*depth)/10000))

# adding density column
density_df <- density_df %>% mutate(density = swRho(salinity_mean, temperature=temp_mean, 
                                                    pressure=pressure, eos = 'unesco'))

# function for density profile
densityProfile <- function(site,data) {
  # setting date bounds
  if(site == 'EI') {
    start <- as.Date('2014-03-05') 
    end <- as.Date('2014-07-17')
  } else if(site == 'KGI') {
    start <- as.Date('2015-02-10') 
    end <- as.Date('2016-01-29')
  } else {
    start <- as.Date('2016-02-04') 
    end <- as.Date('2016-12-02')
  }
  
  # filtering for site dates
  filtered <- data %>% filter(date >= start & date <= end)
  
  
  profile <- ggplot(data = filtered, aes(x=date,y=-depth,fill=density)) + 
    # colored rectangle for each depth range
    geom_rect(aes(ymin = -depth_min, ymax = -depth_max,
                  xmin=date-0.5, xmax=date+0.5)) +
    scale_fill_viridis_c(option = "plasma", name ='Density (kg/m3)',direction = -1) +
    
    # axis labels, plot scales, and theme
    labs(y = 'Depth Below Sea Surface (m)', x = 'Date', title = paste0('Density Profile at ',site)) +
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_date(breaks = seq(start,end, by = '2 month'),
                 date_labels = "%b %Y", expand = c(0, 0)) + theme_bw()
  
  return(profile)
}

EI_density <- densityProfile('EI',density_df)
KGI_density <- densityProfile('KGI',density_df)
CI_density <- densityProfile('CI',density_df)

shallow <- density_df %>% filter(depth < 150)
EI_shallow <- densityProfile('EI',shallow)
KGI_shallow <- densityProfile('KGI',shallow)
CI_shallow <- densityProfile('CI',shallow)