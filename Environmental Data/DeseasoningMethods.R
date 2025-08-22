library(stars)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(tidyverse)
library(gridExtra) # for grid.arrange
library(stats)
library(forecast)
library(timeDate)
library(patchwork)

# ----------- Step 0: Load & Prep Data -----------------------
# Data has been downloaded at a 40 km bounding box for each site, +/- 5 years from start/end dates
# For CI: temperature, salinity, mixed layer will have NAs from 6/30/21 to 12/02/21 because of time range of data

# Opening NCDF files
# temperature, salinity, mixed layer depth 
# downloaded from: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/download
EI_phys <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/EI-deseasoning_cmems_mod_glo_phy_my_0.083deg_P1D-m_1755702365804.nc')
KGI_phys <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/KGI-deseasoning_cmems_mod_glo_phy_my_0.083deg_P1D-m_1755702629576.nc')
CI_phys <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/CI-deseasoning_cmems_mod_glo_phy_my_0.083deg_P1D-m_1755702734202.nc')

# chlorophyll, NPP
# downloaded from: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/download 
EI_bgc <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/EI-desasoning_cmems_mod_glo_bgc_my_0.25deg_P1D-m_1755702469441.nc')
KGI_bgc <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/KGI-deseasoning_cmems_mod_glo_bgc_my_0.25deg_P1D-m_1755702634659.nc')
CI_bgc <- nc_open('/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/CI-deseasoning_cmems_mod_glo_bgc_my_0.25deg_P1D-m_1755702768852.nc')

# Functions to extract variables and create dataframe
physFromNC <- function(data) {
  # Extracting relevant data
  # Dimensions
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  
  # Oceanographic variables
  sst <- ncvar_get(data, 'thetao') # sea surface temperature
  sss <- ncvar_get(data, 'so') # sea surface salinity
  mlayer <- ncvar_get(data, 'mlotst') # mixed layer depth
  
  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   sst = as.vector(sst), sss = as.vector(sss), mlayer = as.vector(mlayer))
  
  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(sst = mean(sst, na.rm=TRUE), sss = mean(sss, na.rm=TRUE),
                                                       mlayer = mean(mlayer, na.rm=TRUE))
  avg_df <- ungroup(avg_df)
  return(avg_df)
}
# Constructing site dataframes
EI_physdf <- physFromNC(EI_phys)
KGI_physdf <- physFromNC(KGI_phys)
CI_physdf <- physFromNC(CI_phys)

# Similar function to above, for biogeochem vars
bgcFromNC <- function(data) {
  # Extracting relevant data
  # Dimensions
  time <- ncvar_get(data, 'time')
  lat <- ncvar_get(data, 'latitude')
  lon <- ncvar_get(data, 'longitude')
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  
  # Biogeochem variables
  chla <- ncvar_get(data, 'chl') # chlorophyll-a concentration
  npp <- ncvar_get(data, 'nppv') # net primary production

  # Creating a dataframe with all variables
  lonlattime <- expand.grid(lon = lon, lat = lat, time = time_obs)
  df <- data.frame(lon = lonlattime$lon, lat = lonlattime$lat, date = lonlattime$time,
                   chla = as.vector(chla), npp = as.vector(npp))
  
  # Creating a dataframe with daily spatial averages
  avg_df <- df %>% group_by(date) %>% summarize(chla = mean(chla, na.rm=TRUE), npp = mean(npp, na.rm=TRUE))
  avg_df <- ungroup(avg_df)
  return(avg_df)
}
# Constructing site dataframes
EI_bgcdf <- bgcFromNC(EI_bgc)
KGI_bgcdf <- bgcFromNC(KGI_bgc)
CI_bgcdf <- bgcFromNC(CI_bgc)

# making final dataframes
EI_df <- left_join(EI_bgcdf, EI_physdf, by = 'date')
KGI_df <- left_join(KGI_bgcdf, KGI_physdf, by = 'date')
CI_df <- left_join(CI_bgcdf, CI_physdf, by = 'date')


# ----------- Step 1: Plot full timeseries ---------
# For each site, plot stacked timeseries of variables
makeTimeseries <- function(data, site) {
  data$date <- as.Date(data$date)

  # salinity
  sss_plot <- ggplot(data = data, mapping = aes(x = date, y = sss)) + 
    geom_line(color = 'gold', linewidth = .6) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Sea Surface Salinity (psu)') + 
    scale_x_date(labels = NULL) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),units = "line"), 
          plot.title = element_text(size=10, margin = margin(t = 0, b = 0), face = 'bold'),
          panel.background = element_rect(fill = 'white',color = 'black'),
          panel.grid.major = element_line(color = 'gray'))
  
  # mixed layer depth
  mlayer_plot <- ggplot(data = data, mapping = aes(x = date, y = mlayer)) + 
    geom_line(color = 'purple', linewidth = .6) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Mixed Layer Depth (m)') + 
    scale_x_date(labels = NULL) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),units = "line"), 
          plot.title = element_text(size=10, margin = margin(t = 0, b = 0), face = 'bold'),
          panel.background = element_rect(fill = 'white',color = 'black'),
          panel.grid.major = element_line(color = 'gray'))
  
  # temperature
  sst_plot <- ggplot(data = data, mapping = aes(x = date, y = sss)) + 
    geom_line(color = 'deeppink', linewidth = .6) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Sea Surface Temperature (°C)') + 
    scale_x_date(labels = NULL) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),units = "line"), 
          plot.title = element_text(size=10, margin = margin(t = 0, b = 0), face = 'bold'),
          panel.background = element_rect(fill = 'white',color = 'black'),
          panel.grid.major = element_line(color = 'gray'))
  
  # chlorophyll
  chla_plot <- ggplot(data = data, mapping = aes(x = date, y = chla)) + 
    geom_line(color = 'forestgreen', linewidth = .6) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Chlorophyll-a Concentration (mg/m^3)') + 
    scale_x_date(labels = NULL) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),units = "line"), 
          plot.title = element_text(size=10, margin = margin(t = 0, b = 0), face = 'bold'),
          panel.background = element_rect(fill = 'white',color = 'black'),
          panel.grid.major = element_line(color = 'gray'))
  
  # productivity
  npp_plot <- ggplot(data = data, mapping = aes(x = date, y = npp)) + 
    geom_line(color = 'royalblue', linewidth = .6) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Net Primary Production (mg/m3/day carbon)') + 
    scale_x_date(date_labels = "%b %Y") +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),units = "line"), 
          plot.title = element_text(size=10, margin = margin(t = 0, b = 0), face = 'bold'),
          panel.background = element_rect(fill = 'white',color = 'black'),
          panel.grid.major = element_line(color = 'gray'))
  
  # Arrange plots to one figure
  final_plot <- grid.arrange(sst_plot, sss_plot, mlayer_plot, chla_plot, npp_plot, nrow = 5, top = site)
  return(final_plot)
  
}
# Construct timeseries for all sites
EI_tsplot <- makeTimeseries(EI_df, "Elephant Island")
KGI_tsplot <- makeTimeseries(KGI_df, "King George Island")
CI_tsplot <- makeTimeseries(CI_df, "Clarence Island")



# ----------- Step 2: Deseasoning Methods ----------
# Will try: additive decomposition, multiplicative decomposition, STL, differencing
decomposeSite <- function(data, site) {
  if(site == 'EI') {
    start <- c(2009, 64)
   # end <- c(2019, 198)
  } else if(site == 'KGI') {
    start <- c(2010, 41)
   # end <- c(2021, 29)
  } else if(site == 'CI') {
    start <- c(2011, 35)
   # end <- c(2021, 336)
  }
  
  chla_df <- decomposeVar('chla', site, start, end, data)
  npp_df <- decomposeVar('npp', site, start, end, data)
  sst_df <- decomposeVar('sst', site, start, end, data)
  sss_df <- decomposeVar('sss', site, start, end, data)
  mlayer_df <- decomposeVar('mlayer', site, start, end, data)
  
  decomposed <- list(chla_df, npp_df, sst_df, sss_df, mlayer_df) %>%
    reduce(full_join, by = "date")
  
  return(decomposed)

}

decomposeVar <- function(varName, site, start, end, data) {
  # making a timeseries object
  var_ts <- ts(data[[varName]], start = start, frequency = 365)
  
  # function to pad vector to meet correct length (because deseasoning often truncates data)
  pad <- function(trend, target_length) {
    n <- length(trend)
    if(n == target_length) return(trend)
    
    # find where the non-NA values start in the trend vector
    na_start <- which(!is.na(trend))[1]
    na_end <- which(!is.na(trend))[length(which(!is.na(trend)))]
    # pad NAs at the beginning
    pad_left <- na_start - 1
    # pad NAs at the end
    pad_right <- target_length - na_end
    padded <- c(rep(NA, pad_left), trend[na_start:na_end], rep(NA, pad_right))
    return(padded)
  }
  
  # additive decomposition
  varDecomposed <- data.frame(date = data$date)
  additive <- decompose(var_ts, type = 'additive')
  additive_trend <- pad(additive[['trend']], length(var_ts))
  varDecomposed[[paste0(varName,'_additive')]] <- additive_trend
  
  # multiplicative decomposition
  multiplicative <- decompose(var_ts, type = 'multiplicative')
  multiplicative_trend <- pad(multiplicative[['trend']], length(var_ts))
  varDecomposed[[paste0(varName,'_multiplicative')]] <- multiplicative_trend
  
  
  # STL decomposition
  stl <- stl(var_ts, s.window = 'periodic',na.action = na.omit)
  stl_trend <- pad(stl$time.series[, "trend"],length(var_ts))
  varDecomposed[[paste0(varName,'_stl')]] <- stl_trend
  
  # keeping version without deseasoning
  varDecomposed[[paste0(varName, "_none")]] <- data[[varName]]
  
  return(varDecomposed)
}

EI_deseasoned <- decomposeSite(EI_df, 'EI')
KGI_deseasoned <- decomposeSite(KGI_df, 'KGI')
CI_deseasoned <- decomposeSite(CI_df, 'CI')

# ---------- Step 3: Plot Deseasoning Full Timeseries  -------------
# First, reformat data to long format for plotting
longFormat <- function(data) {
  # list of methods used
  methods <- c('Additive','Multiplicative','STL','None')
  
  # making chlorophyll dataframe longer
  chla <- data %>% pivot_longer(cols = c('chla_additive','chla_multiplicative','chla_stl','chla_none'), 
                                names_to = 'method', values_to = 'chla') %>% subset(select = c('date','method','chla'))
  method_vals <- c('chla_additive','chla_multiplicative','chla_stl','chla_none')
  chla$method <- replace(chla$method, chla$method %in% method_vals, methods)
  
  # making primary production dataframe longer
  npp <- data %>% pivot_longer(cols = c('npp_additive','npp_multiplicative','npp_stl','npp_none'), 
                                names_to = 'method', values_to = 'npp') %>% subset(select = c('date','method','npp'))
  method_vals <- c('npp_additive','npp_multiplicative','npp_stl','npp_none')
  npp$method <- replace(npp$method, npp$method %in% method_vals, methods)
  
  # making temperature dataframe longer
  sst <- data %>% pivot_longer(cols = c('sst_additive','sst_multiplicative','sst_stl','sst_none'), 
                                names_to = 'method', values_to = 'sst') %>% subset(select = c('date','method','sst'))
  method_vals <- c('sst_additive','sst_multiplicative','sst_stl','sst_none')
  sst$method <- replace(sst$method, sst$method %in% method_vals, methods)
  
  # making salinity dataframe longer
  sss <- data %>% pivot_longer(cols = c('sss_additive','sss_multiplicative','sss_stl','sss_none'), 
                                names_to = 'method', values_to = 'sss') %>% subset(select = c('date','method','sss'))
  method_vals <- c('sss_additive','sss_multiplicative','sss_stl','sss_none')
  sss$method <- replace(sss$method, sss$method %in% method_vals, methods)
  
  # making mixed layer dataframe longer
  mlayer <- data %>% pivot_longer(cols = c('mlayer_additive','mlayer_multiplicative', 'mlayer_stl','mlayer_none'), 
                                names_to = 'method', values_to = 'mlayer') %>% 
    subset(select = c('date','method','mlayer'))
  method_vals <- c('mlayer_additive','mlayer_multiplicative','mlayer_stl','mlayer_none')
  mlayer$method <- replace(mlayer$method, mlayer$method %in% method_vals, methods)
  
  # merging dataframes together
  long <- merge(chla, npp, by = c('date','method'))
  long <- merge(long, sst, by = c('date','method'))
  long <- merge(long, sss, by = c('date','method'))
  long <- merge(long, mlayer, by = c('date','method'))
  
  return(long)
}
EI_long <- longFormat(EI_deseasoned)
KGI_long <- longFormat(KGI_deseasoned)
CI_long <- longFormat(CI_deseasoned)

# For each site, plot stacked timeseries of variables
# theme for plots
base_theme <- theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),units = "line"), 
                    plot.title = element_text(size=10, margin = margin(t = 0, b = 0), face = 'bold'),
                    panel.background = element_rect(fill = 'white',color = 'black'),
                    panel.grid.major = element_line(color = 'gray'),
                    legend.position = 'none', legend.title = element_text(face = 'bold', size = 9),
                    legend.text = element_text(size = 8), legend.margin = margin(c(0,0,0,0)))
cols <- c("Additive" = "#00BFC4", "STL" = "#C77CFF", "Multiplicative" = "#7CAE00", "None" = "#FFE97C")

# function for plotting
deseasonTimeseries <- function(data, site, includeUntransformed) {
  data$date <- as.Date(data$date)

  if(includeUntransformed == FALSE) {
    data <- data %>% filter(method != 'None')
  }
  
  # salinity
  sss_plot <- ggplot(data = data, mapping = aes(x = date, y = sss, color = method)) + 
    geom_line(alpha=0.8) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Sea Surface Salinity (psu)') + 
    scale_y_continuous(limits = c((min(data$sss)*0.95), (max(data$sss)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') + base_theme
  
  # mixed layer depth
  mlayer_plot <- ggplot(data = data, mapping = aes(x = date, y = mlayer, color = method)) + 
    geom_line(alpha=0.8) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Mixed Layer Depth (m)') + 
    scale_y_continuous(limits = c((min(data$mlayer)*0.95), (max(data$mlayer)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') + base_theme
  
  # temperature
  sst_plot <- ggplot(data = data, mapping = aes(x = date, y = sst, color = method)) + 
    geom_line(alpha=0.8) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Sea Surface Temperature (°C)') + 
    scale_y_continuous(limits = c((min(data$sst)*0.95), (max(data$sst)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') + base_theme
  
  # chlorophyll
  chla_plot <- ggplot(data = data, mapping = aes(x = date, y = chla, color = method)) + 
    geom_line(alpha=0.8) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Chlorophyll-a Concentration (mg/m^3)') + 
    scale_y_continuous(limits = c((min(data$chla)*0.95), (max(data$chla)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') +  base_theme
  
  # productivity
  npp_plot <- ggplot(data = data, mapping = aes(x = date, y = npp, color = method)) + 
    geom_line(alpha=0.8) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Net Primary Production (mg/m3/day carbon)') + 
    scale_y_continuous(limits = c((min(data$npp)*0.95), (max(data$npp)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') + base_theme

  
  # Arrange plots to one figure
  final_plot <- (sst_plot | sss_plot | mlayer_plot | chla_plot | npp_plot) +
    plot_layout(ncol=1, nrow=5, guides="collect") +
    plot_annotation(title = site, 
                    theme = theme(plot.title = element_text(hjust=0.5))) &
    theme(legend.position="bottom") &
    scale_colour_manual(values = cols) &
    guides(color = guide_legend(title = "Method")) 
  return(final_plot)
  
}
# Construct timeseries for all sites
EI_all_tsFull <- deseasonTimeseries(EI_long, "Elephant Island",TRUE)
KGI_all_tsFull <- deseasonTimeseries(KGI_long, "King George Island",TRUE)
CI_all_tsFull <- deseasonTimeseries(CI_long, "Clarence Island",TRUE)
EI_deseasoned_tsFull <- deseasonTimeseries(EI_long, "Elephant Island",FALSE)
KGI_deseasoned_tsFull <- deseasonTimeseries(KGI_long, "King George Island",FALSE)
CI_deseasoned_tsFull <- deseasonTimeseries(CI_long, "Clarence Island",FALSE)


# --------------- Step 4: Plot Deseasoning Cropped Timeseries ----------
# same plots as previous steps, but cropped for time
deseasonCropped <- function(data, site, includeUntransformed) {
  data$date <- as.Date(data$date)

  if(includeUntransformed == FALSE) {
    data <- data %>% filter(method != 'None')
  }
  
  # settting date bounds
  if(site == 'Elephant Island') {
    d1 <- as.Date('03-05-2014', format = "%m-%d-%Y")
    d2 <- as.Date('07-17-2014', format = "%m-%d-%Y")
  } else if(site == 'King George Island') {
    d1 <- as.Date('02-10-2015', format = "%m-%d-%Y")
    d2 <- as.Date('01-29-2016', format = "%m-%d-%Y")
  } else if(site == 'Clarence Island') {
    d1 <- as.Date('02-04-2016', format = "%m-%d-%Y")
    d2 <- as.Date('12-02-2016', format = "%m-%d-%Y")
  }
  
  cropped <- data %>% filter(date > d1, date < d2)
  
  # salinity
  sss_plot <- ggplot(data = cropped, mapping = aes(x = date, y = sss, color = method)) + 
    geom_line(linewidth = 0.8, alpha = 0.8) + 
    xlab(NULL) + ylab(NULL) + labs(title = 'Sea Surface Salinity (psu)') + 
    scale_y_continuous(limits = c((min(data$sss)*0.95), (max(data$sss)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  # mixed layer depth
  mlayer_plot <- ggplot(data = cropped, mapping = aes(x = date, y = mlayer, color = method)) + 
    geom_line(linewidth = 0.8, alpha = 0.8) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Mixed Layer Depth (m)') + 
    scale_y_continuous(limits = c((min(data$mlayer)*0.95), (max(data$mlayer)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  # temperature
  sst_plot <- ggplot(data = cropped, mapping = aes(x = date, y = sst, color = method)) + 
    geom_line(linewidth = 0.8, alpha = 0.8) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Sea Surface Temperature (°C)') + 
    scale_y_continuous(limits = c((min(data$sst)*0.95), (max(data$sst)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  # chlorophyll
  chla_plot <- ggplot(data = cropped, mapping = aes(x = date, y = chla, color = method)) + 
    geom_line(linewidth = 0.8, alpha = 0.8) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Chlorophyll-a Concentration (mg/m^3)') + 
    scale_y_continuous(limits = c((min(data$chla)*0.95), (max(data$chla)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  # productivity
  npp_plot <- ggplot(data = cropped, mapping = aes(x = date, y = npp, color = method)) + 
    geom_line(linewidth = 0.8, alpha = 0.8) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Net Primary Production (mg/m3/day carbon)') + 
    scale_y_continuous(limits = c((min(data$npp)*0.95), (max(data$npp)*1.05))) +
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  
  # Arrange plots to one figure
  final_plot <- (sst_plot | sss_plot | mlayer_plot | chla_plot | npp_plot) +
    plot_layout(ncol=1, nrow=5, guides="collect") +
    plot_annotation(title = site, 
                    theme = theme(plot.title = element_text(hjust=0.5))) &
    theme(legend.position="bottom") &
    scale_colour_manual(values = cols) &
    guides(color = guide_legend(title = "Method")) 
  return(final_plot)
  
}
# Construct timeseries for all sites
EI_all_tsCrop <- deseasonCropped(EI_long, "Elephant Island",TRUE)
KGI_all_tsCrop <- deseasonCropped(KGI_long, "King George Island",TRUE)
CI_all_tsCrop <- deseasonCropped(CI_long, "Clarence Island",TRUE)

EI_deseasoned_tsCrop <- deseasonCropped(EI_long, "Elephant Island",FALSE)
KGI_deseasoned_tsCrop <- deseasonCropped(KGI_long, "King George Island",FALSE)
CI_deseasoned_tsCrop <- deseasonCropped(CI_long, "Clarence Island",FALSE)

# --------------- Step 5: Adding and Plotting Anomalies ---------------
anomalySite <- function(data) {
  
  varNames <- c('chla_none','npp_none','sst_none','sss_none','mlayer_none')
  vars <- c('chla','npp','sst','sss', 'mlayer')
  data$julian <- yday(data$date)
  
  for(v in varNames) {
    # differencing (based on average for that day of the year)
    data <- data %>% group_by(julian) %>% mutate(avg_val = mean(!!sym(v), na.rm=TRUE)) %>%
      ungroup()
    anomaly <- data[[v]] - data$avg_val
    
    idx <- match(v,varNames)
    data[[paste0(vars[idx], "_anomaly")]] <- anomaly
  }
  
  return(data)
}

EI_deseasoned <- anomalySite(EI_deseasoned)
KGI_deseasoned <- anomalySite(KGI_deseasoned)
CI_deseasoned <- anomalySite(CI_deseasoned)

anomalyFull <- function(data, site) {
  data$date <- as.Date(data$date)
  
  # temperature
  sst_plot <- ggplot(data = data, mapping = aes(x = date)) + 
    geom_line(mapping = aes(y = sst_anomaly), color = 'pink') +
    xlab(NULL) + ylab(NULL) + labs(title = 'Sea Surface Temperature Anomaly (°C)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') + base_theme
  
  # salinity
  sss_plot <- ggplot(data = data, mapping = aes(x = date)) + 
    geom_line(mapping = aes(y = sss_anomaly), color = 'khaki') +    
    xlab(NULL) + ylab(NULL) + labs(title = ' Sea Surface Salinity Anomaly (psu)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') + base_theme
  
  # mixed layer depth
  mlayer_plot <- ggplot(data = data, mapping = aes(x = date)) + 
    geom_line(mapping = aes(y = mlayer_anomaly), color = 'plum') +    
    xlab(NULL) + ylab(NULL) + labs(title = 'Mixed Layer Depth Anomaly (m)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') + base_theme
  
  # chlorophyll
  chla_plot <- ggplot(data = data, mapping = aes(x = date)) + 
    geom_line(mapping = aes(y = chla_anomaly), color = '#90BC8B') +        
    xlab(NULL) + ylab(NULL) + labs(title = 'Chlorophyll-a Concentration Anomaly (mg/m^3)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') +  base_theme
  
  # productivity
  npp_plot <- ggplot(data = data, mapping = aes(x = date, y = npp, color = method)) + 
    geom_line(mapping = aes(y = npp_anomaly), color = 'lightblue') +            
    xlab(NULL) + ylab(NULL) + labs(title = 'Net Primary Production Anomaly (mg/m3/day carbon)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '1 year') + base_theme
  
  final_plot <- grid.arrange(sst_plot, sss_plot, mlayer_plot, chla_plot, npp_plot, nrow = 5,
                             top = paste0(site,' Anomalies'))

  return(final_plot)
  
}
EI_anomalyFull <- anomalyFull(EI_deseasoned, 'Elephant Island')
KGI_anomalyFull <- anomalyFull(KGI_deseasoned, 'King George Island')
CI_anomalyFull <- anomalyFull(CI_deseasoned, 'Clarence Island')

anomalyCropped <- function(data, site) {
  data$date <- as.Date(data$date)
  
  # settting date bounds
  if(site == 'Elephant Island') {
    d1 <- as.Date('03-05-2014', format = "%m-%d-%Y")
    d2 <- as.Date('07-17-2014', format = "%m-%d-%Y")
  } else if(site == 'King George Island') {
    d1 <- as.Date('02-10-2015', format = "%m-%d-%Y")
    d2 <- as.Date('01-29-2016', format = "%m-%d-%Y")
  } else if(site == 'Clarence Island') {
    d1 <- as.Date('02-04-2016', format = "%m-%d-%Y")
    d2 <- as.Date('12-02-2016', format = "%m-%d-%Y")
  }
  
  cropped <- data %>% filter(date > d1, date < d2)
  
  # temperature
  sst_plot <- ggplot(data = cropped, mapping = aes(x = date)) + 
    geom_line(mapping = aes(y = sst_anomaly), color = 'pink', linewidth = 1) +
    xlab(NULL) + ylab(NULL) + labs(title = 'Sea Surface Temperature Anomaly (°C)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  # salinity
  sss_plot <- ggplot(data = cropped, mapping = aes(x = date)) + 
    geom_line(mapping = aes(y = sss_anomaly), color = 'khaki', linewidth = 1) +    
    xlab(NULL) + ylab(NULL) + labs(title = ' Sea Surface Salinity Anomaly (psu)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  # mixed layer depth
  mlayer_plot <- ggplot(data = cropped, mapping = aes(x = date)) + 
    geom_line(mapping = aes(y = mlayer_anomaly), color = 'plum', linewidth = 1) +    
    xlab(NULL) + ylab(NULL) + labs(title = 'Mixed Layer Depth Anomaly (m)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  # chlorophyll
  chla_plot <- ggplot(data = cropped, mapping = aes(x = date)) + 
    geom_line(mapping = aes(y = chla_anomaly), color = '#90BC8B', linewidth = 1) +        
    xlab(NULL) + ylab(NULL) + labs(title = 'Chlorophyll-a Concentration Anomaly (mg/m^3)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  # productivity
  npp_plot <- ggplot(data = cropped, mapping = aes(x = date, y = npp, color = method)) + 
    geom_line(mapping = aes(y = npp_anomaly), color = 'lightblue', linewidth = 1) +            
    xlab(NULL) + ylab(NULL) + labs(title = 'Net Primary Production Anomaly (mg/m3/day carbon)') + 
    scale_x_date(date_labels = "%b %Y", date_breaks = '2 months', limits = c(d1,d2)) + base_theme
  
  final_plot <- grid.arrange(sst_plot, sss_plot, mlayer_plot, chla_plot, npp_plot, nrow = 5,
                             top = paste0(site,' Anomalies'))
  
  return(final_plot)
  
}
EI_anomalyCrop <- anomalyCropped(EI_deseasoned, 'Elephant Island')
KGI_anomalyCrop <- anomalyCropped(KGI_deseasoned, 'King George Island')
CI_anomalyCrop <- anomalyCropped(CI_deseasoned, 'Clarence Island')

# Saving deseasoned/anomaly dataframes (to add to FinalData.R)
write.csv(EI_deseasoned, '/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/EI_deseasoned.csv')
write.csv(KGI_deseasoned, '/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/KGI_deseasoned.csv')
write.csv(CI_deseasoned, '/Users/trisha/scripps/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/CI_deseasoned.csv')

# --------------- Step 6: Test Model (BW29 at CI) ------------
# Testing deseasoned data and anomalies to see if they are significant for modeling
# Using BW29 at CI as an example

# Using STL decomposition because it agrees exactly with addditive and STL is a more robust method
#   reasoning here: 
#   https://docs.google.com/document/d/1c4EE4yHpgrsjrsNqD4hE1Xo3FPFdulGD5tkxw108lCY/edit?tab=t.0

# loading data (that already has STL deseasoned and anomalies in it)
detections <- read.csv('/Users/trisha/scripps/antarctic-odontocete-habitat/Data/allData_40km.csv')
CI_data <- detections %>% filter(Site == "CI")

# binning by ACF value (known from table)
# only binning deseasoned and anomaly data
deseasonedACF <- function(data, bin, species) {
  # Adding  bin_start date
  data$date <- as.Date(data$date, format = '%Y-%m-%d')
  data <- data %>% mutate(bin_start = floor_date(date, unit = paste(bin, 'days')))
  
  # Function for taking mean with correct syntax for a column
  mean_col <- function(col) expr(mean(!!sym(col), na.rm = TRUE))
  # Expression to summarize species presence
  species_expr <- set_names(list(expr(as.integer(any(!!sym(species) == 1)))), species)
  
  # Taking the average of environnmental variables
  # including deseasoned, anomalized, and normal version of each
  summarize_cols <- list(chla_stl = mean_col("chla_stl"), npp_stl = mean_col("npp_stl"),
                         sst_stl = mean_col("sst_stl"), sss_stl = mean_col("sss_stl"), 
                         mlayer_stl = mean_col("mlayer_stl"), chla_anomaly = mean_col("chla_anomaly"),
                         npp_anomaly = mean_col("npp_anomaly"),
                         sst_anomaly = mean_col("sst_anomaly"), sss_anomaly = mean_col("sss_anomaly"),
                         mlayer_anomaly = mean_col("mlayer_anomaly"), chla = mean_col('chla_0'),
                         npp = mean_col('productivity_0'), sst = mean_col('temperature_0'), 
                         sss = mean_col('salinity_0'), mlayer = mean_col('mixed_layer'))
  
  # Joining together dataframes and grouping to make final binned data
  all_summaries <- c(species_expr, summarize_cols)
  sp_binned <- data %>%
    group_by(bin_start, Site) %>%
    summarise(!!!all_summaries, .groups = "drop")
  
  return(sp_binned)
}
CI_binned <- deseasonedACF(CI_data,7, 'BW29')
# adding weights
CI_binned$weights <- ifelse(CI_binned$BW29 == 1, 2, 1)

# TESTING METHODS FOR EACH VARIABLE
# function to plot correctly
plotGams <- function(gam1,gam2,gam3, title) {
  par(mfrow = c(3, 1), mar = c(4,4,1,1), oma = c(1,1,2,1))
  plot(gam1, trans = plogis, shift = coef(gam1)[1], seWithMean = TRUE)
  plot(gam2, trans = plogis, shift = coef(gam2)[1], seWithMean = TRUE)
  plot(gam3, trans = plogis, shift = coef(gam3)[1], seWithMean = TRUE)
  mtext(title, outer = TRUE)
}

# CHLOROPHYLL
# deseasoned
gam1 <- gam(BW29 ~ s(chla_stl,k=4), weights = weights, family = binomial, data = CI_binned)
# Formula:
#   BW29 ~ s(chla_stl, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1371     0.2519   0.544    0.586
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(chla_stl) 2.15  2.544  5.286  0.0769 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0701   Deviance explained = 9.28%
# UBRE = 0.8405  Scale est. = 1         n = 51

# anomaly
gam2 <- gam(BW29 ~ s(chla_anomaly,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 86.06412
# Formula:
#   BW29 ~ s(chla_anomaly, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1885     0.2847   0.662    0.508
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(chla_anomaly) 2.891   2.99   9.87  0.0223 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.181   Deviance explained = 18.9%
# UBRE = 0.68753  Scale est. = 1         n = 51

# unchanged
gam3 <- gam(BW29 ~ s(chla,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 63.24564
# Formula:
#   BW29 ~ s(chla, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.3315     0.3741  -0.886    0.375
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(chla) 2.411  2.733  20.24   1e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.476   Deviance explained = 41.5%
# UBRE = 0.24011  Scale est. = 1         n = 51

# visualize & conclude
plotGams(gam1, gam2, gam3, 'BW29 at CI - Chlorophyll')
# only anomalized and unchanged chlorophyll are significant, but unchanged has a much lower AIC
# both have relatively interpretable trends
# because ot the much lower AIC, prefer unchanged version of chlorophyll for this model




# PRIMARY PRODUCTION
# deseasoned
gam1 <- gam(BW29 ~ s(npp_stl,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 90.16775
# Formula:
#   BW29 ~ s(npp_stl, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.09123    0.25896   0.352    0.725
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(npp_stl) 2.154  2.535  9.907  0.0181 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.126   Deviance explained = 13.1%
# UBRE =  0.768  Scale est. = 1         n = 51

# anomaly
gam2 <- gam(BW29 ~ s(npp_anomaly,k=4), weights = weights, family = binomial, data = CI_binned)
# Formula:
#   BW29 ~ s(npp_anomaly, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1637     0.2408   0.679    0.497
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(npp_anomaly) 1.534  1.865  0.776   0.611
# 
# R-sq.(adj) =  -0.00879   Deviance explained = 1.66%
# UBRE = 0.96059  Scale est. = 1         n = 51

# unchanged
gam3 <- gam(BW29 ~ s(npp,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 72.54428
# Formula:
#   BW29 ~ s(npp, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.3414     0.3747  -0.911    0.362
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(npp) 2.709  2.938  23.17 7.02e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.465   Deviance explained = 42.9%
# UBRE = 0.22636  Scale est. = 1         n = 51

# visualize & conclude
plotGams(gam1, gam2, gam3, 'BW29 at CI - Primary Production')
# only significant one is unchanged version, use that





# SURFACE TEMPERATURE
# deseasoned
gam1 <- gam(BW29 ~ s(sst_stl,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 73.88487
# Formula:
#   BW29 ~ s(sst_stl, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2395     0.3363  -0.712    0.476
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(sst_stl) 2.807  2.971  18.34 0.00037 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.36   Deviance explained = 31.2%
# UBRE = 0.45088  Scale est. = 1         n = 51

# anomaly
gam2 <- gam(BW29 ~ s(sst_anomaly,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 86.9496
# Formula:
#   BW29 ~ s(sst_anomaly, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.02258    0.27369   0.083    0.934
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(sst_anomaly) 2.844  2.982  12.17 0.00468 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.206   Deviance explained = 17.9%
# UBRE = 0.70489  Scale est. = 1         n = 51

# unchanged
gam3 <- gam(BW29 ~ s(sst,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 57.02501
# Formula:
#   BW29 ~ s(sst, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1287     0.3666  -0.351    0.726
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(sst) 2.393  2.735  22.58 3.29e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.503   Deviance explained =   48%
# UBRE = 0.11814  Scale est. = 1         n = 51

# visualize & conclude
plotGams(gam1, gam2, gam3, 'BW29 at CI - Surface Temperature')
# all three are significant, in order of lowest to highest AIC:
#   unchanged, deseasoned, anomalized
# prefer unchanged because much lower AIC, but clear trends in deseasoned and anomalized too





# SURFACE SALINITY
# deseasoned
gam1 <- gam(BW29 ~ s(sss_stl,k=4), weights = weights, family = binomial, data = CI_binned)
# Formula:
#   BW29 ~ s(sss_stl, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1659     0.2552    0.65    0.516
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(sss_stl) 1.834  2.211  5.411  0.0756 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0666   Deviance explained = 9.36%
# UBRE = 0.82663  Scale est. = 1         n = 51

# anomaly
gam2 <- gam(BW29 ~ s(sss_anomaly,k=4), weights = weights, family = binomial, data = CI_binned)
# Formula:
#   BW29 ~ s(sss_anomaly, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1901     0.2736   0.695    0.487
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(sss_anomaly) 2.79  2.968  7.753  0.0685 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.107   Deviance explained = 14.8%
# UBRE = 0.76193  Scale est. = 1         n = 51

# unchanged
gam3 <- gam(BW29 ~ s(sss,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 87.03119
# Formula:
#   BW29 ~ s(sss, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) 0.008715   0.270321   0.032    0.974
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(sss) 2.301  2.689  11.37 0.00677 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.173   Deviance explained = 16.7%
# UBRE = 0.70649  Scale est. = 1         n = 51

# visualize & conclude
plotGams(gam1, gam2, gam3, 'BW29 at CI - Surface Salinity')
# only significant one is unchanged version, use that




# MIXED LAYER DEPTH
# deseasoned
gam1 <- gam(BW29 ~ s(mlayer_stl,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 89.81315
# Formula:
#   BW29 ~ s(mlayer_stl, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1563     0.2637   0.593    0.553
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(mlayer_stl) 2.121  2.515  7.953  0.0371 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.112   Deviance explained = 13.4%
# UBRE = 0.76104  Scale est. = 1         n = 51

# anomaly
gam2 <- gam(BW29 ~ s(mlayer_anomaly,k=4), weights = weights, family = binomial, data = CI_binned)
# Formula:
#   BW29 ~ s(mlayer_anomaly, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1463     0.2455   0.596    0.551
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(mlayer_anomaly)   1      1  2.771   0.096 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0219   Deviance explained = 3.05%
# UBRE = 0.91331  Scale est. = 1         n = 51

# unchanged
gam3 <- gam(BW29 ~ s(mlayer,k=4), weights = weights, family = binomial, data = CI_binned)
# AIC: 95.53677
# Formula:
# BW29 ~ s(mlayer, k = 4)
# 
# Parametric coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1242     0.2500   0.497    0.619
# 
# Approximate significance of smooth terms:
#           edf Ref.df Chi.sq p-value  
# s(mlayer)   1      1  4.493   0.034 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0505   Deviance explained = 5.17%
# UBRE = 0.87327  Scale est. = 1         n = 51

# visualize & conclude
plotGams(gam1, gam2, gam3, 'BW29 at CI - Mixed Layer Depth')
# plots for unchanged and anomalized look exactly the same
# unchanged and deseasoned are significant, but neither look very significant
# prefer deseasoned: AIC slightly lower & it might have more of a trend (but not very large in either case)
