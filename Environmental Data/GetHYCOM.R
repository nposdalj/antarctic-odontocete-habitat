library(R.matlab)
library(lubridate)
library(tidyverse)
library(gridExtra)

# --------------------------Step 1: read files----------------------------
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
      file_path <- paste0("/Users/trisha/R/antarctic-odontocete-habitat/Environmental Data/HYCOM/Data/W303E304S61N59_",
                          date_str, "T0000Z.mat")
    } else if(site == 'KGI') {
      file_path <- paste0("/Users/trisha/R/antarctic-odontocete-habitat/Environmental Data/HYCOM/Data/W303E304S62N61_",
                          date_str, "T0000Z.mat")
    } else {
      file_path <- paste0("/Users/trisha/R/antarctic-odontocete-habitat/Environmental Data/HYCOM/Data/W306E307S62N61_",
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
  #df$date <- as.POSIXct(df$date, origin = "1970-01-01", tz="UTC") # Changing from seconds to time

  return(df)
  
}
EI_df <- readHYCOM('EI','2014-03-05','2014-07-17')
KGI_df <- readHYCOM('KGI','2015-02-10','2016-01-29') 
CI_df <- readHYCOM('CI','2016-02-04','2016-12-02')


# ----------------------Step 2: Make Timeseries-------------------
plotHYCOM <- function(data, depth_layer, site) {
  filtered <- filter(data, depth == depth_layer)
  filtered <- na.omit(filtered)
  data$date <- as.Date(data$date)
  salinity <- ggplot(data = filtered, mapping = aes(x = date, y = salt)) + geom_line(color = "saddlebrown", size = 0.7) + 
     labs(x = NULL, y = "Salinity (psu)") +
    scale_x_date(date_labels = "%b %Y")
  temp <- ggplot(data = filtered, mapping = aes(x = date, y = temp)) + geom_line(color = "steelblue", size = 0.7) + 
    labs(x = NULL, y = "Temperature (°C)") +
    scale_x_date(date_labels = "%b %Y")
  
  #quartz()
  timeseries <- grid.arrange(salinity, temp, nrow = 2, top = paste(site," at ", as.character(depth_layer), "m Depth"))
  #timeseries <- grid.arrange(temp, nrow = 3, top = paste(site," at ", as.character(depth_layer), "m Depth"))
}
plotHYCOM(EI_df, 0, "Elephant Island")
plotHYCOM(KGI_df, 0, "King George Island")
plotHYCOM(CI_df, 0, "Clarence Island")