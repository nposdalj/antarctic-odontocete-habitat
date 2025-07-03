library(R.matlab)
library(lubridate)
library(tidyverse)

test <- readMat("~/GitHub/antarctic-odontocete-habitat/Environmental Data/HYCOM/Data/W303E304S61N59_20140305T0000Z.mat")

# --------------------------Step 1: read files----------------------------
readHYCOM <- function(start, end) {
  start <- as.Date(start)
  end <- as.Date(end)
  dates <- seq(start, end, by = "day")
  
  df <- data.frame()

  for (current in dates) {
    # Format date into yyyymmdd
    date_str <- format(as.Date(current), "%Y%m%d")
    
    # Construct filename (you can make site-specific filename logic here if needed)
    file_path <- paste0("~/GitHub/antarctic-odontocete-habitat/Environmental Data/HYCOM/Data/W303E304S61N59_",
                        date_str, "T0000Z.mat")
    
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
    grid$date <- current
    
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

# ----------------------step 2: make timeseries-------------------
EI <- readHYCOM('2014-03-05','2014-07-17')
KGI <- readHYCOM('2014-03-05','2014-07-17') # fix dates here
CI <- readHYCOM('2014-03-05','2014-07-17')