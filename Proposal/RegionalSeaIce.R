library(tidyverse)
library(sf)
library(terra)
library(ncdf4)
library(geosphere)
library(patchwork)

# =======================================================================
# OUTPUT DIRECTORY
# =======================================================================

output_dir <- "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Proposal/"

# =======================================================================
# LOAD KML
# =======================================================================

sites_kml <- st_read(
  "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Proposal/Antarctica.kml"
)

coords <- st_coordinates(sites_kml)

sites <- data.frame(
  name = sites_kml$Name,
  lon = coords[,1],
  lat = coords[,2]
)

# =======================================================================
# DEFINE REGIONS
# =======================================================================

regions <- list(
  
  CapeShirreff = c(
    "CS4",
    "CS5",
    "CS6",
    "Proposed_CS"
  ),
  
  BoydStrait = c(
    "Boyd_NE",
    "Boyd_NW",
    "Proposed_Boyd"
  ),
  
  GerlacheStrait = c(
    "Gerlache",
    "Ger2"
  ),
  
  BransfieldStrait = c(
    "F1",
    "F2",
    "F3",
    "Proposed BFS - TBW Flow"
  )
)

# =======================================================================
# LOAD SEA ICE FILES
# =======================================================================

nc_files <- list(
  CapeShirreff = "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Proposal/CS.nc",
  BoydStrait = "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Proposal/BS.nc",
  GerlacheStrait = "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Proposal/GS.nc",
  BransfieldStrait = "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Proposal/BrS.nc"
)

# =======================================================================
# FUNCTION TO EXTRACT SEA ICE
# =======================================================================

extract_seaice_timeseries <- function(nc_path, region_sites){
  
  nc <- nc_open(nc_path)
  
  # Variables
  lon <- ncvar_get(nc, "longitude")
  lat <- ncvar_get(nc, "latitude")
  time <- ncvar_get(nc, "time")
  
  print(range(lon))
  print(range(lat))
  
  # Convert time
  
  time_units <- ncatt_get(nc, "time", "units")$value
  
  print(time_units)
  
  # Extract origin date
  origin_string <- sub(".*since ", "", time_units)
  
  # Remove timezone if present
  origin_string <- strsplit(origin_string, " ")[[1]][1]
  
  # Convert based on units
  if(grepl("days since", time_units)){
    
    dates <- as.Date(time, origin = origin_string)
    
  } else if(grepl("hours since", time_units)){
    
    dates <- as.POSIXct(time * 3600,
                        origin = paste0(origin_string, " 00:00:00"),
                        tz = "UTC")
    
    dates <- as.Date(dates)
    
  } else if(grepl("seconds since", time_units)){
    
    dates <- as.POSIXct(time,
                        origin = paste0(origin_string, " 00:00:00"),
                        tz = "UTC")
    
    dates <- as.Date(dates)
  }
  
  # Last 3 years only
  keep_dates <- dates >= (max(dates) - 365*3)
  
  dates <- dates[keep_dates]
  
  # Sea ice variables
  siconc <- ncvar_get(nc, "siconc")[,,keep_dates]
  sithick <- ncvar_get(nc, "sithick")[,,keep_dates]
  
  nc_close(nc)
  
  # Grid
  grid <- expand.grid(
    lon = lon,
    lat = lat
  )
  
  # -------------------------------------------------------------------
  # Loop through sites
  # -------------------------------------------------------------------
  
  all_site_data <- list()
  
  for(i in 1:nrow(region_sites)){
    
    site_name <- region_sites$name[i]
    site_lon <- region_sites$lon[i]
    site_lat <- region_sites$lat[i]
    
    # Calculate distances
    distances <- distHaversine(
      matrix(c(grid$lon, grid$lat), ncol = 2),
      c(site_lon, site_lat)
    )
    
    # 10 km mask
    mask <- distances <= 25000
    
    # Get indices
    mask_indices <- which(mask)
    
    # Skip site if no grid cells found
    if(length(mask_indices) == 0){
      
      warning(
        paste("No grid cells found for", site_name)
      )
      
      next
    }
    
    # Convert to matrix indices
    grid_indices <- arrayInd(
      mask_indices,
      .dim = c(length(lon), length(lat))
    )
    
    # Force matrix format if only one cell found
    grid_indices <- as.matrix(grid_indices)
    
    # ----------------------------------------------------------------
    # Time series
    # ----------------------------------------------------------------
    
    ts_df <- data.frame(
      date = dates,
      site = site_name,
      concentration = NA,
      thickness = NA
    )
    
    for(t in 1:length(dates)){
      
      conc_vals <- c()
      thick_vals <- c()
      
      for(j in seq_len(nrow(grid_indices))){
        
        x <- grid_indices[j,1]
        y <- grid_indices[j,2]
        
        conc_vals <- c(
          conc_vals,
          siconc[x,y,t]
        )
        
        thick_vals <- c(
          thick_vals,
          sithick[x,y,t]
        )
      }
      
      ts_df$concentration[t] <- mean(conc_vals, na.rm = TRUE) * 100
      ts_df$thickness[t] <- mean(thick_vals, na.rm = TRUE)
    }
    
    all_site_data[[site_name]] <- ts_df
  }
  
  final_df <- bind_rows(all_site_data)
  
  return(final_df)
}

# =======================================================================
# FUNCTION TO PLOT
# =======================================================================

plot_region_seaice <- function(data, region_name){
  
  # Concentration
  p1 <- ggplot(
    data,
    aes(
      x = date,
      y = concentration,
      color = site
    )
  ) +
    geom_line(linewidth = 1) +
    labs(
      title = paste(region_name, "- Sea Ice Concentration"),
      y = "Sea Ice Concentration (%)",
      x = NULL
    ) +
    theme_bw()
  
  # Thickness
  p2 <- ggplot(
    data,
    aes(
      x = date,
      y = thickness,
      color = site
    )
  ) +
    geom_line(linewidth = 1) +
    labs(
      title = paste(region_name, "- Sea Ice Thickness"),
      y = "Sea Ice Thickness (m)",
      x = "Date"
    ) +
    theme_bw()
  
  final_plot <- p1 / p2
  
  # Save
  ggsave(
    filename = file.path(
      output_dir,
      paste0(region_name, "_SeaIceTimeseries.png")
    ),
    plot = final_plot,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(region_name, "_SeaIceTimeseries.pdf")
    ),
    plot = final_plot,
    width = 12,
    height = 8
  )
  
  return(final_plot)
}

# =======================================================================
# RUN ALL REGIONS
# =======================================================================

for(region_name in names(regions)){
  
  cat("\nProcessing:", region_name, "\n")
  
  region_sites <- sites %>%
    filter(name %in% regions[[region_name]])
  
  seaice_df <- extract_seaice_timeseries(
    nc_path = nc_files[[region_name]],
    region_sites = region_sites
  )
  
  assign(
    paste0(region_name, "_SeaIce"),
    seaice_df
  )
  
  plot_region_seaice(
    seaice_df,
    region_name
  )
}