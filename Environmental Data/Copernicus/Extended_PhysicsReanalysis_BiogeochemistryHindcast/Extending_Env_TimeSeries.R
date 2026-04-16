# ===================== LIBRARIES =====================
library(ncdf4)
library(tidyverse)
library(geosphere)

# ===================== USER INPUT =====================
phys_file <- "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Extended_PhysicsReanalysis_BiogeochemistryHindcast/cmems_mod_glo_phy_my_0.083deg_P1D-m_1776348264028.nc"
bio_file  <- "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Extended_PhysicsReanalysis_BiogeochemistryHindcast/cmems_mod_glo_bgc_my_0.25deg_P1D-m_1776348387969.nc"

kgi_lat <- -61.457817
kgi_lon <- -57.941917

radii <- c(40, 60, 100)

# ===================== FUNCTIONS =====================

fix_grid <- function(var) {
  var <- aperm(var, c(1, 2, 3))
  return(as.vector(var))
}

# ------------------ PHYSICS ------------------
physFromNC <- function(data) {
  
  time <- ncvar_get(data, 'time')
  lat  <- ncvar_get(data, 'latitude')
  lon  <- ncvar_get(data, 'longitude')
  
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC")
  
  df <- expand.grid(lon = lon, lat = lat, date = time_obs)
  
  df$temp  <- fix_grid(ncvar_get(data, 'thetao'))
  df$zos   <- fix_grid(ncvar_get(data, 'zos'))
  df$u     <- fix_grid(ncvar_get(data, 'uo'))
  df$v     <- fix_grid(ncvar_get(data, 'vo'))
  df$salinity <- fix_grid(ncvar_get(data, 'so'))
  df$mixed_layer <- fix_grid(ncvar_get(data, 'mlotst'))
  df$siconc  <- fix_grid(ncvar_get(data, 'siconc'))
  df$sithick <- fix_grid(ncvar_get(data, 'sithick'))
  
  return(df)
}

# ------------------ BIO ------------------
bioFromNC <- function(data) {
  
  time <- ncvar_get(data, 'time')
  lat  <- ncvar_get(data, 'latitude')
  lon  <- ncvar_get(data, 'longitude')
  
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="UTC")
  
  df <- expand.grid(lon = lon, lat = lat, date = time_obs)
  
  df$chla <- fix_grid(ncvar_get(data, 'chl'))
  df$o2   <- fix_grid(ncvar_get(data, 'o2'))
  df$productivity <- fix_grid(ncvar_get(data, 'nppv'))
  
  return(df)
}

# ------------------ DISTANCE ------------------
add_distance <- function(df) {
  df$distance_km <- distHaversine(
    cbind(df$lon, df$lat),
    c(kgi_lon, kgi_lat)
  ) / 1000
  return(df)
}

# ===================== LOAD DATA =====================
phys_nc <- nc_open(phys_file)
bio_nc  <- nc_open(bio_file)

phys_df <- physFromNC(phys_nc)
bio_df  <- bioFromNC(bio_nc)

phys_df <- add_distance(phys_df)
bio_df  <- add_distance(bio_df)

# ===================== OUTPUT DIR =====================
output_dir <- dirname(phys_file)

# ===================== MASK LOOP =====================
for (r in radii) {
  
  cat("Processing", r, "km mask...\n")
  
  # Apply mask
  phys_mask <- phys_df %>% filter(distance_km <= r)
  bio_mask  <- bio_df  %>% filter(distance_km <= r)
  
  # ----------- COMPUTE EKE AT GRID LEVEL -----------
  phys_mask <- phys_mask %>%
    mutate(EKE = 0.5 * (u^2 + v^2) * 10000)
  
  # ----------- SPATIAL AVERAGES (PHYSICS) -----------
  phys_avg <- phys_mask %>%
    group_by(date) %>%
    summarize(
      temp_mean = mean(temp, na.rm = TRUE),
      zos_mean  = mean(zos, na.rm = TRUE),
      salinity_mean = mean(salinity, na.rm = TRUE),
      mixed_layer_mean = mean(mixed_layer, na.rm = TRUE),
      siconc_mean = mean(siconc, na.rm = TRUE),
      sithick_mean = mean(sithick, na.rm = TRUE),
      u_mean = mean(u, na.rm = TRUE),
      v_mean = mean(v, na.rm = TRUE),
      EKE_mean = mean(EKE, na.rm = TRUE),
      EKE_mad  = mad(EKE, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ----------- SPATIAL AVERAGES (BIO) -----------
  bio_avg <- bio_mask %>%
    group_by(date) %>%
    summarize(
      chla_mean = mean(chla, na.rm = TRUE),
      o2_mean   = mean(o2, na.rm = TRUE),
      productivity_mean = mean(productivity, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ----------- COMBINE -----------
  combined <- left_join(bio_avg, phys_avg, by = "date")
  
  # ----------- SAVE -----------
  write.csv(
    combined,
    file.path(output_dir, paste0("KGI_", r, "km_mask.csv")),
    row.names = FALSE
  )
}

cat("DONE\n")