# ===================== LIBRARIES =====================
library(ncdf4)
library(tidyverse)
library(geosphere)
library(lubridate)
library(mgcv)
library(forecast)
library(patchwork)

# ===================== USER INPUT =====================
phys_file <- "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Extended_PhysicsReanalysis_BiogeochemistryHindcast/cmems_mod_glo_phy_my_0.083deg_P1D-m_1776348264028.nc"
bio_file  <- "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Extended_PhysicsReanalysis_BiogeochemistryHindcast/cmems_mod_glo_bgc_my_0.25deg_P1D-m_1776348387969.nc"

kgi_lat <- -61.457817
kgi_lon <- -57.941917

radii <- c(40, 60, 100)

output_dir <- dirname(phys_file)

# ===================== FUNCTIONS =====================

fix_grid <- function(var) {
  var <- aperm(var, c(1, 2, 3))
  return(as.vector(var))
}

mean_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

mad_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mad(x, na.rm = TRUE)
}

clean_var_name <- function(v) {
  gsub("_mean$", "", v)
}

# ------------------ PHYSICS ------------------
physFromNC <- function(data) {
  
  time <- ncvar_get(data, "time")
  lat  <- ncvar_get(data, "latitude")
  lon  <- ncvar_get(data, "longitude")
  
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")
  
  df <- expand.grid(lon = lon, lat = lat, date = time_obs)
  
  df$temp        <- fix_grid(ncvar_get(data, "thetao"))
  df$zos         <- fix_grid(ncvar_get(data, "zos"))
  df$u           <- fix_grid(ncvar_get(data, "uo"))
  df$v           <- fix_grid(ncvar_get(data, "vo"))
  df$salinity    <- fix_grid(ncvar_get(data, "so"))
  df$mixed_layer <- fix_grid(ncvar_get(data, "mlotst"))
  df$siconc      <- fix_grid(ncvar_get(data, "siconc"))
  df$sithick     <- fix_grid(ncvar_get(data, "sithick"))
  
  return(df)
}

# ------------------ BIO ------------------
bioFromNC <- function(data) {
  
  time <- ncvar_get(data, "time")
  lat  <- ncvar_get(data, "latitude")
  lon  <- ncvar_get(data, "longitude")
  
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")
  
  df <- expand.grid(lon = lon, lat = lat, date = time_obs)
  
  df$chla         <- fix_grid(ncvar_get(data, "chl"))
  df$o2           <- fix_grid(ncvar_get(data, "o2"))
  df$productivity <- fix_grid(ncvar_get(data, "nppv"))
  
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

# ------------------ SMOOTHED DAILY CLIMATOLOGY ANOMALIES ------------------
add_smoothed_daily_anomalies <- function(df, vars, k_smooth = 30) {
  
  df <- df %>%
    arrange(date) %>%
    mutate(
      date = as.Date(date)
    ) %>%
    filter(!(month(date)==2 & mday(date)==29)) %>%   # 👈 ADD HERE
    mutate(
      yday = yday(date)
    )
  
  for (v in vars) {
    
    base_name <- clean_var_name(v)
    anom_name <- paste0(base_name, "_anom")
    clim_name <- paste0(base_name, "_clim")
    
    temp_df <- df %>%
      select(yday, value = all_of(v)) %>%
      filter(!is.na(value))
    
    if (nrow(temp_df) > 30 && length(unique(temp_df$yday)) > 20) {
      
      # STEP 1: raw climatology (mean across years)
      clim <- temp_df %>%
        group_by(yday) %>%
        summarize(clim = mean(value, na.rm=TRUE), .groups="drop")   # 👈 YOUR LINE
      
      # STEP 2: smooth climatology
      clim_fit <- gam(
        clim ~ s(yday, bs="cc", k = k_smooth),
        data = clim,
        knots = list(yday = c(0.5, 365.5)),
        method = "REML"
      )
      
      # STEP 3: predict smooth climatology
      df[[clim_name]] <- predict(
        clim_fit,
        newdata = data.frame(yday = df$yday)
      )
      
      # STEP 4: anomaly
      df[[anom_name]] <- df[[v]] - df[[clim_name]]
      
    } else {
      df[[clim_name]] <- NA_real_
      df[[anom_name]] <- NA_real_
    }
  }
  
  return(df)
}

# ------------------ STL DE-SEASONING ------------------
add_stl_deseasoned <- function(df, vars) {
  
  df <- df %>%
    arrange(date)
  
  for (v in vars) {
    
    base_name <- clean_var_name(v)
    stl_name <- paste0(base_name, "_STL")
    
    x <- df[[v]]
    
    if (sum(!is.na(x)) > 730) {
      
      x_fill <- forecast::na.interp(x)
      
      stl_fit <- stl(
        ts(x_fill, frequency = 365),
        s.window = "periodic"
      )
      
      seasonal <- stl_fit$time.series[, "seasonal"]
      
      df[[stl_name]] <- x - seasonal
      
    } else {
      
      df[[stl_name]] <- NA_real_
    }
  }
  
  return(df)
}

# ------------------ PLOTTING ------------------
plot_variable_set <- function(df, vars, plot_title) {
  
  plot_list <- list()
  
  for (v in vars) {
    
    base_name <- clean_var_name(v)
    anom_name <- paste0(base_name, "_anom")
    stl_name  <- paste0(base_name, "_STL")
    
    p_raw <- ggplot(df, aes(date)) +
      geom_line(
        aes(y=.data[[v]]),
        linewidth=.3,
        color="black"
      ) +
      geom_line(
        aes(y=.data[[paste0(base_name,"_clim")]]),
        color="red",
        linewidth=.8
      ) +
      labs(y=v) +
      theme_bw()
    
    p_anom <- ggplot(df, aes(x = date, y = .data[[anom_name]])) +
      geom_line(linewidth = 0.3) +
      labs(y = anom_name, x = NULL) +
      theme_bw()
    
    p_stl <- ggplot(df, aes(x = date, y = .data[[stl_name]])) +
      geom_line(linewidth = 0.3) +
      labs(y = stl_name, x = "Date") +
      theme_bw()
    
    plot_list[[v]] <- p_raw / p_anom / p_stl
  }
  
  wrap_plots(plot_list, ncol = 1) +
    plot_annotation(title = plot_title)
}

# ===================== LOAD DATA =====================
phys_nc <- nc_open(phys_file)
bio_nc  <- nc_open(bio_file)

phys_df <- physFromNC(phys_nc)
bio_df  <- bioFromNC(bio_nc)

nc_close(phys_nc)
nc_close(bio_nc)

phys_df <- add_distance(phys_df)
bio_df  <- add_distance(bio_df)

# ===================== VARIABLES TO DE-SEASON =====================
vars_to_deseason <- c(
  "chla_mean",
  "o2_mean",
  "productivity_mean",
  "temp_mean",
  "mixed_layer_mean",
  "salinity_mean",
  "siconc_mean",
  "sithick_mean",
  "zos_mean",
  "EKE_mean",
  "EKE_mad"
)

bio_vars <- c(
  "chla_mean",
  "o2_mean",
  "productivity_mean"
)

phys_vars <- c(
  "temp_mean",
  "mixed_layer_mean",
  "salinity_mean",
  "siconc_mean",
  "sithick_mean"
)

dynamic_vars <- c(
  "zos_mean",
  "EKE_mean",
  "EKE_mad"
)

# ===================== MASK LOOP =====================
for (r in radii) {
  
  cat("Processing", r, "km mask...\n")
  
  # Apply mask
  phys_mask <- phys_df %>%
    filter(distance_km <= r)
  
  bio_mask <- bio_df %>%
    filter(distance_km <= r)
  
  # ----------- COMPUTE EKE AT GRID LEVEL -----------
  phys_mask <- phys_mask %>%
    mutate(EKE = 0.5 * (u^2 + v^2) * 10000)
  
  # ----------- SPATIAL AVERAGES: PHYSICS -----------
  phys_avg <- phys_mask %>%
    group_by(date) %>%
    summarize(
      temp_mean        = mean_na(temp),
      zos_mean         = mean_na(zos),
      salinity_mean    = mean_na(salinity),
      mixed_layer_mean = mean_na(mixed_layer),
      siconc_mean      = mean_na(siconc),
      sithick_mean     = mean_na(sithick),
      u_mean           = mean_na(u),
      v_mean           = mean_na(v),
      EKE_mean         = mean_na(EKE),
      EKE_mad          = mad_na(EKE),
      .groups = "drop"
    )
  
  # ----------- SPATIAL AVERAGES: BIOGEOCHEMISTRY -----------
  bio_avg <- bio_mask %>%
    group_by(date) %>%
    summarize(
      chla_mean         = mean_na(chla),
      o2_mean           = mean_na(o2),
      productivity_mean = mean_na(productivity),
      .groups = "drop"
    )
  
  # ----------- COMBINE -----------
  combined <- left_join(bio_avg, phys_avg, by = "date") %>%
    arrange(date)
  
  # ----------- DE-SEASON -----------
  combined <- add_smoothed_daily_anomalies(
    df = combined,
    vars = vars_to_deseason,
    k_smooth = 20
  )
  
  combined <- add_stl_deseasoned(
    df = combined,
    vars = vars_to_deseason
  )
  
  # ----------- SAVE TABLE -----------
  write.csv(
    combined,
    file.path(output_dir, paste0("KGI_", r, "km_mask_with_deseasoned.csv")),
    row.names = FALSE
  )
  
  # ----------- PLOTS -----------
  p_bio <- plot_variable_set(
    combined,
    bio_vars,
    paste0("KGI ", r, " km mask: Biogeochemistry")
  )
  
  p_phys <- plot_variable_set(
    combined,
    phys_vars,
    paste0("KGI ", r, " km mask: Hydrography and Sea Ice")
  )
  
  p_dynamic <- plot_variable_set(
    combined,
    dynamic_vars,
    paste0("KGI ", r, " km mask: Dynamic Variables")
  )
  
  # ----------- SAVE PLOTS -----------
  ggsave(
    file.path(output_dir, paste0("KGI_", r, "km_biogeochem_raw_anom_STL.png")),
    p_bio,
    width = 12,
    height = 14,
    units = "in",
    dpi = 300
  )
  
  ggsave(
    file.path(output_dir, paste0("KGI_", r, "km_physics_ice_raw_anom_STL.png")),
    p_phys,
    width = 12,
    height = 20,
    units = "in",
    dpi = 300
  )
  
  ggsave(
    file.path(output_dir, paste0("KGI_", r, "km_dynamic_raw_anom_STL.png")),
    p_dynamic,
    width = 12,
    height = 14,
    units = "in",
    dpi = 300
  )
}

cat("DONE\n")