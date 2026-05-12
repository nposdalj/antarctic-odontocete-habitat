# =========================================================
# REVISED EKE MAD VISUALIZATION
# =========================================================
#
# Changes:
# - Only uses year 2016
# - Bottom LEFT = EKE mean timeseries
# - Bottom RIGHT = EKE MAD timeseries
# - Highlights selected example day
# - Adds exact date to figure title
#
# =========================================================

# ===================== LIBRARIES =====================

library(ncdf4)
library(tidyverse)
library(geosphere)
library(lubridate)
library(patchwork)

# ===================== USER INPUT =====================

phys_file <- "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Extended_PhysicsReanalysis_BiogeochemistryHindcast/cmems_mod_glo_phy_my_0.083deg_P1D-m_1776348264028.nc"

output_dir <- file.path(
  dirname(phys_file),
  "EKE_MAD_visualizations"
)

dir.create(output_dir, showWarnings = FALSE)

# KGI coordinates
kgi_lat <- -61.457817
kgi_lon <- -57.941917

# Radius mask (km)
radius_km <- 60

# ===================== FUNCTIONS =====================

fix_grid <- function(var) {
  var <- aperm(var, c(1, 2, 3))
  return(as.vector(var))
}

# ===================== LOAD DATA =====================

phys_nc <- nc_open(phys_file)

time <- ncvar_get(phys_nc, "time")
lat  <- ncvar_get(phys_nc, "latitude")
lon  <- ncvar_get(phys_nc, "longitude")

time_obs <- as.POSIXct(
  time,
  origin = "1970-01-01",
  tz = "UTC"
)

phys_df <- expand.grid(
  lon  = lon,
  lat  = lat,
  date = time_obs
)

phys_df$u <- fix_grid(
  ncvar_get(phys_nc, "uo")
)

phys_df$v <- fix_grid(
  ncvar_get(phys_nc, "vo")
)

nc_close(phys_nc)

# ===================== DISTANCE MASK =====================

phys_df$distance_km <- distHaversine(
  cbind(phys_df$lon, phys_df$lat),
  c(kgi_lon, kgi_lat)
) / 1000

phys_mask <- phys_df %>%
  filter(distance_km <= radius_km)

# ===================== CALCULATE EKE =====================

phys_mask <- phys_mask %>%
  mutate(
    EKE = 0.5 * (u^2 + v^2) * 10000
  )

# ===================== DAILY METRICS =====================

daily_metrics <- phys_mask %>%
  group_by(date) %>%
  summarize(
    EKE_mean   = mean(EKE, na.rm = TRUE),
    EKE_median = median(EKE, na.rm = TRUE),
    EKE_MAD    = mad(EKE, na.rm = TRUE),
    .groups = "drop"
  )

# ===================== KEEP ONLY 2016 =====================

phys_mask <- phys_mask %>%
  filter(year(date) == 2016)

daily_metrics <- daily_metrics %>%
  filter(year(date) == 2016)

# ===================== LOW/HIGH MAD DAYS =====================

low_mad_day <- daily_metrics$date[
  which.min(daily_metrics$EKE_MAD)
]

high_mad_day <- daily_metrics$date[
  which.max(daily_metrics$EKE_MAD)
]

print(low_mad_day)
print(high_mad_day)

# ===================== PLOTTING FUNCTION =====================

plot_EKE_MAD <- function(
    df,
    metrics_df,
    selected_day,
    figure_label
) {
  
  # -----------------------------
  # Subset selected day
  # -----------------------------
  
  day_df <- df %>%
    filter(date == selected_day)
  
  # -----------------------------
  # Statistics
  # -----------------------------
  
  med_EKE <- median(
    day_df$EKE,
    na.rm = TRUE
  )
  
  MAD_EKE <- mad(
    day_df$EKE,
    na.rm = TRUE
  )
  
  day_df <- day_df %>%
    mutate(
      abs_dev = abs(EKE - med_EKE)
    )
  
  # -----------------------------
  # PANEL 1: Raw EKE field
  # -----------------------------
  
  p1 <- ggplot(
    day_df,
    aes(lon, lat, fill = EKE)
  ) +
    
    geom_raster() +
    
    coord_fixed() +
    
    scale_fill_viridis_c() +
    
    labs(
      title = "Raw EKE Field",
      subtitle = paste0(
        "Median EKE = ",
        round(med_EKE, 2)
      ),
      fill = "EKE"
    ) +
    
    theme_bw()
  
  # -----------------------------
  # PANEL 2: Absolute deviations
  # -----------------------------
  
  p2 <- ggplot(
    day_df,
    aes(lon, lat, fill = abs_dev)
  ) +
    
    geom_raster() +
    
    coord_fixed() +
    
    scale_fill_viridis_c() +
    
    labs(
      title = "Spatial Variability in EKE",
      subtitle = paste0(
        "MAD = ",
        round(MAD_EKE, 2)
      ),
      fill = "|EKE - median|"
    ) +
    
    theme_bw()
  
  # -----------------------------
  # PANEL 3: EKE mean timeseries
  # -----------------------------
  
  p3 <- ggplot(
    metrics_df,
    aes(date, EKE_mean)
  ) +
    
    geom_line(
      linewidth = 0.5
    ) +
    
    geom_point(
      data = metrics_df %>%
        filter(date == selected_day),
      color = "red",
      size = 3
    ) +
    
    labs(
      title = "Mean EKE Through Time",
      x = "Date",
      y = "Mean EKE"
    ) +
    
    theme_bw()
  
  # -----------------------------
  # PANEL 4: EKE MAD timeseries
  # -----------------------------
  
  p4 <- ggplot(
    metrics_df,
    aes(date, EKE_MAD)
  ) +
    
    geom_line(
      linewidth = 0.5
    ) +
    
    geom_point(
      data = metrics_df %>%
        filter(date == selected_day),
      color = "red",
      size = 3
    ) +
    
    labs(
      title = "Spatial MAD of EKE Through Time",
      x = "Date",
      y = "EKE MAD"
    ) +
    
    theme_bw()
  
  # -----------------------------
  # Combine
  # -----------------------------
  
  combined_plot <- (
    p1 | p2
  ) / (
    p3 | p4
  ) +
    
    plot_annotation(
      title = paste0(
        figure_label,
        "\n",
        format(selected_day, "%Y-%m-%d")
      )
    )
  
  return(combined_plot)
}

# ===================== CREATE FIGURES =====================

p_low <- plot_EKE_MAD(
  phys_mask,
  daily_metrics,
  low_mad_day,
  "LOW MAD Example"
)

p_high <- plot_EKE_MAD(
  phys_mask,
  daily_metrics,
  high_mad_day,
  "HIGH MAD Example"
)

# ===================== SAVE FIGURES =====================

ggsave(
  file.path(
    output_dir,
    "LOW_MAD_2016.png"
  ),
  p_low,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "HIGH_MAD_2016.png"
  ),
  p_high,
  width = 14,
  height = 10,
  dpi = 300
)

cat("DONE\n")