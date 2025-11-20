library(CopernicusMarine)
library(stars)
library(ncdf4)
library(raster)
library(tidyverse)
library(lubridate)
library(gridExtra)

# ----------------- 1. Load Copernicus data -------------------
APEI <- nc_open("L:/Shared drives/Antarctic Marine Mammals/Krill Data/cmems_mod_glo_phy_my_0.083deg_P1D-m_1763418042285.nc")
APDE <- nc_open("L:/Shared drives/Antarctic Marine Mammals/Krill Data/cmems_mod_glo_phy_my_0.083deg_P1D-m_1763440909642.nc")
output_dir <- "L:/Shared drives/Antarctic Marine Mammals/Krill Data/Sofie's Analysis"

# ----------------- 2. Load krill data -------------------
data_dir   <- "L:/Shared drives/Antarctic Marine Mammals/Krill Data/CCAMLR Statistical Bulletin"
krill_file <- file.path(data_dir, "AggregatedKrillCatch.csv")

krill_raw <- readr::read_csv(krill_file, guess_max = 200000)

# ----------------- 3. Helper functions -------------------
fix_grid <- function(var) {
  d <- length(dim(var))
  if (is.null(d)) return(as.vector(var))
  if (d == 3) var <- aperm(var, c(1, 2, 3))
  if (d == 2) var <- aperm(var, c(1, 2))
  as.vector(var)
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

physFromNC <- function(data) {
  # Dimensions
  time <- ncvar_get(data, "time")
  lat  <- ncvar_get(data, "latitude")
  lon  <- ncvar_get(data, "longitude")
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")
  
  # Oceanographic variables
  temp        <- fix_grid(ncvar_get(data, "thetao"))
  ssh         <- fix_grid(ncvar_get(data, "zos"))
  n_velocity  <- fix_grid(ncvar_get(data, "vo"))
  e_velocity  <- fix_grid(ncvar_get(data, "uo"))
  salinity    <- fix_grid(ncvar_get(data, "so"))
  mixed_layer <- fix_grid(ncvar_get(data, "mlotst"))
  
  # Sea ice variables
  sice_conc    <- fix_grid(ncvar_get(data, "siconc"))
  sice_thick   <- fix_grid(ncvar_get(data, "sithick"))
  sice_e_veloc <- fix_grid(ncvar_get(data, "usi"))
  sice_n_veloc <- fix_grid(ncvar_get(data, "vsi"))
  
  lonlattime <- expand.grid(lon = lon, lat = lat, date = time_obs)
  
  df <- data.frame(
    lon   = lonlattime$lon,
    lat   = lonlattime$lat,
    date  = lonlattime$date,
    ssh   = ssh,
    n_velocity  = n_velocity,
    e_velocity  = e_velocity,
    salinity    = salinity,
    mixed_layer = mixed_layer,
    sice_conc    = sice_conc,
    sice_thick   = sice_thick,
    sice_e_veloc = sice_e_veloc,
    sice_n_veloc = sice_n_veloc,
    temp        = temp
  )
  
  # ----- cut at Jan 1993 and turn sea-ice NA/NaN into 0 -----
  df <- df %>%
    filter(date >= as.Date("1993-01-01")) %>%
    mutate(across(
      .cols = starts_with("sice_"),
      .fns  = ~ ifelse(is.na(.) | is.nan(.), 0, .)
    ))
  
  # Daily spatial averages
  avg_df <- df %>%
    group_by(date) %>%
    summarise(
      ssh_mean         = mean(ssh, na.rm = TRUE),
      ssh_sd           = sd(ssh, na.rm = TRUE),
      n_velocity_mean  = mean(n_velocity, na.rm = TRUE),
      n_velocity_mad   = mad(n_velocity, na.rm = TRUE),
      e_velocity_mean  = mean(e_velocity, na.rm = TRUE),
      e_velocity_mad   = mad(e_velocity, na.rm = TRUE),
      salinity_mean    = mean(salinity, na.rm = TRUE),
      salinity_sd      = sd(salinity, na.rm = TRUE),
      mixed_layer_mean = mean(mixed_layer, na.rm = TRUE),
      mixed_layer_sd   = sd(mixed_layer, na.rm = TRUE),
      sice_conc_mean   = mean(sice_conc, na.rm = TRUE),
      ice_conc_sd      = sd(sice_conc, na.rm = TRUE),
      sice_thick_mean  = mean(sice_thick, na.rm = TRUE),
      sice_e_veloc_mean = mean(sice_e_veloc, na.rm = TRUE),
      sice_n_veloc_mean = mean(sice_n_veloc, na.rm = TRUE),
      temp_mean        = mean(temp, na.rm = TRUE),
      temp_sd          = sd(temp, na.rm = TRUE),
      .groups = "drop"
    )
  
  avg_df
}

make_monthly <- function(df_daily) {
  df_daily %>%
    mutate(
      Calendar_Year = year(date),
      Month         = month(date)
    ) %>%
    group_by(Calendar_Year, Month) %>%
    summarise(
      across(
        .cols  = where(is.numeric),
        .fns   = safe_mean,
        .names = "{.col}_month_mean"
      ),
      .groups = "drop"
    )
}

# ----------------- 4. Build env time series -------------------
EI_daily  <- physFromNC(APEI)
KGI_daily <- physFromNC(APDE)

EI_monthly_all  <- make_monthly(EI_daily)
KGI_monthly_all <- make_monthly(KGI_daily)

# ----------------- 5. Clean krill data -------------------
krill <- krill_raw %>%
  rename(
    year         = Calendar_Year,
    month        = Month,
    group_ssmu   = Group_SSMU_Code,
    ssmu         = SSMU_Code,
    krill_weight = Krill_Green_Weight
  ) %>%
  mutate(
    year         = as.integer(year),
    month        = as.integer(month),
    group_ssmu   = as.character(group_ssmu),
    ssmu         = as.character(ssmu),
    krill_weight = as.numeric(krill_weight)
  )

ap_ssmus <- c("APDPE", "APDPW", "APEI", "APPA", "APW", "APBSE", "APBSW", "APE")

krill_ap <- krill %>%
  filter(group_ssmu == "AP",
         ssmu %in% ap_ssmus,
         !is.na(krill_weight)) %>%
  rename(
    Calendar_Year = year,
    Month         = month
  )

krill_APEI  <- krill_ap %>% filter(ssmu == "APEI")
krill_APDPE <- krill_ap %>% filter(ssmu == "APDPE")

# ----------------- 6. Join: keep ALL env months ≥ 1993 -------------------
# (env on the left, krill on the right; missing krill → NA)

APEI_krill_env <- EI_monthly_all %>%
  left_join(krill_APEI,
            by = c("Calendar_Year", "Month"))

APDPE_krill_env <- KGI_monthly_all %>%
  left_join(krill_APDPE,
            by = c("Calendar_Year", "Month"))

# ---- Monthly joined env + krill ----
readr::write_csv(APEI_krill_env,
                 file.path(output_dir, "APEI_krill_env_monthly.csv"))

readr::write_csv(APDPE_krill_env,
                 file.path(output_dir, "APDPE_krill_env_monthly.csv"))


# ----------------- 7. Make yearly tables -------------------
make_yearly <- function(df) {
  df %>%
    group_by(Calendar_Year) %>%
    summarise(
      # ----- Krill -----
      total_krill = sum(krill_weight, na.rm = TRUE),
      mean_krill  = mean(krill_weight, na.rm = TRUE),
      
      # ----- Sea-ice concentration (monthly mean over grid) -----
      mean_sice_conc   = mean(sice_conc_mean_month_mean, na.rm = TRUE),
      median_sice_conc = median(sice_conc_mean_month_mean, na.rm = TRUE),
      max_sice_conc    = max(sice_conc_mean_month_mean, na.rm = TRUE),
      min_sice_conc    = min(sice_conc_mean_month_mean, na.rm = TRUE),
      sd_sice_conc     = sd(sice_conc_mean_month_mean, na.rm = TRUE),
      
      # ----- Sea-ice thickness -----
      mean_sice_thick   = mean(sice_thick_mean_month_mean, na.rm = TRUE),
      median_sice_thick = median(sice_thick_mean_month_mean, na.rm = TRUE),
      max_sice_thick    = max(sice_thick_mean_month_mean, na.rm = TRUE),
      min_sice_thick    = min(sice_thick_mean_month_mean, na.rm = TRUE),
      sd_sice_thick     = sd(sice_thick_mean_month_mean, na.rm = TRUE),
      
      .groups = "drop"
    )
}

APEI_yearly  <- make_yearly(APEI_krill_env)
APDPE_yearly <- make_yearly(APDPE_krill_env)

## Add anomalies
## ---- APEI climatology & anomalies ----
APEI_clim <- APEI_yearly %>%
  summarise(
    clim_mean_sice_conc   = mean(mean_sice_conc,   na.rm = TRUE),
    clim_mean_sice_thick  = mean(mean_sice_thick,  na.rm = TRUE)
  )

APEI_yearly <- APEI_yearly %>%
  mutate(
    anom_mean_sice_conc  = mean_sice_conc  - APEI_clim$clim_mean_sice_conc,
    anom_mean_sice_thick = mean_sice_thick - APEI_clim$clim_mean_sice_thick
  )

## ---- APDPE climatology & anomalies ----
APDPE_clim <- APDPE_yearly %>%
  summarise(
    clim_mean_sice_conc   = mean(mean_sice_conc,   na.rm = TRUE),
    clim_mean_sice_thick  = mean(mean_sice_thick,  na.rm = TRUE)
  )

APDPE_yearly <- APDPE_yearly %>%
  mutate(
    anom_mean_sice_conc  = mean_sice_conc  - APDPE_clim$clim_mean_sice_conc,
    anom_mean_sice_thick = mean_sice_thick - APDPE_clim$clim_mean_sice_thick
  )

# ---- Yearly summaries (with anomalies) ----
readr::write_csv(APEI_yearly,
                 file.path(output_dir, "APEI_krill_env_yearly.csv"))

readr::write_csv(APDPE_yearly,
                 file.path(output_dir, "APDPE_krill_env_yearly.csv"))
