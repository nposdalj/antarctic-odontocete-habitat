#1 -- load environmental/krill data --

# -- visualization of data--
#   shape
#   normal distribution
#   fit according model

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(readxl)
library(sf)
library(units)
library(car)
library(readr)

# -----------------------------
# CONFIGURATIONS
# -----------------------------

KGI_60km_mask <- read_csv("~/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Extended_PhysicsReanalysis_BiogeochemistryHindcast/KGI_60km_mask.csv")
catch_file <- file.path("D:/CCAMLR Data/756/R", "756_USA_2026-01-02.Rds")

env_var    <- "productivity_mean"# chla_mean, o2_mean, productivity_mean, temp_mean, zos_mean, salinity_mean, mixed_layer_mean, siconc_mean, sithick_mean, u_mean, v_mean, EKE_mean
env_var2   <- "chla_mean"
taxon_code <- "KRI"

env_label <- env_labels[[env_var]]

y_label <- env_labels[env_var]


data      <- read.csv("~/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Extended_PhysicsReanalysis_BiogeochemistryHindcast/KGI_60km_mask.csv")
data_list <- readRDS(catch_file)

C1       <- data_list[["C1"]]
C1_CATCH <- data_list[["C1_CATCH"]]


buffer_m <- as.numeric(set_units(buffer_km, "km") |> set_units("m"))

sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = antarctic_crs)

sites_buffer <- st_buffer(sites_sf, dist = buffer_m)

C1_sf <- C1 %>%
  filter(!is.na(longitude_haul_start), !is.na(latitude_haul_start)) %>%
  st_as_sf(coords = c("longitude_haul_start", "latitude_haul_start"), crs = 4326) %>%
  st_transform(crs = antarctic_crs)

hauls_in_buffer <- st_join(C1_sf, sites_buffer, join = st_within)

haul_sites <- hauls_in_buffer %>%
  st_drop_geometry() %>%
  filter(!is.na(name)) %>%
  select(c1_id, Site = name)

# -----------------------------
# DAILY ENVIRONMENT DATA
# -----------------------------
daily_env_raw <- data %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= start_date_env) %>%
  group_by(date) %>%
  summarise(
    chla_mean  = mean(chla_none,  na.rm = TRUE),
    o2_mean = mean(o2_none,  na.rm = TRUE),
    productivity_mean = mean(productivity_none, na.rm = TRUE),
    temp_mean = mean(temp_none,  na.rm = TRUE),
    zos_mean = mean(zos_none, na.rm = TRUE),
    salinity_mean = mean(salinity_none, na.rm = TRUE),
    mixed_layer_mean = mean(mixed_layer_none, na.rm = TRUE),
    siconc_mean = mean(siconc_none, na.rm = TRUE),
    sithick_mean = mean(sithick_none, na.rm = TRUE),
    u_mean = mean(u_none, na.rm = TRUE),
    v_mean = mean(v_none, na.rm = TRUE),
    EKE_mean = mean(EKE_none, na.rm = TRUE),
    .groups = "drop"
  )

date_range <- seq(min(daily_env_raw$date), max(daily_env_raw$date), by = "day")

daily_env_full <- expand.grid(date = date_range) %>%
  mutate(date = as.Date(date)) %>%
  left_join(daily_env_raw, by = "date") %>%
  mutate(recording_effort = ifelse(is.na(recording_effort), 0, recording_effort))


daily_krill <- C1_CATCH %>%
  filter(taxon_code == taxon_code) %>%
  left_join(haul_sites, by = "c1_id") %>%
  filter(!is.na(Site)) %>%
  left_join(C1 %>% select(c1_id, date_catchperiod_start), by = "c1_id") %>%
  mutate(date = as.Date(date_catchperiod_start)) %>%
  filter(date >= start_date_catch) %>%
  group_by(date) %>%
  summarise(value = sum(greenweight_caught_kg, na.rm = TRUE), .groups = "drop")

daily_combined <- merge(daily_env_full, daily_krill, by = "date", all.x = TRUE)

model_data <- daily_combined %>%
  filter(!is.na(value), !is.na(.data[[env_var]]))
# how to reshape poison into normal dsitribution (log transform and a square root transform) 
# -------- model_data$value_sqrt - sqrt(model_data$value) / hist(model_data$value_sqrt) --------

# -----------------------------
# TRANSFORM KRILL DATA
# -----------------------------

# Square root transformation
# Visualize the original vs transformed distributions

#2 -- Check/adjust for autocorrelation --
#   Find how much autocorrelation -> adjust
#   How many days is data autocorrelated (find avg in autocorrelation bins)
# simple GAM model (need to know the shape of your data) with 'value' to the left of the formula and a few environmental variables to the right
# run function acf() on the residuals of your simple gam model 
# ACFval to average into bins


#3  -- Look for covariants (env_data)
#    Independent variables
#   Find similar variables -> remove similar

#4 -- Build GAM model ---
# side 1: dependent
#side 2: independent
# find p value: stepwise selection (take least significant and remove least significant variable until all variables are significant)

#5 -- Visualize model -- 



