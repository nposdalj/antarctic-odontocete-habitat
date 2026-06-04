# ---------------------------------------
# KRILL CATCH vs EKE_MAD TIMESERIES
# ---------------------------------------

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(readxl)
library(sf)
library(units)
library(zoo)

# -----------------------------
# FILE PATHS
# -----------------------------
env_file <- file.path(
  "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Extended_PhysicsReanalysis_BiogeochemistryHindcast",
  "KGI_60km_mask_with_deseasoned.csv"
)

catch_file <- file.path(
  "D:/CCAMLR Data/756/R",
  "756_USA_2026-01-02.Rds"
)

# -----------------------------
# SETTINGS
# -----------------------------
taxon_code <- "KRI"

start_date_env   <- as.Date("2011-01-05")
start_date_catch <- as.Date("2011-01-01")

buffer_km     <- 60
antarctic_crs <- 3031

sites <- data.frame(
  name = c("KGI"),
  lon  = c(-57.941917),
  lat  = c(-61.457817)
)

# -----------------------------
# LOAD DATA
# -----------------------------
env_data  <- read.csv(env_file)
data_list <- readRDS(catch_file)

C1       <- data_list[["C1"]]
C1_CATCH <- data_list[["C1_CATCH"]]

# -----------------------------
# BUILD SPATIAL BUFFER
# -----------------------------
buffer_m <- as.numeric(set_units(buffer_km, "km") |> set_units("m"))

sites_sf <- st_as_sf(
  sites,
  coords = c("lon", "lat"),
  crs = 4326
) %>%
  st_transform(crs = antarctic_crs)

sites_buffer <- st_buffer(sites_sf, dist = buffer_m)

# -----------------------------
# HAUL LOCATIONS
# -----------------------------
C1_sf <- C1 %>%
  filter(
    !is.na(longitude_haul_start),
    !is.na(latitude_haul_start)
  ) %>%
  st_as_sf(
    coords = c("longitude_haul_start", "latitude_haul_start"),
    crs = 4326
  ) %>%
  st_transform(crs = antarctic_crs)

hauls_in_buffer <- st_join(C1_sf, sites_buffer, join = st_within)

haul_sites <- hauls_in_buffer %>%
  st_drop_geometry() %>%
  filter(!is.na(name)) %>%
  select(c1_id, Site = name)

# -----------------------------
# DAILY EKE DATA
# -----------------------------
daily_env <- env_data %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= start_date_env) %>%
  group_by(date) %>%
  summarise(
    EKE_mad = mean(EKE_mad, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# DAILY KRILL CATCH
# -----------------------------
daily_krill <- C1_CATCH %>%
  filter(taxon_code == "KRI") %>%
  left_join(haul_sites, by = "c1_id") %>%
  filter(!is.na(Site)) %>%
  left_join(
    C1 %>% select(c1_id, date_catchperiod_start),
    by = "c1_id"
  ) %>%
  mutate(date = as.Date(date_catchperiod_start)) %>%
  filter(date >= start_date_catch) %>%
  group_by(date) %>%
  summarise(
    krill_kg = sum(greenweight_caught_kg, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# MERGE DATASETS
# -----------------------------
timeseries_data <- full_join(
  daily_env,
  daily_krill,
  by = "date"
) %>%
  arrange(date)

# Replace NA krill values with 0
timeseries_data$krill_kg[is.na(timeseries_data$krill_kg)] <- 0

# -----------------------------
# OPTIONAL SMOOTHING
# -----------------------------
timeseries_data <- timeseries_data %>%
  mutate(
    EKE_mad_smooth = rollmean(EKE_mad, 30,
                              fill = NA,
                              align = "center"),
    
    krill_smooth = rollmean(krill_kg, 30,
                            fill = NA,
                            align = "center")
  )

# -----------------------------
# SCALE EKE FOR DUAL AXIS
# -----------------------------
scale_factor <- max(timeseries_data$krill_smooth, na.rm = TRUE) /
  max(timeseries_data$EKE_mad_smooth, na.rm = TRUE)

# -----------------------------
# PLOT
# -----------------------------
ggplot(timeseries_data, aes(x = date)) +
  
  # KRILL
  geom_line(
    aes(y = krill_smooth),
    color = "darkgreen",
    linewidth = 1
  ) +
  
  # EKE
  geom_line(
    aes(y = EKE_mad_smooth * scale_factor),
    color = "blue",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name = "Krill Catch (kg)",
    
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "EKE MAD"
    )
  ) +
  
  labs(
    title = "Krill Catch and EKE MAD Through Time",
    x = "Date"
  ) +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.title.y.left = element_text(color = "darkgreen"),
    axis.title.y.right = element_text(color = "blue")
  )
