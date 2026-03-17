library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(readxl)
library(sf)
library(units)

# ------------------------------
# LOAD DATA
# ------------------------------
data <- read.csv("C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning/KGI_deseasoned.csv")

#data_list <- readRDS("D:/CCAMLR Data/756/R/756_USA_2026-01-02.Rds")
data_list <- readRDS("F:/Antarctica/756_USA_2026-01-02.Rds")
C1        <- data_list[["C1"]]
C1_CATCH  <- data_list[["C1_CATCH"]]

# ------------------------------
# SPATIAL BUFFERS FOR KRILL CATCH
# ------------------------------
buffer_km     <- 100
buffer_m      <- set_units(buffer_km, "km") |> set_units("m")
antarctic_crs <- 3031

sites <- data.frame(
  name = c("KGI"),
  lon  = c(-57.941917),
  lat  = c(-61.457817)
)

sites_sf     <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326) %>%
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

# ------------------------------
# DAILY ENV DATA
# ------------------------------
daily_env_raw <- data %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2011-01-05")) %>%
  group_by(date) %>%
  summarise(
    temperature_0    = mean(sst_none,  na.rm = TRUE),
    salinity_0       = mean(sss_none,     na.rm = TRUE),
    chla_0           = mean(chla_none,         na.rm = TRUE),
    productivity_0   = mean(npp_none, na.rm = TRUE),
    mlayer   = mean(mlayer_none, na.rm = TRUE),
    recording_effort = n(),
    .groups = "drop"
  )

date_range  <- seq(min(daily_env_raw$date), max(daily_env_raw$date), by = "day")

daily_env_full <- expand.grid(
  date = date_range,
  stringsAsFactors = FALSE
) %>%
  mutate(date = as.Date(date)) %>%
  left_join(daily_env_raw, by = c("date")) %>%
  mutate(recording_effort = replace_na(recording_effort, 0))

# ------------------------------
# KRILL CATCH — spatially filtered to sites via buffer, aggregated by date
# ------------------------------
daily_krill <- C1_CATCH %>%
  filter(taxon_code == "KRI") %>%
  left_join(haul_sites, by = "c1_id") %>%
  filter(!is.na(Site)) %>%
  left_join(C1 %>% select(c1_id, date_catchperiod_start), by = "c1_id") %>%
  mutate(date = as.Date(date_catchperiod_start)) %>%
  filter(date >= as.Date("2011-01-01")) %>%
  group_by(date) %>%
  summarise(value = sum(greenweight_caught_kg, na.rm = TRUE), .groups = "drop")

daily_combined <- merge(daily_env_full, daily_krill, by = "date", all = TRUE)

### Linear Regression
plot(daily_combined$value,daily_combined$temperature_0)
abline(lm(daily_combined$temperature_0 ~ daily_combined$value), col = "red")

# To find p-value you'll have to fit a linear regression model
