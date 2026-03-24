library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(readxl)
library(sf)
library(units)

# -------------------------
# CONFIGURATON
# -------------------------
env_file   <- file.path("data", "environment", "CI_deseasoned.csv")
catch_file <- file.path("data", "catch", "756_USA.Rds")

env_vars <- c("temperature_0", "salinity_0", "chla_0", "productivity_0", "mlayer_0")

taxon_code <- "KRI"

start_date_env   <- as.Date("2011-01-05")
start_date_catch <- as.Date("2011-01-01")

buffer_km     <- 100
antarctic_crs <- 3031

sites <- data.frame(
  name = "KGI",
  lon  = -57.941917,
  lat  = -61.457817
)

#-------------------------------------------------------------------
# temperature_0, salinity_0 , chla_0, productivity_0, mlayer_0, SI_0
# ------------------------------------------------------------------

env_labels <- c(
  temperature_0  = "Temperature",
  salinity_0     = "Salinity",
  chla_0         = "Chlorophyll",
  productivity_0 = "Productivity",
  mlayer_0       = "Mixed Layer"
)


data      <- read.csv(env_file)
data_list <- readRDS(catch_file)

C1       <- data_list[["C1"]]
C1_CATCH <- data_list[["C1_CATCH"]]


buffer_m <- set_units(buffer_km, "km") |> set_units("m")

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

daily_env_raw <- data %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= start_date_env) %>%
  group_by(date) %>%
  summarise(
    temperature_0  = mean(sst_none,  na.rm = TRUE),
    salinity_0     = mean(sss_none,  na.rm = TRUE),
    chla_0         = mean(chla_none, na.rm = TRUE),
    productivity_0 = mean(npp_none,  na.rm = TRUE),
    mlayer_0       = mean(mlayer_none, na.rm = TRUE),
    recording_effort = n(),
    .groups = "drop"
  )

date_range <- seq(min(daily_env_raw$date), max(daily_env_raw$date), by = "day")

daily_env_full <- expand.grid(
  date = date_range,
  stringsAsFactors = FALSE
) %>%
  mutate(date = as.Date(date)) %>%
  left_join(daily_env_raw, by = "date") %>%
  mutate(recording_effort = replace_na(recording_effort, 0))

# -------------------------
# KRILL DATA
# -------------------------
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


plot_data <- daily_combined %>%
  pivot_longer(
    cols = all_of(env_vars),
    names_to = "env_var",
    values_to = "env_value"
  ) %>%
  filter(!is.na(value), !is.na(env_value))

plot_data$label <- env_labels[plot_data$env_var]



ggplot(plot_data, aes(x = value, y = env_value)) +
  
  geom_point(
    size = 1.2,       
    color = "black",
    alpha = 0.7       
  ) +
  
  geom_smooth(
    method = "lm",
    color = "black",
    linewidth = 1.2,
    se = TRUE,
    fill = "grey80",
    alpha = 0.3,
    fullrange = TRUE
  ) +
  
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, NA),
    breaks = scales::pretty_breaks(n = 5)  
  ) +
  
  facet_wrap(~ label, scales = "free_y") +
  
  labs(
    title = "Krill Catch versus Environmental Variables",
    x = "Krill Catch (kg)",
    y = NULL
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text  = element_text(color = "black", size = 10, angle = 45, hjust = 1), 
    axis.title = element_text(face = "bold"),
    
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.3),
    
    axis.line  = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    
    strip.text = element_text(face = "bold", size = 12)  
  )

