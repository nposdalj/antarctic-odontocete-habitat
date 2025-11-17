## ------------------------------------------------------------
## 0. Packages
## ------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

## ------------------------------------------------------------
## 1. Load & basic cleaning
## ------------------------------------------------------------

# Adjust path as needed
krill_path <- "L:/Shared drives/Antarctic Marine Mammals/Krill Data/krillbase_data.csv"
plot_dir <- "L:/Shared drives/Antarctic Marine Mammals/Krill Data/"

krill_raw <- readr::read_csv(krill_path, guess_max = 200000)

krill <- krill_raw %>%
  # Standardize column names a bit
  rename(
    lat = LATITUDE,
    lon = LONGITUDE,
    season = SEASON,
    krill_n = NUMBER_OF_KRILL_UNDER_1M2,
    krill_std = STANDARDISED_KRILL_UNDER_1M2
  ) %>%
  # Convert date (looks like d/m/Y, e.g. 15/01/2008)
  mutate(
    date = lubridate::dmy(DATE),
    year = year(date),
    # If any longitudes are 0–360, convert to -180–180
    lon = if_else(lon > 180, lon - 360, lon)
  ) %>%
  filter(!is.na(date))

glimpse(krill)

## ------------------------------------------------------------
## 2. Define your study sites & study-area radius
## ------------------------------------------------------------
## >>> EDIT THESE VALUES <<<
sites_df <- tibble::tibble(
  site = c("EI","KGI", "CI"),
  lon = c(-55.95400, -57.941917, -53.483433),
  lat = c(-60.8869, -61.457817, -61.251867))

radius_km <- 600   # “near” radius around EACH site

# Sites as sf
sites_sf <- sites_df %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Krill data as sf (if not already)
krill_sf <- krill %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Project to Antarctic polar stereographic for distance / buffers
krill_proj <- sf::st_transform(krill_sf, 3031)
sites_proj <- sf::st_transform(sites_sf, 3031)

# Create a buffer around each site, then union as study area
study_area <- sites_proj %>%
  sf::st_buffer(dist = radius_km * 1000) %>%   # radius_km in meters
  sf::st_union()

# Filter krill hauls that fall within the union of buffers
krill_near <- krill_proj[ sf::st_within(krill_proj, study_area, sparse = FALSE) , ]

# Quick check of counts
nrow(krill_proj)
nrow(krill_near)


## ------------------------------------------------------------
## 3. Bubble map of data density in study area (all years)
## ------------------------------------------------------------

# Summarize per location in study area
krill_near_summary <- krill_near %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(lon, lat) %>%
  dplyr::summarise(
    n_hauls    = dplyr::n(),
    mean_krill = mean(krill_std, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Get bounding box around sites (with some padding)
lon_buffer <- 9
lat_buffer <- 4
lon_range  <- range(sites_df$lon) + c(-lon_buffer, lon_buffer)
lat_range  <- range(sites_df$lat) + c(-lat_buffer, lat_buffer)

p_map_all = ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  # Optional: show study area boundary (back-transformed to 4326)
  geom_sf(data = sf::st_transform(study_area, 4326),
          fill = NA, color = "black", linetype = "dashed") +
  geom_sf(data = krill_near_summary,
          aes(size = n_hauls, colour = mean_krill),
          alpha = 0.7) +
  geom_sf(data = sites_sf,
          aes(shape = site),
          size = 3, colour = "red") +
  coord_sf(
    xlim = lon_range,
    ylim = lat_range
  ) +
  scale_size_continuous(name = "Number of hauls") +
  scale_colour_viridis_c(name = "Mean std. krill\n(under 1 m²)", na.value = "grey80") +
  labs(
    title = "Krillbase hauls in study area (all years)",
    subtitle = paste0("Within ", radius_km, " km of KGI, CI, and EI"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "krill_map_all_years.png"),
  plot = p_map_all,
  width = 10, height = 8, dpi = 300
)

## ------------------------------------------------------------
## 4. Time series of krill density / abundance (all years)
## ------------------------------------------------------------

krill_ts_all <- krill_near %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(
    n_hauls        = dplyr::n(),
    mean_krill_n   = mean(krill_n, na.rm = TRUE),
    mean_krill_std = mean(krill_std, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# Dot plot of standardized krill density
p_ts_all = ggplot(krill_ts_all, aes(x = date, y = mean_krill_std)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Standardised krill density in study area (all years)",
    subtitle = paste0("Within ", radius_km, " km of KGI, CI, and EI"),
    x = "Date",
    y = "Mean standardised krill (under 1 m² per day)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "krill_ts_all_years.png"),
  plot = p_ts_all,
  width = 10, height = 6, dpi = 300
)

## ------------------------------------------------------------
## 5. Map + time series for 2014–2016 only
## ------------------------------------------------------------

krill_near_1416 <- krill_near %>%
  dplyr::filter(year >= 2014, year <= 2016)

# Map (2014–2016)
krill_near_1416_summary <- krill_near_1416 %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(lon, lat) %>%
  dplyr::summarise(
    n_hauls    = dplyr::n(),
    mean_krill = mean(krill_std, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)

p_map_1416 = ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = sf::st_transform(study_area, 4326),
          fill = NA, color = "black", linetype = "dashed") +
  geom_sf(data = krill_near_1416_summary,
          aes(size = n_hauls, colour = mean_krill),
          alpha = 0.7) +
  geom_sf(data = sites_sf,
          aes(shape = site),
          size = 3, colour = "red") +
  coord_sf(
    xlim = lon_range,
    ylim = lat_range
  ) +
  scale_size_continuous(name = "Number of hauls") +
  scale_colour_viridis_c(name = "Mean std. krill\n(under 1 m²)", na.value = "grey80") +
  labs(
    title = "Krillbase hauls in study area (2014–2016)",
    subtitle = paste0("Within ", radius_km, " km of KGI, CI, and EI"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "krill_map_2014_2016.png"),
  plot = p_map_1416,
  width = 10, height = 8, dpi = 300
)

# Time series 2014–2016 (dot plots)
krill_ts_1416 <- krill_near_1416 %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(
    n_hauls        = dplyr::n(),
    mean_krill_n   = mean(krill_n, na.rm = TRUE),
    mean_krill_std = mean(krill_std, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

p_ts_1416 = ggplot(krill_ts_1416, aes(x = date, y = mean_krill_std)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Standardised krill density in study area (2014–2016)",
    subtitle = paste0("Within ", radius_km, " km of KGI, CI, and EI"),
    x = "Date",
    y = "Mean standardised krill (under 1 m² per day)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "krill_ts_2014_2016.png"),
  plot = p_ts_1416,
  width = 10, height = 6, dpi = 300
)
