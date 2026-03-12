library(dplyr)
library(ggplot2)
library(sf)
library(marmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ggspatial)
library(units)

# Resolve select conflict
select <- dplyr::select

# Load data
EI_CI <- readRDS("D:/CCAMLR Data/756/R/756_USA_2026-01-02.Rds")
C1       <- EI_CI[["C1"]]
C1_CATCH <- EI_CI[["C1_CATCH"]]
OBS_HAUL <- EI_CI[["OBS_HAUL"]]

# Filter krill catch data and join coordinates
krill_data <- OBS_HAUL %>%
  filter(!is.na(greenweight_catch_kg)) %>%
  left_join(
    C1 %>% dplyr::select(obs_haul_id, latitude_haul_start, longitude_haul_start),
    by = "obs_haul_id"
  ) %>%
  filter(!is.na(latitude_haul_start), !is.na(longitude_haul_start))

krill_sf <- st_as_sf(
  krill_data,
  coords = c("longitude_haul_start", "latitude_haul_start"),
  crs = 4326
)

# Crop land to data extent
antarctica_land <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_crop(
    xmin = min(krill_data$longitude_haul_start) - 1,
    xmax = max(krill_data$longitude_haul_start) + 1,
    ymin = min(krill_data$latitude_haul_start)  - 1,
    ymax = max(krill_data$latitude_haul_start)  + 1
  )

# Remove land points
krill_ocean <- krill_sf[!st_intersects(krill_sf, antarctica_land, sparse = FALSE), ]

# Islands — Clarence Island at exactly -54 lon, -61.1 lat
# nudge_y is per-island: Elephant +0.3 (up), King George +0.3 (up),
# Clarence -0.3 (down) to avoid overlap with Elephant Island
islands <- data.frame(
  name    = c("Elephant Island", "King George Island", "Clarence Island"),
  lon     = c(-55.3,  -58.0,  -54.0),
  lat     = c(-61.1,  -62.0,  -61.1),
  nudge_y = c( 0.3,    0.3,   -0.3)   # Clarence nudged DOWN
)

islands_sf <- st_as_sf(islands, coords = c("lon", "lat"), crs = 4326)

# ── Grid averaging (shared across plots) ─────────────────────────────────────
grid <- st_make_grid(krill_ocean, cellsize = 0.5, what = "polygons", square = TRUE)
grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid)

krill_grid <- st_join(krill_ocean, grid_sf, join = st_within)
krill_avg  <- krill_grid %>%
  group_by(grid_id) %>%
  summarise(mean_catch = mean(greenweight_catch_kg, na.rm = TRUE), .groups = "drop") %>%
  st_centroid()

# ── Bathymetry (shared across plots) ─────────────────────────────────────────
bathy <- getNOAA.bathy(
  lon1 = min(krill_data$longitude_haul_start) - 1,
  lon2 = max(krill_data$longitude_haul_start) + 1,
  lat1 = min(krill_data$latitude_haul_start)  - 1,
  lat2 = max(krill_data$latitude_haul_start)  + 1,
  resolution = 1
)
bathy_df <- fortify.bathy(bathy)

# ── Plot 1: Raw haul locations ────────────────────────────────────────────────
ggplot() +
  geom_sf(data = antarctica_land, fill = "grey95", color = "grey70") +
  geom_sf(data = krill_ocean,
          aes(size = greenweight_catch_kg, color = greenweight_catch_kg),
          alpha = 0.7) +
  scale_size_continuous(range = c(2, 8), name = "Catch (kg)") +
  scale_color_viridis_c(option = "plasma", name = "Catch (kg)") +
  labs(title    = "Spatial Distribution of Krill Catch",
       subtitle = "Real haul locations from CCAMLR C1 data",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 14)

# ── Plot 2: Raw haul locations with island labels ─────────────────────────────
ggplot() +
  geom_sf(data = antarctica_land, fill = "grey95", color = "grey70") +
  geom_sf(data = krill_ocean,
          aes(size = greenweight_catch_kg, color = greenweight_catch_kg),
          alpha = 0.7) +
  geom_sf(data = islands_sf, color = "red", size = 3) +
  geom_text(data = islands,
            aes(x = lon, y = lat, label = name, vjust = ifelse(nudge_y > 0, -0.5, 1.5)),
            fontface = "bold", color = "red", size = 3) +
  scale_size_continuous(range = c(2, 8), name = "Catch (kg)") +
  scale_color_viridis_c(option = "plasma", name = "Catch (kg)") +
  labs(title    = "Spatial Distribution of Krill Catch in site 756",
       subtitle = "Real haul locations from CCAMLR C1 data",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 14)

# ── Plot 3: Grid-averaged, no labels ─────────────────────────────────────────
ggplot() +
  geom_sf(data = antarctica_land, fill = "grey95", color = "grey70") +
  geom_sf(data = krill_avg, aes(size = mean_catch, color = mean_catch), alpha = 0.7) +
  scale_size_continuous(range = c(2, 8), name = "Mean Catch (kg)") +
  scale_color_viridis_c(option = "plasma", name = "Mean Catch (kg)") +
  labs(title = "Spatial Distribution of Krill Catch (Grid-Averaged)",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 14)

# ── Plot 4: Grid-averaged with island labels (black) ─────────────────────────
ggplot() +
  geom_sf(data = antarctica_land, fill = "grey95", color = "grey70") +
  geom_sf(data = krill_avg, aes(size = mean_catch, color = mean_catch), alpha = 0.7) +
  geom_sf(data = islands_sf, color = "black", size = 3) +
  geom_text(data = islands,
            aes(x = lon, y = lat, label = name, vjust = ifelse(nudge_y > 0, -0.5, 1.5)),
            fontface = "bold", color = "black", size = 3) +
  scale_size_continuous(range = c(2, 8), name = "Mean Catch (kg)") +
  scale_color_viridis_c(option = "plasma", name = "Mean Catch (kg)") +
  labs(title    = "Spatial Distribution of Krill Catch (Averaged by Grid)",
       subtitle = "Real haul locations from CCAMLR C1 data",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 14)

# ── Plot 5: Bathymetry + grid-averaged catch + island labels ──────────────────
ggplot() +
  geom_raster(data = bathy_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(low = "#cce6ff", high = "#0066cc", name = "Depth (m)") +
  
  geom_sf(data = antarctica_land, fill = "grey95", color = "grey70") +
  
  geom_sf(data = krill_avg, aes(size = mean_catch, color = mean_catch), alpha = 0.8) +
  
  geom_sf(data = islands_sf, color = "black", size = 3) +
  geom_text(data = islands,
            aes(x = lon, y = lat, label = name, vjust = ifelse(nudge_y > 0, -0.5, 1.5)),
            fontface = "bold", color = "black", size = 3) +
  
  scale_size_continuous(range = c(2, 8), name = "Mean Catch (kg)") +
  scale_color_viridis_c(option = "plasma", name = "Mean Catch (kg)") +
  
  labs(title    = "Spatial Distribution of Grid-Averaged Krill Catch",
       subtitle = "CCAMLR site 756",
       x = "Longitude", y = "Latitude") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position    = "right",
    legend.title       = element_text(size = 6),
    legend.text        = element_text(size = 5),
    legend.key.size    = unit(0.3, "cm"),
    legend.box.spacing = unit(0.1, "cm"),
    plot.title         = element_text(size = 16, face = "bold"),
    plot.subtitle      = element_text(size = 12),
    panel.grid         = element_line(color = "transparent")
  )

# ── Plot 6: Density hotspots ──────────────────────────────────────────────────
ggplot() +
  geom_sf(data = antarctica_land, fill = "grey95", color = "grey70") +
  
  stat_density_2d(
    data = krill_data,
    aes(x = longitude_haul_start, y = latitude_haul_start,
        fill = after_stat(level), alpha = after_stat(level)),
    geom    = "polygon",
    contour = TRUE
  ) +
  scale_fill_viridis_c(option = "plasma", name = "Krill Density") +
  scale_alpha(range = c(0.3, 0.7), guide = "none") +
  
  stat_density_2d(
    data = krill_data,
    aes(x = longitude_haul_start, y = latitude_haul_start),
    geom      = "contour",
    color     = "red",
    linewidth = 0.5
  ) +
  
  geom_sf(data = islands_sf, color = "black", size = 3) +
  geom_text(data = islands,
            aes(x = lon, y = lat, label = name, vjust = ifelse(nudge_y > 0, -0.5, 1.5)),
            fontface = "bold", color = "black", size = 3) +
  
  labs(title    = "Krill Density Hotspots",
       subtitle = "CCAMLR site 756",
       x = "Longitude", y = "Latitude") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    panel.grid      = element_line(color = "transparent")
  )