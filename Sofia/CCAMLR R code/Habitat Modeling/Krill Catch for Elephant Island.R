
library(dplyr)
library(ggplot2)
library(sf)
library(marmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ggspatial)
library(units)

EI_CI <- readRDS("D:/CCAMLR Data/756/R/756_USA_2026-01-02.Rds")
C1 <- EI_CI[["C1"]]        
OBS_HAUL <- EI_CI[["OBS_HAUL"]]  


krill_data <- OBS_HAUL %>%
  filter(!is.na(greenweight_catch_kg)) %>%   
  left_join(
    C1 %>% select(obs_haul_id, latitude_haul_start, longitude_haul_start),
    by = "obs_haul_id"
  ) %>%
  filter(!is.na(latitude_haul_start), !is.na(longitude_haul_start))

krill_sf <- st_as_sf(
  krill_data,
  coords = c("longitude_haul_start", "latitude_haul_start"),
  crs = 4326
)


antarctica_land <- ne_countries(scale = "medium", returnclass = "sf")


elephant_island <- st_sfc(st_point(c(-55.3, -61.1)), crs = 4326)


plot_area <- st_buffer(st_transform(elephant_island, 3031), dist = set_units(70000, m)) %>%
  st_transform(4326)


antarctica_land_crop <- st_crop(antarctica_land, plot_area)


krill_ocean <- krill_sf[!st_intersects(krill_sf, antarctica_land_crop, sparse = FALSE), ]


grid <- st_make_grid(
  krill_ocean,
  cellsize = 0.03,  
  what = "polygons",
  square = TRUE
)
grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid)

krill_grid <- st_join(krill_ocean, grid_sf, join = st_within)
krill_avg <- krill_grid %>%
  group_by(grid_id) %>%
  summarise(mean_catch = mean(greenweight_catch_kg, na.rm = TRUE), .groups = "drop") %>%
  st_centroid()


bathy <- getNOAA.bathy(
  lon1 = st_bbox(plot_area)["xmin"] - 0.05,
  lon2 = st_bbox(plot_area)["xmax"] + 0.05,
  lat1 = st_bbox(plot_area)["ymin"] - 0.05,
  lat2 = st_bbox(plot_area)["ymax"] + 0.05,
  resolution = 1
)
bathy_df <- fortify.bathy(bathy)


islands <- data.frame(name = "Elephant Island", lon = -55.3, lat = -61.1)
islands_sf <- st_as_sf(islands, coords = c("lon", "lat"), crs = 4326)


ggplot() +
  geom_raster(data = bathy_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(low = "#cce6ff", high = "#003366", name = "Depth (m)") +
  
  geom_sf(data = antarctica_land_crop, fill = "grey95", color = "grey70") +
  
  geom_sf(data = krill_avg, aes(size = mean_catch, color = mean_catch), alpha = 0.8) +
  scale_size_continuous(range = c(2, 6), name = "Mean Catch (kg)") +
  scale_color_viridis_c(option = "plasma", name = "Mean Catch (kg)") +
  
  geom_sf(data = islands_sf, color = "red", size = 4) +
  geom_text(
    data = islands,
    aes(x = lon, y = lat, label = name),
    nudge_y = 0.02,  
    fontface = "bold",
    color = "red",
    size = 4
  ) +
  
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  
  coord_sf(
    xlim = c(st_bbox(plot_area)["xmin"], st_bbox(plot_area)["xmax"]),
    ylim = c(st_bbox(plot_area)["ymin"], st_bbox(plot_area)["ymax"])
  ) +
  
  labs(
    title = "Krill Catch around Elephant Island",
    subtitle = "CCAMLR Site 756 – Grid-Averaged Krill Catch",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid = element_line(color = "transparent")
  )

#---------------------------------------------------
# RAW DATA
#--------------------------------------------------


library(dplyr)
library(ggplot2)
library(sf)
library(marmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ggspatial)
library(units)

EI_CI    <- readRDS("D:/CCAMLR Data/756/R/756_USA_2026-01-02.Rds")
C1       <- EI_CI[["C1"]]
OBS_HAUL <- EI_CI[["OBS_HAUL"]]

krill_data <- OBS_HAUL %>%
  filter(!is.na(greenweight_catch_kg)) %>%
  left_join(
    C1 %>% select(obs_haul_id, latitude_haul_start, longitude_haul_start),
    by = "obs_haul_id"
  ) %>%
  filter(!is.na(latitude_haul_start), !is.na(longitude_haul_start))

krill_sf <- st_as_sf(
  krill_data,
  coords = c("longitude_haul_start", "latitude_haul_start"),
  crs = 4326
)

antarctica_land <- ne_countries(scale = "medium", returnclass = "sf")

elephant_island <- st_sfc(st_point(c(-55.3, -61.1)), crs = 4326)

plot_area <- st_buffer(st_transform(elephant_island, 3031), dist = set_units(70000, m)) %>%
  st_transform(4326)

antarctica_land_crop <- st_crop(antarctica_land, plot_area)


krill_ocean <- krill_sf[!st_intersects(krill_sf, antarctica_land_crop, sparse = FALSE), ]

bathy <- getNOAA.bathy(
  lon1 = st_bbox(plot_area)["xmin"] - 0.05,
  lon2 = st_bbox(plot_area)["xmax"] + 0.05,
  lat1 = st_bbox(plot_area)["ymin"] - 0.05,
  lat2 = st_bbox(plot_area)["ymax"] + 0.05,
  resolution = 1
)
bathy_df <- fortify.bathy(bathy)

islands    <- data.frame(name = "Elephant Island", lon = -55.3, lat = -61.1)
islands_sf <- st_as_sf(islands, coords = c("lon", "lat"), crs = 4326)

ggplot() +
  geom_raster(data = bathy_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(low = "#cce6ff", high = "#003366", name = "Depth (m)") +
  
  geom_sf(data = antarctica_land_crop, fill = "grey95", color = "grey70") +
  
  geom_sf(data = krill_ocean, aes(size = greenweight_catch_kg,
                                  color = greenweight_catch_kg), alpha = 0.7) +
  scale_size_continuous(range = c(1, 6), name = "Catch (kg)") +
  scale_color_viridis_c(option = "plasma", name = "Catch (kg)") +
  
  geom_sf(data = islands_sf, color = "red", size = 4) +
  geom_text(
    data = islands,
    aes(x = lon, y = lat, label = name),
    nudge_y = 0.02,
    fontface = "bold",
    color = "red",
    size = 4
  ) +
  
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  
  coord_sf(
    xlim = c(st_bbox(plot_area)["xmin"], st_bbox(plot_area)["xmax"]),
    ylim = c(st_bbox(plot_area)["ymin"], st_bbox(plot_area)["ymax"])
  ) +
  
  labs(
    title = "Krill Catch around Elephant Island",
    subtitle = "CCAMLR Site 756 – Individual Haul Krill Catch",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(size = 16, face = "bold"),
    plot.subtitle    = element_text(size = 12),
    panel.grid       = element_line(color = "transparent")
  )
