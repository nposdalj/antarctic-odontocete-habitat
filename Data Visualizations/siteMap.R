library(tidyverse)
library(raster)
library(sp)
library(sf)
library(terra)
library(ggspatial)
library(graticule)
# add a 40 kmsquare around sites
# ----------------------- Method 1: WGS84 Projection Map  --------------------------------
# site map with default map projection
bathymetry <- raster('/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/GMRT/GMRTbathymetry_lowerres_wideView.grd')

# crop bathymetry grid to bounding box
#bbox <- extent(-59, -52, -62, -60) #original
bbox <- extent(-61, -51, -63.5, -60)
bathymetry <- crop(bathymetry, bbox)

# site coordinates
sites <- data.frame(
  name = c("EI", "KGI", "CI"),
  lon = c(-55.95400, -57.941917, -53.483433),
  lat = c(-60.8869, -61.457817, -61.251867))
sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

# 500 m contour lines
contours <- rasterToContour(bathymetry, levels = c(-500, -1000, -1500, -2000, -2500, -3000, -3500,
                                                   -4000, -4500))
contours_sf <- st_as_sf(contours)


# dataframe for bathymetry data
bathy_df <- as.data.frame(bathymetry, xy = TRUE, na.rm = TRUE)
colnames(bathy_df) <- c("x", "y", "depth")
#border_box <- st_as_sfc(st_bbox(c(xmin = -59, xmax = -52, ymin = -62, ymax = -60), crs = 4326)) #original
border_box <- st_as_sfc(st_bbox(c(xmin = -61, xmax = -51, ymin = -63.5, ymax = -60), crs = 4326))


# create map
ggplot() +
  # ocean bathymetry
  geom_tile(data = bathy_df, aes(x = x, y = y, fill = depth)) +
  # making land gray
  geom_tile(data = filter(bathy_df, depth >= 0), aes(x = x, y = y), fill = "#66616B") +
  # contour lines
  geom_sf(data = contours_sf, color = "#66616B", size = 0.3, linetype = "solid") +
  # Sites
  geom_sf(data = sites_sf, shape = 21, fill = "darkmagenta", size = 5, color = "#22192d") +
  geom_text(data = sites, aes(x = lon, y = lat, label = name),
            nudge_y = -0.12, nudge_x = 0.18, color = "#22192d", size = 8) +
  
  # Border around bathymetry extent
  geom_sf(data = border_box, fill = NA, color = "#22192d", size = 1) +
  
  # Color scale for depth
  scale_fill_viridis_c(option = "mako", name = "Depth (m)", limits = c(-5500, 0)) +
  
  # # Site points - original
  # coord_sf(xlim = c(-59,-52),
  #          ylim = c(-62,-60), expand = FALSE) 
# Site points
coord_sf(xlim = c(-61,-51),
         ylim = c(-63.5,-60), expand = FALSE) 
  
  # # Oceanographic data bounding boxes
  # # EI
  # annotate('rect', xmin = -56.69255, xmax = -55.21545, ymin = -61.24778, ymax = -60.52602, 
  #               color = 'darkmagenta', linewidth = 1, fill = 'darkmagenta', alpha = 0) +
  # # KGI
  # annotate('rect', xmin = -58.69396, xmax = -57.18988, ymin = -61.8187, ymax = -61.09694, 
  #               color = 'darkmagenta', linewidth = 1, fill = 'darkmagenta', alpha = 0) +
  # # CI
  # annotate('rect', xmin = -54.23054, xmax = -52.73632, ymin = -61.61275, ymax = -60.89099, 
  #               color = 'darkmagenta', linewidth = 1, fill = 'darkmagenta', alpha = 0) +
  # 
  # # Theming
  # theme_minimal() + 
  # annotation_scale(bar_cols = c('#F5FFF8','#22192d')) +
  # theme(axis.title = element_blank(), axis.text = element_text(size = 11,vjust = -0.2),
  #       legend.title = element_text(size=11,color='#22192d'),
  #       legend.text = element_text(size=10,color='#22192d'),
  #       legend.key.height = unit(1.5,'cm'))

# ------------- Method 2: Lambert Azimuthal Equal-Area Projection Map -------------
# NOT USING THIS CODE/METHOD - IGNORE BELOW
# LAEA projection to minimize polar size distortion

# reproject bathymetry data
laea_proj <- "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +units=m +datum=WGS84"
bathymetry_laea <- projectRaster(bathymetry, crs = laea_proj, method = "bilinear")

# LAEA dataframe
laea_df <- as.data.frame(bathymetry_laea, xy = TRUE, na.rm = TRUE)
colnames(laea_df) <- c("x", "y", "depth")

# reproject site points to LAEA
sites_laea <- st_transform(sites_sf, crs = laea_proj)

# replot in new projection
ggplot() +
  geom_tile(data = laea_df, aes(x = x, y = y, fill = depth)) +
  geom_tile(data = subset(laea_df, depth >= 0), aes(x = x, y = y), fill = "gray40") +
  geom_sf(data = sites_laea, shape = 21, fill = "darkmagenta", color = "black", size = 3) +
  geom_text(data = st_coordinates(sites_laea) %>% as.data.frame() %>% mutate(name = sites$name),
            aes(X, Y, label = name), nudge_y = -25000, color = "black", size = 5) +
  scale_fill_viridis_c(option = "mako", name = "Depth (m)") +
  coord_sf(crs = laea_proj, expand = FALSE) +
  theme_minimal() +
  theme(axis.title = element_blank())