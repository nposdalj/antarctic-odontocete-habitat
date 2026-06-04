library(tidyverse)
library(raster)
library(sf)
library(terra)
library(ggspatial)
library(ggrepel)
library(Cairo)

# ======================================================================
# OUTPUT DIRECTORY
# ======================================================================

output_dir <- "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Proposal/"

# ======================================================================
# LOAD BATHYMETRY
# ======================================================================

bathymetry <- raster(
  "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/GMRT/GMRTbathymetry_highestres_widerView.grd"
)

# ======================================================================
# PROPOSED SITES
# ======================================================================

sites <- data.frame(
  
  name = c(
    "CS01",
    "CS02",
    "BRS01",
    "BRS02",
    "BS01",
    "BS02",
    "GS01",
    "GS02"
  ),
  
  lon = c(
    -61.0768,
    -60.6446,
    -60.0100,
    -59.3671,
    -61.9696,
    -62.1779,
    -61.4524,
    -61.2667
  ),
  
  lat = c(
    -62.1876,
    -62.3195,
    -62.9085,
    -63.2416,
    -62.8671,
    -63.0019,
    -63.7361,
    -63.7500
  )
)

# Convert to sf object
sites_sf <- st_as_sf(
  sites,
  coords = c("lon", "lat"),
  crs = 4326
)

# ======================================================================
# MAP EXTENT
# ======================================================================

xmin <- -64
xmax <- -58

ymin <- -64.5
ymax <- -61.5

# ======================================================================
# CROP + DOWNSAMPLE BATHYMETRY
# ======================================================================

bbox <- extent(
  xmin,
  xmax,
  ymin,
  ymax
)

bathy_crop <- crop(
  bathymetry,
  bbox
)

# Downsample raster for faster plotting/export
bathy_crop <- aggregate(
  bathy_crop,
  fact = 4,
  fun = mean
)

# Convert raster to dataframe
bathy_df <- as.data.frame(
  bathy_crop,
  xy = TRUE,
  na.rm = TRUE
)

colnames(bathy_df) <- c(
  "x",
  "y",
  "depth"
)

# ======================================================================
# BATHYMETRIC CONTOURS
# ======================================================================

contours <- rasterToContour(
  bathy_crop,
  levels = c(
    -500,
    -1000,
    -2000,
    -3000,
    -4000
  )
)

contours_sf <- st_as_sf(contours)

# ======================================================================
# PLOT
# ======================================================================

p <- ggplot() +
  
  # --------------------------------------------------------------------
# Bathymetry
# --------------------------------------------------------------------

geom_tile(
  data = bathy_df,
  aes(
    x = x,
    y = y,
    fill = depth
  )
) +
  
  scale_fill_viridis_c(
    option = "mako",
    limits = c(-5000, 0),
    name = "Depth (m)"
  ) +
  
  # --------------------------------------------------------------------
# Land Mask
# --------------------------------------------------------------------

geom_tile(
  data = filter(bathy_df, depth >= 0),
  aes(
    x = x,
    y = y
  ),
  fill = "grey55"
) +
  
  # --------------------------------------------------------------------
# Bathymetric Contours
# --------------------------------------------------------------------

geom_sf(
  data = contours_sf,
  color = "grey40",
  linewidth = 0.3
) +
  
  # --------------------------------------------------------------------
# Site Locations
# --------------------------------------------------------------------

geom_sf(
  data = sites_sf,
  shape = 21,
  size = 4,
  fill = "#D81B60",
  color = "black",
  stroke = 0.8
) +
  
  # --------------------------------------------------------------------
# Site Labels
# --------------------------------------------------------------------

geom_text_repel(
  data = sites,
  aes(
    x = lon,
    y = lat,
    label = name
  ),
  
  size = 4,
  color = "black",
  
  box.padding = 0.5,
  point.padding = 0.3,
  
  min.segment.length = 0,
  
  seed = 123
) +
  
  # --------------------------------------------------------------------
# Map Extent
# --------------------------------------------------------------------

coord_sf(
  xlim = c(xmin, xmax),
  ylim = c(ymin, ymax),
  expand = FALSE
) +
  
  # --------------------------------------------------------------------
# Scale Bar
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Scale Bar
# --------------------------------------------------------------------

annotation_scale(
  location = "br",
  
  bar_cols = c("white", "black"),
  
  text_cex = 0.9,
  
  line_width = 0.8,
  
  pad_x = unit(0.25, "in"),
  pad_y = unit(0.25, "in"),
  
  style = "bar"
)
  
  # --------------------------------------------------------------------
# Theme
# --------------------------------------------------------------------

theme_minimal() +
  
  theme(
    
    axis.title = element_blank(),
    
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    
    plot.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 12
    ),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.background = element_rect(
      fill = "white",
      color = "black"
    ),
    
    panel.grid.major = element_line(
      color = "grey85",
      linewidth = 0.3
    )
  )

# ======================================================================
# SAVE FIGURES
# ======================================================================

# High-resolution PNG for proposals/manuscripts
ggsave(
  filename = file.path(
    output_dir,
    "Antarctic_Proposed_Sites_Map.png"
  ),
  
  plot = p,
  
  width = 10,
  height = 8,
  
  dpi = 600
)

# Optimized PDF
ggsave(
  filename = file.path(
    output_dir,
    "Antarctic_Proposed_Sites_Map.pdf"
  ),
  
  plot = p,
  
  width = 10,
  height = 8,
  
  dpi = 300,
  
  device = cairo_pdf
)

# ======================================================================
# DISPLAY PLOT
# ======================================================================

p