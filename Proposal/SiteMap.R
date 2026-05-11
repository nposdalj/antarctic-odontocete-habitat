library(tidyverse)
library(raster)
library(sp)
library(sf)
library(terra)
library(ggspatial)
library(graticule)
library(ggrepel)
library(ggnewscale)

# ----------------------- OUTPUT DIRECTORY -------------------------------

output_dir <- "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Proposal/"

# ----------------------- LOAD BATHYMETRY -------------------------------

bathymetry <- raster(
  'C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/GMRT/GMRTbathymetry_highestres_widerView.grd'
)

# ----------------------- LOAD KML --------------------------------------

sites_kml <- st_read(
  "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Proposal/Antarctica.kml"
)

# Extract coordinates
coords <- st_coordinates(sites_kml)

# Build dataframe
sites <- data.frame(
  name = sites_kml$Name,
  lon = coords[,1],
  lat = coords[,2]
)

# ----------------------- SITE DEPTHS -----------------------------------

site_depths <- data.frame(
  name = c(
    "CS4","CS5","CS6","Proposed_CS",
    "Boyd_NE","Boyd_NW","Proposed_Boyd",
    "Gerlache","Ger2",
    "F1","F2","F3","Proposed BFS - TBW Flow"
  ),
  
  depth_m = c(
    500, 750, 1000, 900,
    1200, 1100, 950,
    800, 850,
    700, 750, 900, 1000
  )
)

# Merge depths
sites <- left_join(sites, site_depths, by = "name")

# Convert to sf
sites_sf <- st_as_sf(
  sites,
  coords = c("lon", "lat"),
  crs = 4326
)

# ----------------------- FUNCTION --------------------------------------

make_regional_map <- function(
    site_names,
    output_name,
    ymin_override = NULL,
    ymax_override = NULL,
    xmin_override = NULL,
    xmax_override = NULL,
    buffer_deg = 1.5,
    include_depths = TRUE,
    legend_cols = 1
){
  
  # Subset sites
  region_sites <- sites %>%
    filter(name %in% site_names)
  
  region_sites_sf <- sites_sf %>%
    filter(name %in% site_names)
  
  # Legend labels
  if(include_depths){
    
    region_sites$plot_legend <- paste0(
      region_sites$name,
      " (",
      region_sites$depth_m,
      " m)"
    )
    
  } else {
    
    region_sites$plot_legend <- region_sites$name
  }
  
  region_sites_sf$plot_legend <- region_sites$plot_legend
  
  # Dynamic extent
  xmin <- min(region_sites$lon) - buffer_deg
  xmax <- max(region_sites$lon) + buffer_deg
  ymin <- min(region_sites$lat) - buffer_deg
  ymax <- max(region_sites$lat) + buffer_deg
  
  # Override limits if provided
  if(!is.null(xmin_override)){ xmin <- xmin_override }
  if(!is.null(xmax_override)){ xmax <- xmax_override }
  if(!is.null(ymin_override)){ ymin <- ymin_override }
  if(!is.null(ymax_override)){ ymax <- ymax_override }
  
  # Crop bathymetry
  bbox <- extent(xmin, xmax, ymin, ymax)
  bathy_crop <- crop(bathymetry, bbox)
  
  # Bathy dataframe
  bathy_df <- as.data.frame(
    bathy_crop,
    xy = TRUE,
    na.rm = TRUE
  )
  
  colnames(bathy_df) <- c("x", "y", "depth")
  
  # Contours
  contours <- rasterToContour(
    bathy_crop,
    levels = c(
      -500, -1000, -1500, -2000,
      -2500, -3000, -3500, -4000
    )
  )
  
  contours_sf <- st_as_sf(contours)
  
  # ----------------------- PLOT ----------------------------------------
  
  p <- ggplot() +
    
    # Bathymetry
    geom_tile(
      data = bathy_df,
      aes(x = x, y = y, fill = depth)
    ) +
    
    scale_fill_viridis_c(
      option = "mako",
      name = "Depth (m)",
      limits = c(-5500, 0)
    ) +
    
    # Reset fill scale
    new_scale_fill() +
    
    # Land
    geom_tile(
      data = filter(bathy_df, depth >= 0),
      aes(x = x, y = y),
      fill = "#66616B"
    ) +
    
    # Contours
    geom_sf(
      data = contours_sf,
      color = "#66616B",
      size = 0.3
    ) +
    
    # Sites
    geom_sf(
      data = region_sites_sf,
      aes(fill = plot_legend),
      shape = 21,
      size = 5,
      color = "#22192d"
    ) +
    
    scale_fill_viridis_d(
      option = "plasma",
      name = "Sites",
      guide = guide_legend(ncol = legend_cols)
    ) +
    
    # Labels
    geom_text_repel(
      data = region_sites,
      aes(
        x = lon,
        y = lat,
        label = name
      ),
      color = "#22192d",
      size = 4,
      box.padding = 0.5,
      point.padding = 0.3,
      min.segment.length = 0
    ) +
    
    # Extent
    coord_sf(
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      expand = TRUE
    ) +
    
    # Scale bar
    annotation_scale(
      bar_cols = c('#F5FFF8','#22192d')
    ) +
    
    # Theme
    theme_minimal() +
    
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 11),
      legend.title = element_text(
        size = 11,
        color = '#22192d'
      ),
      legend.text = element_text(
        size = 10,
        color = '#22192d'
      ),
      legend.key.height = unit(1.2,'cm')
    )
  
  # ----------------------- SAVE ----------------------------------------
  
  ggsave(
    filename = file.path(output_dir, paste0(output_name, ".png")),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(output_name, ".pdf")),
    plot = p,
    width = 10,
    height = 8
  )
  
  return(p)
}

# =======================================================================
# FULL MAP
# =======================================================================

full_sites <- sites$name

FullMap <- make_regional_map(
  site_names = full_sites,
  output_name = "FullMap",
  ymin_override = -64,
  ymax_override = -60,
  buffer_deg = 2,
  include_depths = FALSE,
  legend_cols = 3
)

# =======================================================================
# CAPE SHIRREFF
# =======================================================================

CapeShirreff_Map <- make_regional_map(
  site_names = c(
    "CS4",
    "CS5",
    "CS6",
    "Proposed_CS"
  ),
  output_name = "CapeShirreff_Map",
  ymin_override = -63,
  ymax_override = -61.5,
  buffer_deg = 0.8,
  include_depths = FALSE,
  legend_cols = 1
)

# =======================================================================
# BOYD STRAIT
# =======================================================================

# =======================================================================
# BOYD STRAIT
# =======================================================================

BoydStrait_Map <- make_regional_map(
  site_names = c(
    "Boyd_NE",
    "Boyd_NW",
    "Proposed_Boyd"
  ),
  output_name = "BoydStrait_Map",
  
  ymin_override = -63.8,
  ymax_override = -62.4,
  
  xmin_override = -63.5,
  xmax_override = -60.5,
  
  include_depths = FALSE,
  legend_cols = 1
)

# =======================================================================
# GERLACHE STRAIT
# =======================================================================

Gerlache_Map <- make_regional_map(
  site_names = c(
    "Gerlache",
    "Ger2"
  ),
  output_name = "Gerlache_Map",
  buffer_deg = 0.4,
  include_depths = FALSE,
  legend_cols = 1
)

# =======================================================================
# BRANSFIELD STRAIT
# =======================================================================

Bransfield_Map <- make_regional_map(
  site_names = c(
    "F1",
    "F2",
    "F3",
    "Proposed BFS - TBW Flow"
  ),
  output_name = "Bransfield_Map",
  ymin_override = -64,
  ymax_override = -62,
  buffer_deg = 1,
  include_depths = FALSE,
  legend_cols = 1
)