# =============================================================================
# Script  : map_helpers.R
# Section : 2.1 Study area (shared helper functions)
# Purpose : Shared settings and ggplot building blocks for the study-area maps
#           (site_map.R, currents_map.R): GMRT bathymetry and land, isobaths,
#           HARP site markers, and the common map frame (axes with degree
#           labels, scale bar, north arrow, theme).
# Inputs  : Environmental Data/GMRT/GMRTbathymetry_highres_widerView.grd
#             (GMRT MapTool grid, legacy GMT/netCDF format)
#           rnaturalearth / rnaturalearthdata (inset coastlines)
#           final_scripts_odontocetes/config.R (must be sourced first: ENV_DIR,
#             SITES, SITE_ORDER, SITE_COLORS)
# Outputs : none (functions and settings only)
# Author  : Natalie Posdaljian
# Date    : 2026-11-17
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Bathymetry source citation:
#           [INSERT GMRT citation, e.g., Ryan et al. (2009) GMRT version and
#           download date]
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(ncdf4)
  library(ggspatial)
})

# -----------------------------------------------------------------------------
# 1. Shared settings
# -----------------------------------------------------------------------------
# Map extent (degrees) chosen to frame all three sites, the northern South
# Shetland shelf, Bransfield Strait, and the tip of the Antarctic Peninsula.
MAP_XLIM <- c(-61, -51)
MAP_YLIM <- c(-63.5, -60)

# Isobaths every 500 m outline the shelf break, the South Shetland Trench,
# and the Powell Basin.
ISOBATH_LEVELS_M <- seq(500, 5500, by = 500)

LAND_COLOR <- "grey55"
INK_COLOR  <- "#22192d"

DEPTH_LIMITS  <- c(0, 5500)
DEPTH_PALETTE <- viridisLite::mako(256, direction = -1)

BATHY_FILE <- file.path(ENV_DIR, "GMRT", "GMRTbathymetry_highres_widerView.grd")

# -----------------------------------------------------------------------------
# 2. Bathymetry
# -----------------------------------------------------------------------------
# GMRT MapTool .grd files use the legacy GMT netCDF layout (a flat `z` vector
# plus x_range / y_range / dimension attributes) that GDAL/terra cannot open
# directly, so the grid is rebuilt from its metadata. `node_offset = 0` means
# values sit on grid nodes, so the raster extent is padded by half a cell.
# Rows are stored north to south, which is terra's native order.
read_gmt_grid <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  x_range <- ncvar_get(nc, "x_range")
  y_range <- ncvar_get(nc, "y_range")
  spacing <- ncvar_get(nc, "spacing")
  dims    <- ncvar_get(nc, "dimension")
  z       <- as.vector(ncvar_get(nc, "z"))

  rast(
    nrows = dims[2], ncols = dims[1],
    xmin = x_range[1] - spacing[1] / 2, xmax = x_range[2] + spacing[1] / 2,
    ymin = y_range[1] - spacing[2] / 2, ymax = y_range[2] + spacing[2] / 2,
    crs = "EPSG:4326", vals = z
  )
}

# GMRT elevation (m, positive up) cropped to the map extent.
load_bathymetry <- function(xlim = MAP_XLIM, ylim = MAP_YLIM) {
  read_gmt_grid(BATHY_FILE) |>
    crop(ext(c(xlim, ylim)), snap = "out")
}

# ggplot layers for bathymetry, land, and isobaths. The grid is converted to
# a colour image (ocean: depth on the mako scale, light = shallow, dark =
# deep; land: grey) and drawn with annotation_raster, which renders a true
# raster: geom_tile/geom_raster fall back to one rectangle per cell under
# coord_sf, leaving hairline seams and producing a very heavy PDF. Two
# invisible points carry the depth range so ggplot still draws a colourbar.
# `ocean_alpha` < 1 fades the ocean colours, and `palette` can replace the
# default mako ramp (both used when vectors are drawn on top).
bathymetry_layers <- function(bathy, ocean_alpha = 1, palette = DEPTH_PALETTE) {
  depth_to_color <- scales::col_numeric(palette, domain = DEPTH_LIMITS,
                                        na.color = "transparent")
  elevation <- as.matrix(bathy, wide = TRUE)
  ocean_col <- depth_to_color(pmin(pmax(-elevation, 0), DEPTH_LIMITS[2]))
  bathy_img <- matrix(
    ifelse(elevation >= 0, LAND_COLOR, alpha(ocean_col, ocean_alpha)),
    nrow = nrow(elevation)
  )
  bathy_ext <- as.vector(ext(bathy))   # xmin, xmax, ymin, ymax

  legend_anchor <- tibble(x = mean(bathy_ext[1:2]), y = mean(bathy_ext[3:4]),
                          depth = DEPTH_LIMITS)

  # Isobaths are computed on the elevation grid (negative depths)
  isobaths <- as.contour(bathy, levels = -ISOBATH_LEVELS_M) |>
    st_as_sf()

  list(
    annotation_raster(bathy_img, xmin = bathy_ext[1], xmax = bathy_ext[2],
                      ymin = bathy_ext[3], ymax = bathy_ext[4],
                      interpolate = FALSE),
    geom_point(data = legend_anchor, aes(x = x, y = y, fill = depth),
               alpha = 0),
    scale_fill_gradientn(colours = alpha(palette, ocean_alpha),
                         limits = DEPTH_LIMITS,
                         breaks = seq(0, 5000, by = 1000),
                         name = "Depth (m)",
                         guide = guide_colourbar(order = 1)),
    geom_sf(data = isobaths, color = "grey35", linewidth = 0.15, alpha = 0.8)
  )
}

# -----------------------------------------------------------------------------
# 3. Sites
# -----------------------------------------------------------------------------
# HARP sites as white points outlined in the site palette used in all
# manuscript figures, labelled with the site code only (no depths). Labels
# sit 0.32 deg east and 0.13 deg south of each point unless `label_offset`
# (a named list of c(dx, dy) in degrees per site) moves them, e.g., away
# from other map elements.
site_layers <- function(label_offset = list()) {
  sites_df <- SITES |>
    select(site, lon, lat) |>
    mutate(site = factor(site, levels = SITE_ORDER))
  sites_sf <- st_as_sf(sites_df, coords = c("lon", "lat"), crs = 4326)

  offsets <- lapply(setNames(SITE_ORDER, SITE_ORDER),
                    function(s) label_offset[[s]] %||% c(0.32, -0.13))
  labels_df <- sites_df |>
    mutate(lon = lon + map_dbl(as.character(site), \(s) offsets[[s]][1]),
           lat = lat + map_dbl(as.character(site), \(s) offsets[[s]][2]))

  list(
    geom_sf(data = sites_sf, aes(color = site), shape = 21, fill = "white",
            size = 3.2, stroke = 1.4, show.legend = FALSE),
    scale_color_manual(values = SITE_COLORS),
    geom_label(data = labels_df, aes(x = lon, y = lat, label = site),
               size = 4, fontface = "bold",
               color = INK_COLOR, fill = alpha("white", 0.8),
               linewidth = 0, label.padding = unit(0.12, "lines"))
  )
}

# -----------------------------------------------------------------------------
# 4. Map frame
# -----------------------------------------------------------------------------
# Axis labels use an explicit Unicode degree sign (drawn with the ragg /
# cairo devices) instead of ggplot's default graticule labels, which rendered
# as a missing-glyph box on the Windows PNG device in the previous map.
degree_sign <- intToUtf8(176)
label_lon <- function(x) paste0(abs(x), degree_sign, "W")
label_lat <- function(x) paste0(formatC(abs(x), format = "f", digits = 1),
                                degree_sign, "S")

# Scale bar and north arrow sit in the open Powell Basin (bottom right,
# `corner = "br"`) unless another corner is requested. The scale
# bar is computed at its own latitude; map scale changes by ~10%
# across the 3.5 degrees of latitude shown.
map_frame <- function(xlim = MAP_XLIM, ylim = MAP_YLIM, corner = "br") {
  list(
    annotation_scale(location = corner, width_hint = 0.2,
                     bar_cols = c("white", INK_COLOR), text_col = INK_COLOR,
                     line_col = INK_COLOR, pad_y = unit(0.4, "cm")),
    annotation_north_arrow(location = corner, which_north = "true",
                           height = unit(1, "cm"), width = unit(0.8, "cm"),
                           pad_x = unit(0.6, "cm"), pad_y = unit(1.1, "cm"),
                           style = north_arrow_fancy_orienteering(
                             line_col = INK_COLOR,
                             fill = c("white", INK_COLOR),
                             text_col = INK_COLOR)),
    scale_x_continuous(breaks = seq(-60, -52, by = 2), labels = label_lon),
    scale_y_continuous(breaks = seq(-63.5, -60, by = 0.5), labels = label_lat),
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 4326,
             default_crs = 4326),
    theme_bw(base_size = 11),
    theme(axis.title = element_blank(),
          axis.text = element_text(color = INK_COLOR),
          panel.grid = element_blank(),
          legend.title = element_text(color = INK_COLOR),
          legend.text = element_text(color = INK_COLOR),
          legend.key.height = unit(1.4, "cm"))
  )
}

# -----------------------------------------------------------------------------
# 5. Locator inset
# -----------------------------------------------------------------------------
# Antarctic Peninsula / Drake Passage / southern South America in a Lambert
# azimuthal equal-area projection centred on the study area (north up), with
# the main-map extent outlined. Natural Earth 1:50m land (bundled with
# rnaturalearthdata) is sufficient at this scale. Same design as the inset
# built in site_map.R.
locator_inset <- function(xlim = MAP_XLIM, ylim = MAP_YLIM) {
  inset_crs <- "+proj=laea +lat_0=-61 +lon_0=-60 +datum=WGS84 +units=m"

  land_inset <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf") |>
    st_geometry() |>
    st_make_valid() |>
    st_crop(xmin = -90, xmax = -30, ymin = -80, ymax = -45) |>
    st_union() |>
    st_transform(inset_crs)

  # Study-area rectangle, densified so its edges follow the projection
  study_rect <- st_bbox(c(xmin = xlim[1], xmax = xlim[2],
                          ymin = ylim[1], ymax = ylim[2]), crs = 4326) |>
    st_as_sfc() |>
    st_segmentize(units::set_units(10, km)) |>
    st_transform(inset_crs)

  inset_labels <- tibble(
    label = c("Antarctic\nPeninsula", "South\nAmerica", "Drake\nPassage"),
    lon   = c(-55, -62.5, -64),
    lat   = c(-65.5, -54.6, -58.6)
  ) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    st_transform(inset_crs)

  ggplot() +
    geom_sf(data = land_inset, fill = LAND_COLOR, color = "grey30",
            linewidth = 0.1) +
    geom_sf(data = study_rect, fill = NA, color = "red3", linewidth = 0.7) +
    geom_sf_text(data = inset_labels, aes(label = label), size = 2.2,
                 color = INK_COLOR, lineheight = 0.85, fontface = "italic") +
    coord_sf(crs = inset_crs, xlim = c(-1.3e6, 1.1e6),
             ylim = c(-0.85e6, 0.85e6), expand = FALSE, datum = NA) +
    theme_void() +
    theme(panel.background = element_rect(fill = "aliceblue", color = NA),
          panel.border = element_rect(fill = NA, color = INK_COLOR,
                                      linewidth = 0.5))
}

# PNG (600 dpi, ragg) and PDF (cairo, for Unicode glyphs) of one figure.
save_map <- function(plot, file_stem, width_in, height_in) {
  ggsave(out_path("01_study_area", paste0(file_stem, ".png")), plot,
         width = width_in, height = height_in, dpi = 600,
         device = ragg::agg_png, bg = "white")
  ggsave(out_path("01_study_area", paste0(file_stem, ".pdf")), plot,
         width = width_in, height = height_in, device = cairo_pdf)
}
