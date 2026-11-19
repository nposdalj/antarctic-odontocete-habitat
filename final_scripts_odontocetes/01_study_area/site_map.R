# =============================================================================
# Script  : site_map.R
# Section : 2.1 Study area
# Purpose : Map of the three HARP recording sites (EI, KGI, CI) north of the
#           South Shetland Islands on GMRT bathymetry, with isobaths, land,
#           scale bar, north arrow, and a locator inset showing the study area
#           in its Antarctic Peninsula / Southern Ocean context.
# Inputs  : Environmental Data/GMRT/GMRTbathymetry_highres_widerView.grd
#             (GMRT MapTool grid, legacy GMT/netCDF format; read via
#             map_helpers.R)
#           final_scripts_odontocetes/config.R (SITES, SITE_ORDER, SITE_COLORS)
#           final_scripts_odontocetes/01_study_area/map_helpers.R (shared map
#             settings and layers)
#           rnaturalearth / rnaturalearthdata (coastlines for the inset;
#             bundled with the packages, no download needed)
# Outputs : outputs/01_study_area/site_map.png (600 dpi)
#           outputs/01_study_area/site_map.pdf
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Mean surface currents and ACC fronts are drawn on a separate map
#           (currents_map.R). The +/- 40 km environmental boxes are
#           intentionally not drawn.
#           Bathymetry source citation:
#           [INSERT GMRT citation, e.g., Ryan et al. (2009) GMRT version and
#           download date]
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "01_study_area", "map_helpers.R"))

suppressPackageStartupMessages({
  library(rnaturalearth)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# 1. Main map
# -----------------------------------------------------------------------------
# Bathymetry, land, isobaths, sites, scale bar, and north arrow, built from
# the shared layers in map_helpers.R.
bathy <- load_bathymetry()

main_map <- ggplot() +
  bathymetry_layers(bathy) +
  site_layers() +
  map_frame()

# -----------------------------------------------------------------------------
# 2. Locator inset
# -----------------------------------------------------------------------------
# Lambert azimuthal equal-area projection centred on the study area, so north
# is up and the Antarctic Peninsula, Drake Passage, and southern South America
# keep their familiar shapes. Natural Earth 1:50m land (bundled with
# rnaturalearthdata) is sufficient at this scale.
inset_crs <- "+proj=laea +lat_0=-61 +lon_0=-60 +datum=WGS84 +units=m"

land_inset <- ne_countries(scale = 50, returnclass = "sf") |>
  st_geometry() |>
  st_make_valid() |>
  st_crop(xmin = -90, xmax = -30, ymin = -80, ymax = -45) |>
  st_union() |>
  st_transform(inset_crs)

# Study-area rectangle, densified so its edges follow the projected graticule
study_rect <- st_bbox(c(xmin = MAP_XLIM[1], xmax = MAP_XLIM[2],
                        ymin = MAP_YLIM[1], ymax = MAP_YLIM[2]), crs = 4326) |>
  st_as_sfc() |>
  st_segmentize(units::set_units(10, km)) |>
  st_transform(inset_crs)

# Inset window (projected metres around the centre): Tierra del Fuego to
# ~68 S on the Peninsula
inset_xlim <- c(-1.3e6, 1.1e6)
inset_ylim <- c(-0.85e6, 0.85e6)

# Place-name labels positioned in lon/lat and projected with the map
inset_labels <- tibble(
  label = c("Antarctic\nPeninsula", "South\nAmerica", "Drake\nPassage"),
  lon   = c(-55, -62.5, -64),
  lat   = c(-65.5, -54.6, -58.6)
) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(inset_crs)

inset_map <- ggplot() +
  geom_sf(data = land_inset, fill = LAND_COLOR, color = "grey30",
          linewidth = 0.1) +
  geom_sf(data = study_rect, fill = NA, color = "red3", linewidth = 0.7) +
  geom_sf_text(data = inset_labels, aes(label = label), size = 2.2,
               color = INK_COLOR, lineheight = 0.85, fontface = "italic") +
  coord_sf(crs = inset_crs, xlim = inset_xlim, ylim = inset_ylim,
           expand = FALSE, datum = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "aliceblue", color = NA),
        panel.border = element_rect(fill = NA, color = INK_COLOR,
                                    linewidth = 0.5))

# -----------------------------------------------------------------------------
# 3. Combine and save
# -----------------------------------------------------------------------------
# Inset sits over deep water of the Drake Passage (top left), away from the
# sites and the shelf.
site_map <- main_map +
  inset_element(inset_map, left = 0.005, bottom = 0.69, right = 0.3,
                top = 0.995, align_to = "panel")

save_map(site_map, "site_map", width_in = 9, height_in = 5.6)
