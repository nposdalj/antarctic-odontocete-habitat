# =============================================================================
# Script  : study_area_currents_map.R
# Section : 2.1 Study area
# Purpose : Study-area figure (main text): the three HARP sites (EI, KGI, CI)
#           on GMRT bathymetry with the long-term mean surface currents from
#           the Copernicus GLORYS12V1 reanalysis, the climatological ACC
#           fronts of Orsi et al. (1995), and the schematic circulation of
#           Damini et al. (2025, Fig. 1) digitized by
#           digitize_schematic_pathways.R (mCDW intrusion, Bransfield Current
#           / TBW and its recirculation, Antarctic Coastal Current / TWW,
#           Antarctic Slope Front and its bifurcation eddy, Peninsula Front),
#           plus a locator inset, scale bar, and north arrow.
# Inputs  : DOWNLOAD_DIR/copernicus_regional/phy_surface_currents_P1M_<years>.nc
#             (monthly surface uo/vo, GLOBAL_MULTIYEAR_PHY_001_030; written
#             and validated by currents_map.R, which must be run first)
#           01_study_area/schematic_pathways_digitized.csv (written by
#             digitize_schematic_pathways.R)
#           Environmental Data/GMRT/GMRTbathymetry_highres_widerView.grd
#           R packages orsifronts (0.2.0; Orsi et al. 1995 fronts) and
#             rnaturalearth / rnaturalearthdata (inset coastlines)
#           final_scripts_odontocetes/config.R, 01_study_area/map_helpers.R
# Outputs : outputs/01_study_area/study_area_currents_map_schematic.png
#             (600 dpi) / .pdf  -- the manuscript figure
#           outputs/01_study_area/study_area_currents_map.png / .pdf
#             (same map without the schematic pathways)
# Author  : Natalie Posdaljian
# Date    : 2026-11-18
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Mean currents average all months of the climatology window
#             (earliest deployment year - CLIM_HALF_WINDOW_YEARS to latest
#             deployment year + CLIM_HALF_WINDOW_YEARS, i.e., 2009-2021) at
#             the 0.494 m model level.
#           - The script stops if any schematic pathway crosses land (GMRT
#             elevation >= 0 m) within the map.
#           - Data citations: E.U. Copernicus Marine Service Information,
#             GLOBAL_MULTIYEAR_PHY_001_030 (GLORYS12V1),
#             doi:10.48670/moi-00021; Orsi et al. (1995) Deep-Sea Res. I
#             42:641-673; Damini et al. (2025) Deep-Sea Res. I 222, 104516,
#             doi:10.1016/j.dsr.2025.104516;
#             [INSERT GMRT citation, e.g., Ryan et al. (2009), version and
#             download date].
#           - site_map.R and currents_map.R remain the stand-alone versions.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "01_study_area", "map_helpers.R"))

suppressPackageStartupMessages({
  library(orsifronts)
  library(patchwork)
  library(ggnewscale)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
clim_years <- c(min(year(SITES$rec_start)) - CLIM_HALF_WINDOW_YEARS,
                max(year(SITES$rec_end)) + CLIM_HALF_WINDOW_YEARS)
currents_file <- file.path(DOWNLOAD_DIR, "copernicus_regional",
                           sprintf("phy_surface_currents_P1M_%d-%d.nc",
                                   clim_years[1], clim_years[2]))
if (!file.exists(currents_file)) {
  stop("Regional current file not found: ", currents_file,
       "\nRun 01_study_area/currents_map.R first (downloads and validates it).")
}
pathways_file <- file.path(FINAL_DIR, "01_study_area",
                           "schematic_pathways_digitized.csv")

# Lighter bathymetry ramp (mako from its mid-tones to white) so vectors,
# fronts, sites, and labels stay legible on top of it.
light_depth_palette <- viridisLite::mako(256, begin = 0.35, direction = -1)

# Vectors: one arrow every `arrow_stride` model cells (1/12 deg); arrow
# length = `arrow_deg_per_ms` degrees of latitude per 1 m/s, with the
# longitude component stretched by 1/cos(latitude) so direction and length
# are true on the map. Dark grey and thin, so they stay visible but read
# differently from the coloured schematic pathways.
arrow_stride     <- 4
arrow_deg_per_ms <- 1.5
reference_speed  <- 0.2   # m/s, legend arrow
vector_color     <- "grey20"
vector_width     <- 0.35

# Also write the version with the schematic pathways (manuscript figure)
SHOW_SCHEMATICS <- TRUE

en_dash <- intToUtf8(8211)

# -----------------------------------------------------------------------------
# 2. Long-term mean surface currents
# -----------------------------------------------------------------------------
# ncdf4 returns (lon, lat, time); as.vector runs lon fastest, then lat, then
# time, matching expand_grid(time, lat, lon) (last column fastest).
read_currents <- function(file) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  expand_grid(time = ncvar_get(nc, "time"),
              lat  = ncvar_get(nc, "latitude"),
              lon  = ncvar_get(nc, "longitude")) |>
    mutate(uo = as.vector(ncvar_get(nc, "uo")),
           vo = as.vector(ncvar_get(nc, "vo")))
}

# Mean over all months; cells that are land in any month are dropped.
mean_currents <- read_currents(currents_file) |>
  group_by(lon, lat) |>
  summarise(n_months = sum(!is.na(uo)),
            uo = mean(uo), vo = mean(vo), .groups = "drop") |>
  filter(n_months == max(n_months))

thin <- function(v) sort(unique(v))[c(TRUE, rep(FALSE, arrow_stride - 1))]
arrows_df <- mean_currents |>
  filter(lon %in% thin(lon), lat %in% thin(lat)) |>
  mutate(lon_end = lon + uo * arrow_deg_per_ms / cos(lat * pi / 180),
         lat_end = lat + vo * arrow_deg_per_ms)

# Reference-arrow key (top-left corner of the map), citing the data source
key_box <- c(xmin = -60.97, xmax = -58.95, ymin = -60.47, ymax = -60.03)
ref_len_lon <- reference_speed * arrow_deg_per_ms / cos(-60.38 * pi / 180)
ref_arrow <- tibble(lon = mean(key_box[c("xmin", "xmax")]) - ref_len_lon / 2,
                    lat = -60.38, lon_end = lon + ref_len_lon, lat_end = lat)
ref_label <- paste0("Mean surface current ", clim_years[1], en_dash,
                    clim_years[2], "\n(GLORYS12V1): ", reference_speed,
                    " m s", intToUtf8(c(8315, 185)))

# -----------------------------------------------------------------------------
# 3. ACC fronts (Orsi et al. 1995)
# -----------------------------------------------------------------------------
# Only the fronts that cross the map: Southern ACC Front (SACCF) and Southern
# Boundary of the ACC (SB). The Polar Front lies north of the map. Charcoal,
# distinguished by line type, so they are not confused with the coloured
# schematic pathways.
front_names <- c(saccf = "SACCF (Orsi et al. 1995)",
                 sbdy  = "SB (Orsi et al. 1995)")
front_colors <- setNames(c("grey20", "grey20"), front_names)

fronts <- st_as_sf(orsifronts) |>
  st_make_valid() |>
  filter(front %in% names(front_names)) |>
  st_set_agr("constant") |>
  st_crop(xmin = MAP_XLIM[1] - 1, ymin = MAP_YLIM[1] - 1,
          xmax = MAP_XLIM[2] + 1, ymax = MAP_YLIM[2] + 1) |>
  mutate(front_label = factor(front_names[as.character(front)],
                              levels = front_names))

# -----------------------------------------------------------------------------
# 4. Schematic pathways (digitized from Damini et al. 2025, Fig. 1)
# -----------------------------------------------------------------------------
# Style per digitized feature. `legend` text names the pathway and the
# source of its digitized position; features sharing a legend entry share
# its colour and line type. The Bransfield Current recirculation north of
# the islands is drawn as in Damini et al. (2025); Moffat & Meredith (2018)
# describe a southwestward "Polar Slope Current" (200-600 m) on the same
# northern slope.
pathway_styles <- tribble(
  ~feature,             ~legend,                                  ~color,    ~linetype,
  "mCDW intrusion",     "mCDW intrusion (Damini et al. 2025)",    "#C62828", "dashed",
  "Bransfield Current", "BC / TBW (Damini et al. 2025)",          "#6A1B9A", "solid",
  "BC recirculation",   "BC recirculation (Damini et al. 2025)",  "#6A1B9A", "22",
  "TWW / ACoC",         "ACoC / TWW (Damini et al. 2025)",        "#E65100", "solid",
  "ASF",                "ASF and eddy (Damini et al. 2025)",      "#5D4037", "solid",
  "ASF eddy",           "ASF and eddy (Damini et al. 2025)",      "#5D4037", "solid",
  "Peninsula Front",    "Peninsula Front (Damini et al. 2025)",   "#00796B", "42"
)
legend_levels <- unique(pathway_styles$legend)
legend_colors <- setNames(pathway_styles$color, pathway_styles$legend)[
  legend_levels]
legend_types  <- setNames(pathway_styles$linetype, pathway_styles$legend)[
  legend_levels]

pathways <- read_csv(pathways_file, show_col_types = FALSE) |>
  inner_join(pathway_styles, by = "feature") |>
  arrange(feature, segment, order)

# Clip every pathway to the map (inset slightly so arrowheads are not cut by
# the frame). Parts outside the map are dropped; a pathway that leaves and
# re-enters the map becomes several pieces. The arrowhead is drawn at the
# end of the last piece of each pathway (flow direction = point order).
clip_margin <- 0.06
inside <- function(lon, lat) {
  lon > MAP_XLIM[1] + clip_margin & lon < MAP_XLIM[2] - clip_margin &
    lat > MAP_YLIM[1] + clip_margin & lat < MAP_YLIM[2] - clip_margin
}

pathway_pieces <- pathways |>
  group_by(feature, segment) |>
  mutate(is_in = inside(lon, lat),
         piece = cumsum(is_in != lag(is_in, default = first(is_in)))) |>
  filter(is_in) |>
  group_by(feature, segment, piece) |>
  filter(n() >= 2) |>
  ungroup() |>
  group_by(feature, segment) |>
  mutate(has_arrow = piece == max(piece)) |>
  ungroup() |>
  mutate(path_id = paste(feature, segment, piece),
         legend = factor(legend, levels = legend_levels))

# Land check of what is drawn: no schematic pathway may cross land.
bathy <- load_bathymetry()
drawn_on_land <- pathway_pieces |>
  group_by(path_id) |>
  group_modify(function(d, key) {
    n <- ceiling(sum(sqrt(diff(d$lon)^2 + diff(d$lat)^2)) / 0.005) + 1
    xy <- cbind(approx(seq_along(d$lon), d$lon, n = n)$y,
                approx(seq_along(d$lat), d$lat, n = n)$y)
    tibble(n_on_land = sum(terra::extract(bathy, xy)[, 1] >= 0,
                           na.rm = TRUE))
  }) |>
  ungroup()
if (any(drawn_on_land$n_on_land > 0)) {
  print(filter(drawn_on_land, n_on_land > 0))
  stop("Schematic pathway(s) cross land; see digitize_schematic_pathways.R")
}

# Short on-map labels (water-mass / current acronyms), positioned next to
# their pathway and away from sites and other lines.
pathway_labels <- tribble(
  ~label,      ~lon,    ~lat,    ~angle, ~legend,
  "mCDW",      -56.85,  -61.45,  0,      "mCDW intrusion (Damini et al. 2025)",
  "BC (TBW)",  -59.10,  -62.45,  33,     "BC / TBW (Damini et al. 2025)",
  "ACoC (TWW)", -59.0,  -63.33,  17,     "ACoC / TWW (Damini et al. 2025)",
  "ASF",       -52.25,  -62.65,  0,      "ASF and eddy (Damini et al. 2025)",
  "PenF",      -57.75,  -62.75,  27,     "Peninsula Front (Damini et al. 2025)"
)

# -----------------------------------------------------------------------------
# 5. Map
# -----------------------------------------------------------------------------
# Layer order (bottom to top): bathymetry, mean current vectors, ACC fronts,
# schematic pathways (optional), sites, labels, reference key.
vector_head    <- arrow(length = unit(0.07, "cm"), type = "closed", angle = 25)
schematic_head <- arrow(length = unit(0.2, "cm"), type = "closed", angle = 25)

# Site labels moved off the schematic pathways where needed
schematic_site_labels <- list(KGI = c(-0.52, 0.15), EI = c(0.30, 0.10),
                              CI = c(0.32, 0.13))

build_map <- function(show_schematics) {
  schematic_layers <- if (show_schematics) {
    list(
      # ACC fronts share the schematic colour/linetype scales so that fronts
      # and pathways appear in ONE untitled legend (fronts listed first).
      geom_sf(data = fronts, aes(color = front_label, linetype = front_label),
              linewidth = 0.8),
      # One layer per piece so each keeps its own arrowhead
      lapply(split(pathway_pieces, pathway_pieces$path_id), function(p) {
        geom_path(data = p, aes(x = lon, y = lat, color = legend,
                                linetype = legend),
                  linewidth = 0.75, lineend = "round", linejoin = "round",
                  arrow = if (p$has_arrow[1] &&
                              p$feature[1] != "Peninsula Front") {
                    schematic_head
                  })
      }),
      scale_color_manual(
        values = c(SITE_COLORS, front_colors, legend_colors),
        breaks = unname(c(front_names, legend_levels)), name = NULL, drop = FALSE,
        guide = guide_legend(order = 2, override.aes = list(linewidth = 0.9))
      ),
      scale_linetype_manual(
        values = c(legend_types, setNames(c("longdash", "dotdash"),
                                          front_names)),
        breaks = unname(c(front_names, legend_levels)), name = NULL,
        drop = FALSE, guide = guide_legend(order = 2)
      ),
      new_scale("linetype")
    )
  }

  label_layers <- if (show_schematics) {
    lapply(seq_len(nrow(pathway_labels)), function(i) {
      l <- pathway_labels[i, ]
      geom_label(data = l, aes(x = lon, y = lat, label = label),
                 angle = l$angle, color = legend_colors[[l$legend]],
                 fill = alpha("white", 0.8), linewidth = 0, size = 2.9,
                 fontface = "bold", label.padding = unit(0.08, "lines"))
    })
  }

  ggplot() +
    bathymetry_layers(bathy, palette = light_depth_palette) +
    # Mean surface current vectors (data)
    geom_segment(data = arrows_df,
                 aes(x = lon, y = lat, xend = lon_end, yend = lat_end),
                 color = vector_color, linewidth = vector_width,
                 arrow = vector_head) +
    schematic_layers +
    # ACC fronts (Orsi et al. 1995): drawn here only on the plain map; on the
    # schematic map they are part of schematic_layers (single legend).
    (if (!show_schematics) list(
      geom_sf(data = fronts, aes(linetype = front_label),
              color = front_colors[as.character(fronts$front_label)],
              linewidth = 0.8),
      scale_linetype_manual(
        values = setNames(c("longdash", "dotdash"), front_names),
        name = "ACC fronts",
        guide = guide_legend(order = 2, override.aes = list(
          color = unname(front_colors), linewidth = 0.8)))
    )) +
    # Sites (the site colour scale comes from site_layers() unless the
    # schematic colour scale above already maps both)
    (if (show_schematics) site_layers(schematic_site_labels)[-2] else
       site_layers()) +
    label_layers +
    # Reference-arrow key
    annotate("rect", xmin = key_box[["xmin"]], xmax = key_box[["xmax"]],
             ymin = key_box[["ymin"]], ymax = key_box[["ymax"]],
             fill = "white", color = "grey40", linewidth = 0.2) +
    annotate("text", x = mean(key_box[c("xmin", "xmax")]), y = -60.19,
             label = ref_label, size = 2.5, color = INK_COLOR,
             lineheight = 0.9) +
    geom_segment(data = ref_arrow,
                 aes(x = lon, y = lat, xend = lon_end, yend = lat_end),
                 color = vector_color, linewidth = vector_width + 0.1,
                 arrow = vector_head) +
    # Scale bar and north arrow move to the top right on the schematic map,
    # where they do not overlap the Antarctic Slope Front.
    map_frame(corner = if (show_schematics) "tr" else "br") +
    # Compact legends (top right) so the locator inset fits below them
    guides(fill = guide_colourbar(order = 1, theme = theme(
      legend.key.height = unit(2, "cm")))) +
    theme(legend.position = "right",
          legend.justification = "top",
          legend.key.height = unit(0.42, "cm"),
          legend.key.width = unit(1, "cm"),
          legend.text = element_text(size = 8.5),
          legend.title = element_text(size = 10))
}

# -----------------------------------------------------------------------------
# 6. Locator inset and save
# -----------------------------------------------------------------------------
# The inset sits below the legends, outside the map panel, so it hides no
# current vectors.
add_inset <- function(map) {
  map + inset_element(locator_inset(), left = 0.79, bottom = 0.02,
                      right = 0.995, top = 0.30, align_to = "full")
}

save_map(add_inset(build_map(show_schematics = FALSE)),
         "study_area_currents_map", width_in = 10.5, height_in = 5.8)

if (SHOW_SCHEMATICS) {
  save_map(add_inset(build_map(show_schematics = TRUE)),
           "study_area_currents_map_schematic", width_in = 11, height_in = 6)
}
