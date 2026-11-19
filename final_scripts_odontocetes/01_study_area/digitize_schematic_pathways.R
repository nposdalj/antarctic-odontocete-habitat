# =============================================================================
# Script  : digitize_schematic_pathways.R
# Section : 2.1 Study area
# Purpose : Digitize the schematic circulation pathways and fronts of
#           Damini et al. (2025, Fig. 1) for the study-area map: download the
#           source figure, georeference it from its longitude/latitude ticks,
#           trace each pathway by colour segmentation along a rough guide
#           line, convert to lon/lat, check that no pathway crosses land, and
#           save the coordinates as a CSV read by study_area_currents_map.R.
# Inputs  : Damini et al. (2025) accepted manuscript (NERC Open Research
#             Archive, https://nora.nerc.ac.uk/id/eprint/539533), page 3
#             (Fig. 1); downloaded on first run to
#             DOWNLOAD_DIR/schematic_figures/ and its embedded figure image
#             extracted with `pdfimages` (poppler; e.g., bundled with MiKTeX)
#           Environmental Data/GMRT/GMRTbathymetry_highres_widerView.grd
#             (land check, via map_helpers.R)
#           final_scripts_odontocetes/config.R, 01_study_area/map_helpers.R
# Outputs : final_scripts_odontocetes/01_study_area/schematic_pathways_digitized.csv
#             feature, segment, order, lon, lat, source_figure, citation
#           outputs/01_study_area/schematic_georeference_fit.csv
#             tick positions, fitted positions, residuals (pixels, degrees)
#           outputs/01_study_area/schematic_land_check.csv
#           outputs/01_study_area/schematic_digitizing_check.png
#             source figure with the traced points and the GMRT coastline
#             (projected into the figure) overlaid, for visual audit
# Author  : Natalie Posdaljian
# Date    : 2026-11-18
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Source figure: Damini B.Y., Brum A.L., Hall R.A., Dotto T.S.,
#             Azevedo J.L.L., Heywood K.J., Mata M.M., Garcia C.A.E., Kerr R.
#             (2025) Deep-Sea Res. I 222, 104516,
#             doi:10.1016/j.dsr.2025.104516, Fig. 1 (modified from Damini et
#             al. 2022, doi:10.1016/j.dsr.2021.103667). The figure itself is
#             not redistributed; only the digitized coordinates are.
#           - The map in Fig. 1 is a Mercator projection (latitude tick
#             spacing grows with latitude; a Mercator fit leaves residuals
#             < 0.4 px, a linear fit up to 5 px).
#           - Guide lines below are rough pixel positions read from the
#             figure; the traced coordinates come from the coloured pixels
#             within `radius` pixels of the guide, so guide placement only
#             needs to be within a few pixels of the drawn line.
#           - Arrowheads are excluded from the traces (guides stop short of
#             them); flow direction = order of the points.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "01_study_area", "map_helpers.R"))

suppressPackageStartupMessages({
  library(jpeg)
  library(magick)
})

# -----------------------------------------------------------------------------
# 1. Source figure
# -----------------------------------------------------------------------------
source_citation <- "Damini et al. (2025)"
source_figure   <- "Damini et al. (2025) Fig. 1"
pdf_url  <- paste0("https://nora.nerc.ac.uk/id/eprint/539533/1/",
                   "Damini_et_al_2025_revised.pdf")
fig_dir  <- file.path(DOWNLOAD_DIR, "schematic_figures")
pdf_file <- file.path(fig_dir, "Damini_et_al_2025_revised.pdf")
fig_file <- file.path(fig_dir, "damini2025_fig1.jpg")

# Download the accepted manuscript and extract the embedded Fig. 1 image
# (page 3) once; later runs reuse the extracted file.
if (!file.exists(fig_file)) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(pdf_file)) download.file(pdf_url, pdf_file, mode = "wb")
  prefix <- file.path(fig_dir, "damini2025_p3")
  status <- system2("pdfimages", c("-f", "3", "-l", "3", "-j",
                                   shQuote(pdf_file), shQuote(prefix)))
  if (status != 0) stop("pdfimages failed; install poppler (or MiKTeX)")
  file.rename(paste0(prefix, "-000.jpg"), fig_file)
}

fig <- readJPEG(fig_file)
if (!identical(dim(fig), c(815L, 1285L, 3L))) {
  stop("Unexpected figure size ", paste(dim(fig), collapse = " x "),
       ": the guides below were read from the 1285 x 815 px image")
}
n_row <- dim(fig)[1]
n_col <- dim(fig)[2]

# -----------------------------------------------------------------------------
# 2. Georeference from the axis ticks
# -----------------------------------------------------------------------------
# Tick marks are short dark segments just outside the map frame. Their
# centres are located automatically in a strip outside the bottom frame
# (longitude ticks) and the left frame (latitude ticks), then matched to the
# tick labels printed in the figure.
dark <- fig[, , 1] < 0.35 & fig[, , 2] < 0.35 & fig[, , 3] < 0.35

run_centres <- function(v) {
  r <- rle(v)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  ((starts + ends) / 2)[r$values]
}

lon_tick_px <- run_centres(colSums(dark[627:632, ]) >= 5)
lon_tick_px <- lon_tick_px[lon_tick_px > 90]          # drop text left of map
lat_tick_px <- run_centres(rowSums(dark[, 87:92]) >= 5)

lon_ticks <- tibble(axis = "lon", value = c(-66, -62, -58, -54, -50),
                    px = lon_tick_px)
lat_ticks <- tibble(axis = "lat", value = c(-61, -62, -63, -64, -65),
                    px = lat_tick_px)
stopifnot(nrow(lon_ticks) == 5, nrow(lat_ticks) == 5)

# Longitude: linear in x. Latitude: linear in the Mercator ordinate.
mercator <- function(lat) log(tan(pi / 4 + lat * pi / 360))
inv_mercator <- function(m) (2 * atan(exp(m)) - pi / 2) * 180 / pi

fit_lon <- lm(px ~ value, data = lon_ticks)
fit_lat <- lm(px ~ mercator(value), data = lat_ticks)

px_to_lon <- function(x) (x - coef(fit_lon)[1]) / coef(fit_lon)[2]
px_to_lat <- function(y) inv_mercator((y - coef(fit_lat)[1]) / coef(fit_lat)[2])
lon_to_px <- function(lon) coef(fit_lon)[1] + coef(fit_lon)[2] * lon
lat_to_px <- function(lat) coef(fit_lat)[1] + coef(fit_lat)[2] * mercator(lat)

georef_fit <- bind_rows(
  lon_ticks |> mutate(fitted_px = fitted(fit_lon),
                      residual_deg = px_to_lon(px) - value),
  lat_ticks |> mutate(fitted_px = fitted(fit_lat),
                      residual_deg = px_to_lat(px) - value)
) |>
  mutate(residual_px = px - fitted_px)
write_csv(georef_fit, out_path("01_study_area",
                               "schematic_georeference_fit.csv"))
message(sprintf("Georeference: max |residual| %.2f px, %.4f deg",
                max(abs(georef_fit$residual_px)),
                max(abs(georef_fit$residual_deg))))

# -----------------------------------------------------------------------------
# 3. Colour classes
# -----------------------------------------------------------------------------
# Pathway colours in Fig. 1: red dashed = CDW, orange = Bransfield Current
# and its recirculation, purple = Dense Shelf Water / Antarctic Coastal
# Current, blue = ASF bifurcation eddy, black solid = Antarctic Slope Front,
# black dashed = Peninsula Front. Classes are hue ranges on saturated pixels.
hsv <- rgb2hsv(rbind(c(fig[, , 1]), c(fig[, , 2]), c(fig[, , 3])),
               maxColorValue = 1)
hue <- matrix(hsv[1, ], n_row)
sat <- matrix(hsv[2, ], n_row)
val <- matrix(hsv[3, ], n_row)

hue_class <- function(h_min, h_max, s_min = 0.55, v_min = 0.35) {
  in_hue <- if (h_min < h_max) hue >= h_min & hue <= h_max else
    hue >= h_min | hue <= h_max
  sat > s_min & val > v_min & in_hue
}

colour_masks <- list(
  red    = hue_class(0.96, 0.03),
  orange = hue_class(0.04, 0.11),
  purple = hue_class(0.70, 0.80, s_min = 0.40),
  blue   = hue_class(0.57, 0.68, s_min = 0.50, v_min = 0.50) & !dark,
  black  = fig[, , 1] < 0.3 & fig[, , 2] < 0.3 & fig[, , 3] < 0.3
)

# -----------------------------------------------------------------------------
# 4. Guide lines (pixel coordinates, flow order)
# -----------------------------------------------------------------------------
# feature names match the feature table in study_area_currents_map.R.
# Mercator figure pixels: x right, y down.
guides <- tribble(
  ~feature,           ~segment, ~colour,  ~radius, ~guide,
  "mCDW intrusion",   1, "red",    5, "613,97; 641,94; 661,102; 676,123; 685,147; 686,160",
  "Bransfield Current", 1, "orange", 9, "445,336; 479,319; 510,304; 551,275; 599,241; 649,225; 678,210; 715,183; 752,157; 781,136; 805,117",
  "BC recirculation", 1, "orange", 9, "676,170; 662,154; 640,147; 610,154; 563,175; 518,196; 495,215; 479,228; 452,230; 425,228; 402,240; 383,264; 378,296; 392,318; 412,327",
  "TWW / ACoC",       1, "purple", 8, "799,617; 802,570; 810,545; 830,512; 850,490; 890,452; 925,427; 950,407; 965,390; 967,367; 962,350; 947,320; 925,295; 904,272; 872,252; 856,245",
  "TWW / ACoC",       2, "purple", 8, "712,322; 692,318; 660,330; 640,344; 600,358; 570,366; 510,386; 472,401",
  "TWW / ACoC",       3, "purple", 8, "800,555; 792,530; 792,500; 780,485; 752,465; 745,445; 742,430",
  "ASF",              1, "black",  7, "1050,617; 1012,575; 962,550; 945,525; 941,482; 950,460; 975,445; 1000,437; 1050,412; 1084,401; 1148,390; 1168,366; 1152,350; 1084,340; 1066,338; 1050,302; 1020,270; 1000,258; 970,252; 944,240; 926,224; 922,204; 930,190; 950,170; 960,162",
  "Peninsula Front",  1, "black",  5, "452,376; 500,362; 544,342; 592,316; 640,286; 680,266; 704,258; 760,258; 786,252; 800,234; 812,204; 832,182; 860,162; 866,142; 870,110; 866,68"
)

parse_px <- function(s) {
  xy <- str_split_fixed(str_split(s, "\\s*;\\s*")[[1]], "\\s*,\\s*", 2)
  tibble(x = as.numeric(xy[, 1]), y = as.numeric(xy[, 2]))
}

# -----------------------------------------------------------------------------
# 5. Trace each pathway along its guide
# -----------------------------------------------------------------------------
# The guide is densified to 1-px steps; at each step the centroid of the
# class pixels within `radius` px is taken (steps with fewer than 3 class
# pixels, e.g., gaps of dashed lines, are skipped). The centroids are then
# smoothed with a running mean (window 7) and thinned to every 4th point.
densify <- function(g, step = 1) {
  seg_len <- sqrt(diff(g$x)^2 + diff(g$y)^2)
  s_in  <- c(0, cumsum(seg_len))
  s_out <- seq(0, max(s_in), by = step)
  tibble(s = s_out,
         x = approx(s_in, g$x, s_out)$y,
         y = approx(s_in, g$y, s_out)$y)
}

trace_guide <- function(guide, mask, radius) {
  pts <- densify(parse_px(guide))
  hits <- map_dfr(seq_len(nrow(pts)), function(i) {
    x0 <- round(pts$x[i])
    y0 <- round(pts$y[i])
    xs <- max(1, x0 - radius):min(n_col, x0 + radius)
    ys <- max(1, y0 - radius):min(n_row, y0 + radius)
    win <- which(mask[ys, xs, drop = FALSE], arr.ind = TRUE)
    if (nrow(win) == 0) return(NULL)
    px <- xs[win[, 2]]
    py <- ys[win[, 1]]
    keep <- (px - pts$x[i])^2 + (py - pts$y[i])^2 <= radius^2
    if (sum(keep) < 3) return(NULL)
    tibble(s = pts$s[i], x = mean(px[keep]), y = mean(py[keep]))
  })
  smooth <- function(v) {
    out <- stats::filter(v, rep(1 / 7, 7), sides = 2)
    ifelse(is.na(out), v, out)
  }
  hits |>
    mutate(x = smooth(x), y = smooth(y)) |>
    filter(row_number() %% 4 == 1 | row_number() == n())
}

traces <- guides |>
  mutate(trace = pmap(list(guide, colour, radius), function(g, col, r) {
    trace_guide(g, colour_masks[[col]], r)
  })) |>
  select(feature, segment, trace) |>
  unnest(trace)

# ASF bifurcation eddy (blue): the two arcs are traced by angle around the
# centre of the blue pixels in the eddy window (mean radius per 6-degree
# bin). Arrowheads show anticyclonic (counterclockwise) rotation in Fig. 1:
# westward on the north side, eastward on the south side.
eddy_window <- list(x = 880:960, y = 90:145)
eddy_px <- which(colour_masks$blue[eddy_window$y, eddy_window$x],
                 arr.ind = TRUE)
eddy_x <- eddy_window$x[eddy_px[, 2]]
eddy_y <- eddy_window$y[eddy_px[, 1]]
eddy_cx <- (min(eddy_x) + max(eddy_x)) / 2
eddy_cy <- (min(eddy_y) + max(eddy_y)) / 2

eddy_trace <- tibble(x = eddy_x, y = eddy_y) |>
  # Angle in geographic sense (y axis flipped so north is up)
  mutate(angle = atan2(eddy_cy - y, x - eddy_cx) * 180 / pi,
         bin = floor(angle / 6) * 6 + 3) |>
  group_by(bin) |>
  summarise(x = mean(x), y = mean(y), n = n(), .groups = "drop") |>
  filter(n >= 2) |>
  arrange(bin)

# Split into the two drawn arcs at the largest angular gaps and order each
# arc counterclockwise (increasing geographic angle).
gaps <- diff(c(eddy_trace$bin, eddy_trace$bin[1] + 360))
breaks <- sort(order(gaps, decreasing = TRUE)[1:2])
arc_id <- rep(2, nrow(eddy_trace))
arc_id[(breaks[1] + 1):breaks[2]] <- 1
eddy_traces <- eddy_trace |>
  mutate(segment = arc_id, feature = "ASF eddy",
         # unwrap the arc that crosses +/-180 degrees
         s = if_else(segment == 2 & bin <= eddy_trace$bin[breaks[1]],
                     bin + 360, bin)) |>
  arrange(segment, s) |>
  select(feature, segment, s, x, y)

# -----------------------------------------------------------------------------
# 6. Convert to lon/lat
# -----------------------------------------------------------------------------
pathways_raw <- bind_rows(traces, eddy_traces) |>
  mutate(lon = px_to_lon(x), lat = px_to_lat(y))

# -----------------------------------------------------------------------------
# 7. Land check against GMRT, end trimming, and save
# -----------------------------------------------------------------------------
# Every segment between consecutive traced points is sampled at ~0.005 deg
# and checked against GMRT (elevation >= 0 m = land). Where a drawn line in
# the source figure overlaps an island outline at the START or END of a
# pathway (a line-width/generalisation conflict, e.g., the Bransfield
# Current passing Deception Island in Fig. 1), the land-touching end
# segments are trimmed. A land crossing in the INTERIOR of a pathway stops
# the script: it would mean a tracing or georeferencing error that must be
# fixed in the guides, not hidden.
bathy_check <- load_bathymetry(xlim = c(-62, -50), ylim = c(-65.2, -59.5))
check_step_deg <- 0.005

segment_on_land <- function(lon, lat) {
  map_lgl(seq_len(length(lon) - 1), function(i) {
    len <- sqrt((lon[i + 1] - lon[i])^2 + (lat[i + 1] - lat[i])^2)
    n <- max(2, ceiling(len / check_step_deg) + 1)
    xy <- cbind(seq(lon[i], lon[i + 1], length.out = n),
                seq(lat[i], lat[i + 1], length.out = n))
    any(terra::extract(bathy_check, xy)[, 1] >= 0, na.rm = TRUE)
  })
}

# "End" = the first or last `end_fraction` of the pathway's points.
end_fraction <- 0.1

trim_land_ends <- function(d, key) {
  on_land <- segment_on_land(d$lon, d$lat)
  n_seg <- length(on_land)
  land <- which(on_land)
  near_start <- land[land <= ceiling(end_fraction * n_seg)]
  near_end   <- land[land > floor((1 - end_fraction) * n_seg)]
  if (length(setdiff(land, c(near_start, near_end))) > 0) {
    stop(key$feature, " segment ", key$segment,
         ": interior crossing of land; revise its guide line")
  }
  first_seg <- if (length(near_start)) max(near_start) + 1 else 1
  last_seg  <- if (length(near_end)) min(near_end) - 1 else n_seg
  d[first_seg:(last_seg + 1), ] |>
    mutate(trimmed_start_pts = first_seg - 1,
           trimmed_end_pts   = n_seg - last_seg)
}

pathways <- pathways_raw |>
  group_by(feature, segment) |>
  group_modify(trim_land_ends) |>
  mutate(order = row_number()) |>
  ungroup()

land_check <- pathways |>
  group_by(feature, segment) |>
  summarise(n_points = n(),
            trimmed_start_pts = first(trimmed_start_pts),
            trimmed_end_pts   = first(trimmed_end_pts),
            segments_on_land_after_trim =
              sum(segment_on_land(lon, lat)),
            .groups = "drop")
write_csv(land_check, out_path("01_study_area", "schematic_land_check.csv"))
print(land_check)
stopifnot(all(land_check$segments_on_land_after_trim == 0))

write_csv(pathways |>
            transmute(feature, segment, order,
                      lon = round(lon, 4), lat = round(lat, 4),
                      source_figure = source_figure,
                      citation = source_citation),
          file.path(FINAL_DIR, "01_study_area",
                    "schematic_pathways_digitized.csv"))

# -----------------------------------------------------------------------------
# 8. Audit image
# -----------------------------------------------------------------------------
# Traced points (cyan) and the GMRT 0-m contour projected into the figure's
# pixel space (magenta) over the source figure: the coastline overlay checks
# the georeferencing, the points check the tracing.
coast <- as.contour(bathy_check, levels = 0) |>
  st_as_sf() |>
  st_geometry() |>
  st_cast("LINESTRING") |>
  st_coordinates() |>
  as_tibble()

audit <- image_draw(image_read(fig_file))
for (l in unique(coast$L1)) {
  cl <- filter(coast, L1 == l)
  lines(lon_to_px(cl$X), lat_to_px(cl$Y), col = "magenta", lwd = 0.6)
}
points(pathways$x, pathways$y, col = "cyan", pch = 16, cex = 0.35)
dev.off()
image_write(audit, out_path("01_study_area", "schematic_digitizing_check.png"))
