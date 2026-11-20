# =============================================================================
# Script  : 04b_eke_per_cell.R
# Section : 2.5 Environmental data - eddy kinetic energy (EKE) by depth layer
# Purpose : Compute daily EKE in the surface (0-100 m), mid-water (100-500 m)
#           and deep (500-1000 m) layers over the +/-20 km box at each site.
#           EKE (used in the models; Ducet et al. 2000): at each grid cell
#             and model level, u' and v' are the daily velocities minus that
#             cell's 10-year mean; EKE = 0.5 (u'^2 + v'^2) per cell, averaged
#             over the box cells with water, then over the layer with
#             thickness weights.
#           EKE_box (previous definition, kept for comparison only): EKE of
#             the box-mean velocity, 0.5 (u'^2 + v'^2) with u' the daily
#             box-mean u minus its 10-year mean. Eddies inside the box partly
#             cancel, so this measures variability of the net flow.
# Inputs  : DOWNLOAD_DIR/copernicus/<site>_phy_{uo,vo}_<year>.nc
#           (GLOBAL_MULTIYEAR_PHY_001_030 daily velocities, download region,
#           deployment years +/- 5 years; from 01_download_copernicus.R)
# Outputs : outputs/05_environmental_data/eke_layers.csv
#             site, date, layer, EKE (per-cell), EKE_box  [cm2 s-2]
#             (read by 05_compile_dataset.R and 08_eke_per_cell_comparison.R)
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Run after 04_extract_copernicus.R and before 05_compile_dataset.R.
#           The 10-year mean is taken over the whole downloaded record at
#           each site, as for the previous EKE.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages({ library(tidyverse); library(ncdf4) })

BOX_TOL_DEG <- 1e-6   # as in 04_extract_copernicus.R
boxes <- site_boxes(BOX_HALF_KM)

# Thickness weights of each model level in each layer (as in 05_compile_dataset.R)
level_weights <- function(levels) {
  edges <- c(0, (head(levels, -1) + levels[-1]) / 2, Inf)
  expand_grid(i = seq_along(levels), l = seq_len(nrow(LAYERS))) |>
    mutate(depth_m = levels[i], layer = LAYERS$layer[l],
           weight = pmax(0, pmin(edges[i + 1], LAYERS$bottom_m[l]) -
                            pmax(edges[i], LAYERS$top_m[l]))) |>
    filter(weight > 0) |>
    select(depth_m, layer, weight)
}

# Read one velocity file, cropped to the main box: array cells x depth x time
read_box <- function(file, var, bx) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  lon <- as.vector(ncvar_get(nc, "longitude"))
  lat <- as.vector(ncvar_get(nc, "latitude"))
  i_lon <- which(lon >= bx$lon_min - BOX_TOL_DEG & lon <= bx$lon_max + BOX_TOL_DEG)
  i_lat <- which(lat >= bx$lat_min - BOX_TOL_DEG & lat <= bx$lat_max + BOX_TOL_DEG)
  x <- ncvar_get(nc, var, collapse_degen = FALSE)[i_lon, i_lat, , , drop = FALSE]
  units <- ncatt_get(nc, "time", "units")$value
  step <- c(seconds = 1, hours = 3600, days = 86400)[[word(units, 1)]]
  origin <- as.POSIXct(str_remove(units, "^\\w+ since "), tz = "UTC")
  d <- dim(x)
  list(x = array(x, dim = c(d[1] * d[2], d[3], d[4])),
       depth = round(as.vector(ncvar_get(nc, "depth")), 2),
       dates = as.Date(origin + ncvar_get(nc, "time") * step))
}

# Join yearly arrays (cells x depth x time) along time
bind_time <- function(arrs) {
  d <- dim(arrs[[1]])
  array(unlist(arrs), dim = c(d[1], d[2], sum(map_int(arrs, \(a) dim(a)[3]))))
}

eke_levels <- map_dfr(SITE_ORDER, function(s) {
  bx <- filter(boxes, site == s)
  files_u <- sort(list.files(file.path(DOWNLOAD_DIR, "copernicus"),
                             pattern = paste0("^", s, "_phy_uo_[0-9]{4}[.]nc$"),
                             full.names = TRUE))
  u <- map(files_u, read_box, var = "uo", bx = bx)
  v <- map(sub("_uo_", "_vo_", files_u), read_box, var = "vo", bx = bx)
  U <- bind_time(map(u, "x"))
  V <- bind_time(map(v, "x"))
  dates <- do.call(c, map(u, "dates"))
  depth <- u[[1]]$depth

  # Per-cell EKE: anomalies from each cell's 10-year mean at each level
  # (the cells x depth mean recycles along the time dimension)
  u_mean <- apply(U, c(1, 2), mean, na.rm = TRUE)
  v_mean <- apply(V, c(1, 2), mean, na.rm = TRUE)
  eke_cell <- 0.5 * ((U - as.vector(u_mean))^2 + (V - as.vector(v_mean))^2) * 1e4
  eke_cell <- apply(eke_cell, c(2, 3), mean, na.rm = TRUE)         # depth x time

  # Box-mean EKE (previous definition, for comparison)
  ub <- apply(U, c(2, 3), mean, na.rm = TRUE)
  vb <- apply(V, c(2, 3), mean, na.rm = TRUE)
  eke_box <- 0.5 * ((ub - rowMeans(ub, na.rm = TRUE))^2 +
                    (vb - rowMeans(vb, na.rm = TRUE))^2) * 1e4

  message(s, ": ", length(files_u), " years, ", dim(U)[1], " cells, ",
          length(dates), " days")
  expand_grid(t = seq_along(dates), k = seq_along(depth)) |>
    mutate(site = s, date = dates[t], depth_m = depth[k],
           EKE = eke_cell[cbind(k, t)], EKE_box = eke_box[cbind(k, t)]) |>
    select(site, date, depth_m, EKE, EKE_box)
})

eke_layers <- eke_levels |>
  filter(is.finite(EKE)) |>
  inner_join(level_weights(sort(unique(eke_levels$depth_m))), by = "depth_m",
             relationship = "many-to-many") |>
  group_by(site, date, layer) |>
  summarise(EKE = sum(EKE * weight) / sum(weight),
            EKE_box = sum(EKE_box * weight) / sum(weight), .groups = "drop")
write_csv(eke_layers, out_path("05_environmental_data", "eke_layers.csv"))
message("Wrote eke_layers.csv: ", nrow(eke_layers), " rows")
