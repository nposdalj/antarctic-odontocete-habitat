# =============================================================================
# Script  : amsr2_site_series.R
# Section : 2.3 Sea ice - daily sea-ice metrics at each site (AMSR2)
# Purpose : From the daily AMSR2 ASI sea-ice concentration fields (University
#           of Bremen, 6.25 km), compute for each site and day:
#             ice_conc          mean concentration over the +/-20 km box
#                               (fraction 0-1; the sea-ice variable used in
#                               all analyses)
#             open_water_frac   fraction of box cells with < 15% ice
#             dist_ice_edge_km  signed distance from the recorder to the
#                               nearest ice edge (15% threshold): positive =
#                               recorder in open water, distance to the
#                               nearest cell with >= 15% ice; negative =
#                               recorder within the ice, distance to the
#                               nearest cell with < 15% ice; NA = no ice in
#                               the mapped region (edge > ~400 km away)
#             n_cells           AMSR2 cells in the box (land-flagged excluded)
# Inputs  : DOWNLOAD_DIR/bremen_asi/asi_s6250_<year>_region.tif
#           (bremen_01_download.R)
# Outputs : outputs/03_sea_ice/amsr2_daily_site.csv (read by
#           05_environmental_data/05_compile_dataset.R)
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Decision N. Posdaljian 2026-11-20: sea ice from AMSR2 replaces
#           the Copernicus reanalysis ice everywhere (the reanalysis
#           overestimated concentration near the ice edge relative to AMSR2
#           and MODIS imagery). AMSR2 values > 100 (land/missing flags) -> NA.
#           Box cells = AMSR2 cell centres inside the +/-20 km lon/lat box.
#           Distances are geodesic between cell centres; the recorder's cell
#           is the cell nearest to the recorder position.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages({ library(tidyverse); library(terra) })

EDGE_PCT <- 15
asi <- rast(sort(list.files(file.path(DOWNLOAD_DIR, "bremen_asi"), pattern = "_region[.]tif$",
                            full.names = TRUE)))
dates <- as.Date(names(asi), format = "%Y%m%d")
vals <- values(asi)                       # cells x days
vals[vals > 100] <- NA
cells_ll <- project(as.points(asi[[1]], na.rm = FALSE), "EPSG:4326")
ll <- crds(cells_ll)
boxes <- site_boxes(BOX_HALF_KM)

# Cells used for the ice-edge distance (not for the box means):
#  - land = cells flagged/missing on >= 95% of days; cells within
#    COAST_BUFFER_KM of land are excluded, because passive-microwave pixels
#    next to land read falsely as ice (land spill-over);
#  - an ice cell must have >= MIN_ICE_NEIGHBOURS ice cells among its 8
#    neighbours (isolated ice pixels are treated as noise).
COAST_BUFFER_KM <- 12.5
MIN_ICE_NEIGHBOURS <- 2
land <- rowMeans(is.na(vals)) >= 0.95
land_r <- setValues(asi[[1]], ifelse(land, 1, NA))
coast_km <- values(distance(land_r))[, 1] / 1000          # distance to nearest land cell
edge_ok <- !land & coast_km > COAST_BUFFER_KM
ice_bin <- setValues(asi, (vals >= EDGE_PCT) * 1)
ice_nb <- values(focal(ice_bin, w = matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), 3), fun = "sum",
                       na.rm = TRUE))                     # ice neighbours per cell and day

series <- map_dfr(SITE_ORDER, function(st) {
  b <- boxes[boxes$site == st, ]
  s <- SITES[SITES$site == st, ]
  in_box <- which(ll[, 1] >= b$lon_min & ll[, 1] <= b$lon_max &
                  ll[, 2] >= b$lat_min & ll[, 2] <= b$lat_max)
  dist_km <- as.vector(distance(vect(data.frame(x = s$lon, y = s$lat), geom = c("x", "y"),
                                     crs = "EPSG:4326"), cells_ll)) / 1000
  site_cell <- which.min(dist_km)
  map_dfr(seq_along(dates), function(k) {
    v <- vals[, k]
    vb <- v[in_box]
    here <- v[site_cell]
    ice <- !is.na(v) & v >= EDGE_PCT & edge_ok & ice_nb[, k] >= MIN_ICE_NEIGHBOURS
    open <- !is.na(v) & v < EDGE_PCT & edge_ok
    d_edge <- if (is.na(here)) NA_real_ else if (here < EDGE_PCT) {
      if (any(ice)) min(dist_km[ice]) else NA_real_
    } else if (any(open)) -min(dist_km[open]) else NA_real_
    tibble(site = st, date = dates[k], n_cells = sum(!is.na(vb)),
           ice_conc = mean(vb, na.rm = TRUE) / 100,
           open_water_frac = mean(vb < EDGE_PCT, na.rm = TRUE),
           dist_ice_edge_km = d_edge)
  })
})
series <- mutate(series, across(where(is.double), \(x) round(x, 4)))
write_csv(series, out_path("03_sea_ice", "amsr2_daily_site.csv"))
message("Wrote amsr2_daily_site.csv: ", nrow(series), " rows; days with NA ice: ",
        sum(is.na(series$ice_conc)))
