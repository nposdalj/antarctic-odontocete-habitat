# =============================================================================
# Script  : 05_download_models_at_profiles.R
# Section : Supplementary (2.5.1) - model validation against in-situ profiles
# Purpose : Extract the Copernicus and HYCOM model profiles that match each
#           in-situ cast from 04_download_insitu_profiles.R: for every cast,
#           the model water column at the nearest ocean grid cell on the
#           cast's day (Copernicus daily mean; HYCOM 3-hourly snapshot closest
#           to the cast time). Temperature and salinity from both models;
#           dissolved oxygen from the Copernicus biogeochemistry hindcast only
#           (HYCOM has no biogeochemistry).
# Inputs  : outputs/05b_hycom_comparison/insitu_profiles.csv
#           config.R (COPERNICUS products/layers, DOWNLOAD_DIR)
#           05_environmental_data/copernicus_arco_helpers.R
#           Copernicus Marine credentials (environment variables
#             COPERNICUSMARINE_SERVICE_USERNAME / _PASSWORD)
#           HYCOM THREDDS OPeNDAP server, https://tds.hycom.org/thredds/dodsC/
# Outputs : DOWNLOAD_DIR/models_at_profiles/<model>/<YYYY-MM-DD>.rds (cache)
#           outputs/05b_hycom_comparison/models_at_profiles.csv
#             cast_id, model, variable, model_depth_m, value, model_lon,
#             model_lat, distance_km, model_time
# Author  : Natalie Posdaljian
# Date    : 2026-11-19
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Matching cell: the grid cell nearest to the cast position that
#             has a value at the model surface (i.e., is ocean in the model).
#             Its distance from the cast is saved (distance_km).
#           - Copernicus: one request per day and product, covering the
#             positions of that day's casts plus a margin, single time step
#             (multi-step requests can come back scrambled with
#             CopernicusMarine 0.3.2; see copernicus_arco_helpers.R).
#           - HYCOM experiments by date, as in 01_download_hycom.R for the
#             recording periods: expt_53.X (reanalysis) for 2014, expt_56.3
#             (analysis) from 2015-01-01 to 2016-09-30, expt_57.2 (analysis)
#             from 2016-10-01.
#           - Depth range: surface to MAX_DEPTH_M (covers the deepest species
#             depth, 1665 m).
#           - Resumable: days already cached are skipped.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "05_environmental_data", "copernicus_arco_helpers.R"))

suppressPackageStartupMessages({
  library(ncdf4)
  library(tidyverse)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
SECTION      <- "05b_hycom_comparison"
CACHE_DIR    <- file.path(DOWNLOAD_DIR, "models_at_profiles")
MAX_DEPTH_M  <- 1750
MARGIN_DEG   <- c(phy = 0.2, bgc = 0.5)   # added around the day's casts
GRID_RES_DEG <- c(phy = 1 / 12, bgc = 0.25)
COP_VARS     <- list(phy = c(temperature = "thetao", salinity = "so"),
                     bgc = c(oxygen = "o2"))
PLAUSIBLE    <- list(thetao = c(-3, 15), so = c(30, 36), o2 = c(100, 500))

HYCOM_OPENDAP <- "https://tds.hycom.org/thredds/dodsC/GLBv0.08/"
HYCOM_VARS    <- c(temperature = "water_temp", salinity = "salinity")
HYCOM_SEGMENTS <- tribble(
  ~from,        ~to,          ~dataset,
  "2014-01-01", "2014-12-31", "expt_53.X/data/2014",
  "2015-01-01", "2016-09-30", "expt_56.3",
  "2016-10-01", "2016-12-31", "expt_57.2"
) |>
  mutate(from = as.Date(from), to = as.Date(to))
HYCOM_HALF_WINDOW <- 3      # grid cells around a cast read from HYCOM
MAX_TRIES    <- 5
RETRY_WAIT_S <- 20

# -----------------------------------------------------------------------------
# 2. Casts
# -----------------------------------------------------------------------------
casts <- read_csv(out_path(SECTION, "insitu_profiles.csv"),
                  show_col_types = FALSE,
                  col_types = cols(platform = col_character(),
                                   .default = col_guess())) |>
  distinct(cast_id, time, lat, lon) |>
  mutate(date = as.Date(time, tz = "UTC"))
message(nrow(casts), " casts on ", n_distinct(casts$date), " days")

# -----------------------------------------------------------------------------
# 3. Helpers
# -----------------------------------------------------------------------------
with_retries <- function(f, what) {
  for (i in seq_len(MAX_TRIES)) {
    out <- tryCatch(f(), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    message(sprintf("  %s: attempt %d/%d failed (%s)", what, i, MAX_TRIES,
                    conditionMessage(out)))
    if (i < MAX_TRIES) Sys.sleep(i * RETRY_WAIT_S)
  }
  NULL
}

# Index (i, j) of the cell nearest to (lon0, lat0) among cells where `ocean`
# is TRUE; distance in km on a local flat approximation.
nearest_ocean_cell <- function(lon, lat, ocean, lon0, lat0) {
  grid <- expand_grid(j = seq_along(lat), i = seq_along(lon)) |>
    mutate(ok = ocean[cbind(i, j)],
           d = sqrt(((lon[i] - lon0) * cos(lat0 * pi / 180))^2 +
                      (lat[j] - lat0)^2) * 111.32) |>
    filter(ok)
  if (nrow(grid) == 0) return(NULL)
  grid[which.min(grid$d), ]
}

# Long table of one cast's model profile from a [lon, lat, depth] array.
cast_profile <- function(arr, lon, lat, depth, cast, model, variable,
                         model_time) {
  ocean <- !is.na(arr[, , 1])
  cell <- nearest_ocean_cell(lon, lat, ocean, cast$lon, cast$lat)
  if (is.null(cell)) return(NULL)
  tibble(cast_id = cast$cast_id, model = model, variable = variable,
         model_depth_m = depth, value = arr[cell$i, cell$j, ],
         model_lon = lon[cell$i], model_lat = lat[cell$j],
         distance_km = cell$d, model_time = model_time) |>
    filter(!is.na(value))
}

# -----------------------------------------------------------------------------
# 4. Copernicus (one request per day and product)
# -----------------------------------------------------------------------------
copernicus_day <- function(day_casts, attrs) {
  date <- day_casts$date[1]
  t_day <- as.POSIXct(rep(paste(format(date), "00:00:00"), 2), tz = "UTC")
  map(names(COP_VARS), function(p) {
    m <- MARGIN_DEG[[p]]
    region <- c(min(day_casts$lon) - m, min(day_casts$lat) - m,
                max(day_casts$lon) + m, max(day_casts$lat) + m)
    imap(COP_VARS[[p]], function(v, name) {
      x <- with_retries(function() {
        download_subset(COPERNICUS[[p]]$product, COPERNICUS[[p]]$layer, v,
                        region, t_day, verticalrange = c(0, -MAX_DEPTH_M))
      }, sprintf("Copernicus %s %s", format(date), v))
      if (is.null(x)) stop("Copernicus download failed: ", format(date), " ", v)
      a <- subset_to_array(x, v, attrs[[p]][[v]], GRID_RES_DEG[[p]],
                           plausible_range = PLAUSIBLE[[v]])
      map(seq_len(nrow(day_casts)), function(k) {
        cast_profile(a$values[, , , 1], a$lon, a$lat, a$depth,
                     day_casts[k, ], "Copernicus", name, a$time[1])
      }) |>
        list_rbind()
    }) |>
      list_rbind()
  }) |>
    list_rbind()
}

# -----------------------------------------------------------------------------
# 5. HYCOM (one request per cast: small window, snapshot nearest the cast)
# -----------------------------------------------------------------------------
open_hycom <- function(dataset) {
  nc <- nc_open(paste0(HYCOM_OPENDAP, dataset))
  stopifnot(grepl("hours since 2000-01-01",
                  ncatt_get(nc, "time", "units")$value))
  lon <- as.vector(ncvar_get(nc, "lon"))
  list(nc = nc,
       lon = ifelse(lon > 180, lon - 360, lon),
       lat = as.vector(ncvar_get(nc, "lat")),
       depth = as.vector(ncvar_get(nc, "depth")),
       time = as.POSIXct("2000-01-01", tz = "UTC") +
         as.vector(ncvar_get(nc, "time")) * 3600)
}

hycom_cast <- function(ds, cast) {
  it <- which.min(abs(as.numeric(difftime(ds$time, cast$time, units = "hours"))))
  if (abs(as.numeric(difftime(ds$time[it], cast$time, units = "hours"))) > 3) {
    return(NULL)                     # no HYCOM output within 3 h of the cast
  }
  i0 <- which.min(abs(ds$lon - cast$lon))
  j0 <- which.min(abs(ds$lat - cast$lat))
  ii <- max(1, i0 - HYCOM_HALF_WINDOW):min(length(ds$lon), i0 + HYCOM_HALF_WINDOW)
  jj <- max(1, j0 - HYCOM_HALF_WINDOW):min(length(ds$lat), j0 + HYCOM_HALF_WINDOW)
  kk <- which(ds$depth <= MAX_DEPTH_M + 250)  # one level below MAX_DEPTH_M
  start <- c(min(ii), min(jj), min(kk), it)
  count <- c(length(ii), length(jj), length(kk), 1)
  imap(HYCOM_VARS, function(v, name) {
    a <- with_retries(function() {
      ncvar_get(ds$nc, v, start = start, count = count, collapse_degen = FALSE)
    }, sprintf("HYCOM cast %s %s", cast$cast_id, v))
    if (is.null(a)) stop("HYCOM download failed for cast ", cast$cast_id)
    cast_profile(a[, , , 1], ds$lon[ii], ds$lat[jj], ds$depth[kk], cast,
                 "HYCOM", name, ds$time[it])
  }) |>
    list_rbind()
}

# -----------------------------------------------------------------------------
# 6. Run (cached per model and day)
# -----------------------------------------------------------------------------
for (model in c("Copernicus", "HYCOM")) {
  dir.create(file.path(CACHE_DIR, model), recursive = TRUE, showWarnings = FALSE)
}

days <- sort(unique(casts$date))

# Copernicus
check_copernicus_setup()
cop_attrs <- imap(COP_VARS, \(v, p) layer_zarr_attrs(
  COPERNICUS[[p]]$product, COPERNICUS[[p]]$layer, unname(v)))
for (k in seq_along(days)) {
  cache <- file.path(CACHE_DIR, "Copernicus", paste0(format(days[k]), ".rds"))
  if (file.exists(cache)) next
  saveRDS(copernicus_day(filter(casts, date == days[k]), cop_attrs), cache)
  if (k %% 20 == 0) message(sprintf("Copernicus: %d/%d days", k, length(days)))
}

# HYCOM (one dataset opened per segment)
for (s in seq_len(nrow(HYCOM_SEGMENTS))) {
  seg <- HYCOM_SEGMENTS[s, ]
  seg_days <- days[days >= seg$from & days <= seg$to]
  todo <- seg_days[!file.exists(file.path(CACHE_DIR, "HYCOM",
                                          paste0(format(seg_days), ".rds")))]
  if (length(todo) == 0) next
  message(sprintf("HYCOM %s: %d days", seg$dataset, length(todo)))
  ds <- with_retries(function() open_hycom(seg$dataset), seg$dataset)
  if (is.null(ds)) stop("Could not open HYCOM dataset ", seg$dataset)
  for (k in seq_along(todo)) {
    day_casts <- filter(casts, date == todo[k])
    out <- map(seq_len(nrow(day_casts)), \(i) hycom_cast(ds, day_casts[i, ])) |>
      list_rbind()
    saveRDS(out, file.path(CACHE_DIR, "HYCOM", paste0(format(todo[k]), ".rds")))
    if (k %% 20 == 0) message(sprintf("  %d/%d days", k, length(todo)))
  }
  nc_close(ds$nc)
}

# -----------------------------------------------------------------------------
# 7. Combine
# -----------------------------------------------------------------------------
models_at_profiles <- list.files(CACHE_DIR, pattern = "\\.rds$",
                                 recursive = TRUE, full.names = TRUE) |>
  map(readRDS) |>
  list_rbind()
write_csv(models_at_profiles, out_path(SECTION, "models_at_profiles.csv"))

models_at_profiles |>
  group_by(model, variable) |>
  summarise(casts = n_distinct(cast_id),
            median_distance_km = median(distance_km), .groups = "drop") |>
  print()
