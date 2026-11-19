# =============================================================================
# Script  : 01_download_hycom.R
# Section : Supplementary Fig. S3 (2.5.1) — HYCOM vs Copernicus comparison
# Purpose : Download daily HYCOM GOFS 3.1 (GLBv0.08) water temperature and
#           salinity profiles (0-2000 m) over each site's +/- 40 km box for the
#           site's recording period, so that HYCOM can be compared with the
#           Copernicus reanalysis over exactly the same area and days.
# Inputs  : config.R (SITES: box corners and recording periods; DOWNLOAD_DIR)
#           HYCOM THREDDS OPeNDAP server, https://tds.hycom.org/thredds/dodsC/
# Outputs : DOWNLOAD_DIR/hycom/<site>/hycom_<site>_<YYYY-MM-DD>.nc
#             one file per site and day: water_temp (degC) and salinity (psu)
#             on the native HYCOM grid (lon, lat, depth), one time step
#           DOWNLOAD_DIR/hycom/hycom_download_log.csv
#             site, date, experiment, snapshot time, status for every day
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   :
#   * Why a new download: the HYCOM .mat files used originally
#     (Environmental Data/HYCOM/Data/) were requested for +/- 0.5 degree
#     longitude windows that do not match the +/- 40 km boxes (the "KGI" files
#     span 57.44-56.48 W and do not contain the KGI site at 57.94 W).
#   * Experiments. No single GLBv0.08 experiment spans 2014-2016:
#       expt_53.X  GOFS 3.1 reanalysis   1994-01-01 to 2015-12-31 (yearly files)
#       expt_56.3  GOFS 3.1 analysis     2014-07-01 to 2016-09-30
#       expt_57.2  GOFS 3.1 analysis     2016-05-01 to 2017-01-31
#     HYCOM_SEGMENTS below follows the original download script
#     (Environmental Data/HYCOM/get_hycom.m): EI from expt_53.X, KGI from
#     expt_56.3, CI from expt_56.3 until 2016-09-30 and expt_57.2 afterwards.
#     [CONFIRM which experiment the original .mat files came from — get_hycom.m
#     contains these URLs, but GetHYCOM_sites.m also lists expt_56.3 for EI.]
#     All three endpoints were tested on 2026-11-16 and returned data.
#   * Time sampling: HYCOM output is 3-hourly. One snapshot per day is
#     downloaded — the time step closest to 12:00 UTC on that date — because
#     each 3-hourly step costs ~1 s per variable on the server, and a full daily
#     mean would take ~8x longer (~3 h vs ~25 min). The original .mat files were
#     00:00 UTC snapshots. Copernicus values are daily means. The 12:00 step
#     is missing on some days (88 of 793 site-days in the 2026-11-16 run);
#     the 09:00 step is then used (see the log).
#     [CONFIRM 12:00 UTC snapshot vs daily mean of the 3-hourly steps]
#   * Depth: levels down to 2000 m (36 of 40 levels) so that the deepest
#     species depth (1665 m) is bracketed for vertical interpolation in
#     03_compare_hycom_copernicus.R.
#   * The NetcdfSubset (NCSS) service at ncss.hycom.org was also tested; it
#     took 78 s for one day and returned HTTP 500 for a one-week request, so
#     OPeNDAP via ncdf4 is used.
#   * The server is shared and occasionally unstable: one request per variable
#     and day, retries with increasing waits, and days already on disk are
#     skipped, so the script can simply be re-run after an interruption.
#   * Usage: Rscript 01_download_hycom.R [site] [first_date] [last_date]
#     With no arguments all sites and full recording periods are downloaded;
#     the optional arguments restrict the run (e.g. for a quick test).
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

suppressPackageStartupMessages({
  library(ncdf4)
  library(tidyverse)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
HYCOM_DIR         <- file.path(DOWNLOAD_DIR, "hycom")
HYCOM_OPENDAP     <- "https://tds.hycom.org/thredds/dodsC/GLBv0.08/"
HYCOM_VARS        <- c("water_temp", "salinity")
HYCOM_MAX_DEPTH_M <- 2000   # deepest level kept (brackets 1665 m)
SNAPSHOT_HOUR_UTC <- 12     # daily snapshot: time step closest to this hour
MAX_TRIES         <- 5      # attempts per request before giving up on a day
RETRY_WAIT_S      <- 20     # wait before retry i is i * RETRY_WAIT_S seconds

# Experiment used for each site and date range (see Notes). Dates are
# inclusive; NA means open-ended.
HYCOM_SEGMENTS <- tribble(
  ~site, ~from,        ~to,          ~experiment, ~dataset,
  "EI",  NA,           NA,           "expt_53.X", "expt_53.X/data/2014",
  "KGI", NA,           NA,           "expt_56.3", "expt_56.3",
  "CI",  NA,           "2016-09-30", "expt_56.3", "expt_56.3",
  "CI",  "2016-10-01", NA,           "expt_57.2", "expt_57.2"
) %>%
  mutate(from = as.Date(from), to = as.Date(to))

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
# Run f(); on error wait and try again. HYCOM's server drops connections now
# and then, so a transient failure should not abort the whole download.
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

# Open an OPeNDAP dataset once and keep its coordinate vectors, so each daily
# request only transfers the data subset.
open_dataset <- function(dataset) {
  nc   <- nc_open(paste0(HYCOM_OPENDAP, dataset))
  time <- as.POSIXct("2000-01-01", tz = "UTC") +
    ncvar_get(nc, "time") * 3600            # units: hours since 2000-01-01
  stopifnot(grepl("hours since 2000-01-01", ncatt_get(nc, "time", "units")$value))
  list(nc    = nc,
       lon   = as.vector(ncvar_get(nc, "lon")),
       lat   = as.vector(ncvar_get(nc, "lat")),
       depth = as.vector(ncvar_get(nc, "depth")),
       time  = time)
}

# Grid indices of the cells whose centres fall inside a site's box. HYCOM
# longitudes are -180..180 for these GLBv0.08 datasets but 0..360 for some
# other HYCOM products, so both conventions are handled.
box_indices <- function(ds, s) {
  lon <- ds$lon
  if (max(lon) > 180) lon <- ifelse(lon > 180, lon - 360, lon)
  list(i_lon   = which(lon >= s$box40_lon_min & lon <= s$box40_lon_max),
       i_lat   = which(ds$lat >= s$box40_lat_min & ds$lat <= s$box40_lat_max),
       i_depth = which(ds$depth <= HYCOM_MAX_DEPTH_M),
       lon     = lon)
}

# Index of the time step on `date` closest to SNAPSHOT_HOUR_UTC (NA if the
# dataset has no output on that date; e.g. expt_53.X has occasional gaps).
snapshot_index <- function(ds, date) {
  target  <- as.POSIXct(date, tz = "UTC") + SNAPSHOT_HOUR_UTC * 3600
  on_date <- which(as.Date(ds$time) == date)
  if (length(on_date) == 0) return(NA_integer_)
  on_date[which.min(abs(difftime(ds$time[on_date], target, units = "hours")))]
}

# Write one day of box data to a compressed NetCDF-4 file.
write_day <- function(file, idx, ds, values, snap_time, experiment, s) {
  dims <- list(
    ncdim_def("lon",   "degrees_east",  idx$lon[idx$i_lon]),
    ncdim_def("lat",   "degrees_north", ds$lat[idx$i_lat]),
    ncdim_def("depth", "m",             ds$depth[idx$i_depth]),
    ncdim_def("time",  "hours since 2000-01-01 00:00:00",
              as.numeric(difftime(snap_time, as.POSIXct("2000-01-01", tz = "UTC"),
                                  units = "hours")), unlim = TRUE)
  )
  vars <- list(
    ncvar_def("water_temp", "degC", dims, missval = -999, prec = "float",
              longname = "sea water potential temperature", compression = 4),
    ncvar_def("salinity",   "psu",  dims, missval = -999, prec = "float",
              longname = "sea water salinity", compression = 4)
  )
  tmp <- paste0(file, ".part")      # write to a temporary name, then rename,
  nc  <- nc_create(tmp, vars, force_v4 = TRUE)  # so partial files are never kept
  for (v in HYCOM_VARS) ncvar_put(nc, v, values[[v]])
  ncatt_put(nc, 0, "source", paste0(HYCOM_OPENDAP, experiment))
  ncatt_put(nc, 0, "experiment", experiment)
  ncatt_put(nc, 0, "site", s$site)
  ncatt_put(nc, 0, "box", sprintf("lat %.5f to %.5f, lon %.5f to %.5f",
                                  s$box40_lat_min, s$box40_lat_max,
                                  s$box40_lon_min, s$box40_lon_max))
  ncatt_put(nc, 0, "sampling",
            sprintf("single 3-hourly snapshot closest to %02d:00 UTC", SNAPSHOT_HOUR_UTC))
  nc_close(nc)
  file.rename(tmp, file)
}

# Download one site-day; returns a one-row log tibble.
download_day <- function(s, date, seg, ds) {
  out_file <- file.path(HYCOM_DIR, s$site,
                        sprintf("hycom_%s_%s.nc", s$site, format(date)))
  log_row <- function(status, snap = as.POSIXct(NA, tz = "UTC")) {
    tibble(site = s$site, date = date, experiment = seg$experiment,
           snapshot_utc = format(snap, "%Y-%m-%d %H:%M"), status = status)
  }
  if (file.exists(out_file)) {
    # Re-read the snapshot time so the log stays complete across re-runs.
    nc   <- nc_open(out_file)
    snap <- as.POSIXct("2000-01-01", tz = "UTC") + ncvar_get(nc, "time") * 3600
    nc_close(nc)
    return(log_row("already on disk", snap))
  }

  it <- snapshot_index(ds, date)
  if (is.na(it)) return(log_row("no HYCOM output on this date"))

  idx    <- box_indices(ds, s)
  start  <- c(min(idx$i_lon), min(idx$i_lat), min(idx$i_depth), it)
  count  <- c(length(idx$i_lon), length(idx$i_lat), length(idx$i_depth), 1)
  values <- map(set_names(HYCOM_VARS), function(v) {
    with_retries(function() ncvar_get(ds$nc, v, start = start, count = count,
                                      collapse_degen = FALSE),
                 sprintf("%s %s %s", s$site, format(date), v))
  })
  if (any(map_lgl(values, is.null))) return(log_row("FAILED", ds$time[it]))

  write_day(out_file, idx, ds, values, ds$time[it], seg$experiment, s)
  log_row("downloaded", ds$time[it])
}

# -----------------------------------------------------------------------------
# 3. Download
# -----------------------------------------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
run_sites  <- if (length(args) >= 1) args[1] else SITE_ORDER
date_first <- if (length(args) >= 2) as.Date(args[2]) else as.Date(-Inf)
date_last  <- if (length(args) >= 3) as.Date(args[3]) else as.Date(Inf)

# One row per site-day of the recording period, with its HYCOM experiment.
# (as.Date() on rec_start/rec_end gives the UTC calendar days with recording.)
plan <- SITES %>%
  filter(site %in% run_sites) %>%
  rowwise() %>%
  reframe(site, date = seq(as.Date(rec_start), as.Date(rec_end), by = "day")) %>%
  filter(date >= date_first, date <= date_last) %>%
  inner_join(HYCOM_SEGMENTS, by = "site", relationship = "many-to-many") %>%
  filter(is.na(from) | date >= from, is.na(to) | date <= to)

# Every requested site-day must map to exactly one experiment.
stopifnot(!anyDuplicated(plan[c("site", "date")]))

t_start <- Sys.time()
log <- list()
for (seg_id in unique(plan$dataset)) {
  seg_days <- filter(plan, dataset == seg_id)
  message(sprintf("Opening %s (%d site-days)", seg_id, nrow(seg_days)))
  ds <- with_retries(function() open_dataset(seg_id), seg_id)
  if (is.null(ds)) stop("Could not open HYCOM dataset ", seg_id)

  for (k in seq_len(nrow(seg_days))) {
    s   <- site_info(seg_days$site[k])
    dir.create(file.path(HYCOM_DIR, s$site), recursive = TRUE, showWarnings = FALSE)
    log[[length(log) + 1]] <- download_day(s, seg_days$date[k], seg_days[k, ], ds)
    if (k %% 25 == 0) message(sprintf("  %s: %d/%d days", seg_id, k, nrow(seg_days)))
  }
  nc_close(ds$nc)
}

# -----------------------------------------------------------------------------
# 4. Log
# -----------------------------------------------------------------------------
# Append to the existing log (keeping the latest entry per site-day) so that
# repeated or partial runs build up one complete record.
log <- bind_rows(log)
log_file <- file.path(HYCOM_DIR, "hycom_download_log.csv")
if (file.exists(log_file)) {
  old <- read_csv(log_file, col_types = cols(date = col_date(), .default = col_character()))
  log <- bind_rows(anti_join(old, log, by = c("site", "date")), log)
}
log <- arrange(log, factor(site, SITE_ORDER), date)
write_csv(log, log_file)

message(sprintf("Done in %.1f min.", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
print(count(log, site, experiment, status))
if (any(log$status == "FAILED")) {
  warning("Some days failed; re-run the script to retry them.")
}
