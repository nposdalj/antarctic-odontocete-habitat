# =============================================================================
# Script  : 02_download_modis_sst.R
# Section : Supplementary Fig. S3 (2.5.1) — HYCOM vs Copernicus comparison
# Purpose : Download Aqua MODIS daily sea surface temperature over each site's
#           +/- 40 km box for the site's recording period. MODIS is the
#           observational reference for the surface comparison of HYCOM and
#           Copernicus in 03_compare_hycom_copernicus.R.
# Inputs  : config.R (SITES: box corners and recording periods; DOWNLOAD_DIR)
#           NOAA CoastWatch ERDDAP, dataset erdMH1sstd1dayR20190SQ_Lon0360
#           (Aqua MODIS, NPP, L3SMI, 4 km, daily, daytime, R2019.0),
#           field sstMasked (quality-masked SST, degC)
# Outputs : DOWNLOAD_DIR/modis/modis_sst_<site>.nc
#             sstMasked (time, latitude, longitude); longitude in 0..360
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   :
#   * Why a new download: the original MODIS files (Environmental
#     Data/SST_ERDDAP/) were requested for +/- 0.5 degree longitude windows
#     that do not match the +/- 40 km boxes (the "KGI" request did not contain
#     the KGI site).
#   * Dataset checked on 2026-11-16: still served at
#     https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1sstd1dayR20190SQ_Lon0360
#     (coverage 2003-01-01 to 2022-06-29); a one-month test request took ~5 s.
#   * The dataset uses 0..360 longitudes, so the config box longitudes
#     (-180..180) are shifted by +360.
#   * One request per site (whole recording period, ~1-2 MB each). Existing
#     files are kept; delete a file to download it again.
#   * MODIS is daytime skin SST and is often missing under cloud and sea ice
#     (winter), so daily box coverage is low; this is handled in
#     03_compare_hycom_copernicus.R.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

suppressPackageStartupMessages({
  library(rerddap)
  library(tidyverse)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
MODIS_DIR     <- file.path(DOWNLOAD_DIR, "modis")
ERDDAP_URL    <- "https://coastwatch.pfeg.noaa.gov/erddap/"
MODIS_DATASET <- "erdMH1sstd1dayR20190SQ_Lon0360"
MODIS_FIELD   <- "sstMasked"

dir.create(MODIS_DIR, recursive = TRUE, showWarnings = FALSE)
modis_info <- info(MODIS_DATASET, url = ERDDAP_URL)

# -----------------------------------------------------------------------------
# 2. Download one NetCDF per site
# -----------------------------------------------------------------------------
download_site <- function(code) {
  s        <- site_info(code)
  out_file <- file.path(MODIS_DIR, sprintf("modis_sst_%s.nc", code))
  if (file.exists(out_file)) {
    message(code, ": already on disk")
    return(invisible(out_file))
  }
  message(code, ": requesting ", format(as.Date(s$rec_start)), " to ",
          format(as.Date(s$rec_end)))

  # rerddap caches the response under a hashed name in a scratch folder; the
  # file is then copied to a stable, human-readable name.
  cache <- file.path(tempdir(), "modis_cache")
  res <- griddap(modis_info,
                 latitude  = c(s$box40_lat_min, s$box40_lat_max),
                 longitude = c(s$box40_lon_min, s$box40_lon_max) + 360,
                 time      = format(c(as.Date(s$rec_start), as.Date(s$rec_end))),
                 fields    = MODIS_FIELD,
                 store     = disk(path = cache))
  file.copy(res$summary$filename, out_file, overwrite = TRUE)
  message(sprintf("%s: %d days, %d grid cells per day", code,
                  n_distinct(res$data$time),
                  nrow(res$data) / n_distinct(res$data$time)))
  invisible(out_file)
}

walk(SITE_ORDER, download_site)
