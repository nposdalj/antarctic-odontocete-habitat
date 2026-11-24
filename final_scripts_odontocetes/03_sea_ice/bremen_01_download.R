# =============================================================================
# Script  : bremen_01_download.R
# Section : 2.3 Sea ice - satellite sea-ice concentration (comparison, maps)
# Purpose : Download daily AMSR2 ASI sea-ice concentration (University of
#           Bremen, version 5.4, 6.25 km polar stereographic grid, EPSG:3412)
#           for 2013-2021 and keep only a regional crop around the South
#           Shetland Islands, saved as one multi-layer GeoTIFF per year.
# Inputs  : https://data.seaice.uni-bremen.de/amsr2/asi_daygrid_swath/s6250/
#             <year>/<mon>/Antarctic/asi-AMSR2-s6250-<yyyymmdd>-v5.4.tif
#           (open access; Spreen et al. 2008 [CONFIRM reference])
# Outputs : DOWNLOAD_DIR/bremen_asi/asi_s6250_<year>_region.tif (one layer per
#           day, layer names = dates; missing days are skipped)
#           DOWNLOAD_DIR/bremen_asi/missing_dates.csv
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Daily files are deleted after cropping. Values: 0-100 % ice;
#           values above 100 are flags (land/missing) and are set to NA in
#           bremen_02_compare.R.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages({ library(tidyverse); library(terra) })

YEARS <- 2013:2021   # recording years 2014-2016 plus context (AMSR2 starts July 2012)
REGION_LONLAT <- ext(-66, -48, -65, -58)   # lon/lat box for maps and site boxes
BASE_URL <- "https://data.seaice.uni-bremen.de/amsr2/asi_daygrid_swath/s6250"
out_dir <- file.path(DOWNLOAD_DIR, "bremen_asi")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 120, HTTPUserAgent = "research-client")

region_3412 <- NULL
missing <- list()
for (yr in YEARS) {
  out_file <- file.path(out_dir, paste0("asi_s6250_", yr, "_region.tif"))
  if (file.exists(out_file)) { message(yr, " already done"); next }
  dates <- seq(as.Date(paste0(yr, "-01-01")), as.Date(paste0(yr, "-12-31")), by = "day")
  layers <- list()
  for (d in as.list(dates)) {
    tag <- format(d, "%Y%m%d")
    url <- sprintf("%s/%s/%s/Antarctic/asi-AMSR2-s6250-%s-v5.4.tif", BASE_URL, yr,
                   tolower(format(d, "%b")), tag)
    tmp <- file.path(out_dir, paste0("tmp_", tag, ".tif"))
    ok <- tryCatch({ download.file(url, tmp, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok || !file.exists(tmp)) { missing[[tag]] <- tibble(date = d); next }
    r <- rast(tmp)
    if (is.null(region_3412)) region_3412 <- ext(project(as.polygons(REGION_LONLAT, crs = "EPSG:4326"), crs(r)))
    layers[[tag]] <- setNames(crop(r, region_3412) * 1, as.character(d))   # *1 drops the colour table
    unlink(tmp)
  }
  writeRaster(rast(layers), out_file, overwrite = TRUE)
  message(yr, ": ", length(layers), " days")
}
write_csv(bind_rows(missing), file.path(out_dir, "missing_dates.csv"))
