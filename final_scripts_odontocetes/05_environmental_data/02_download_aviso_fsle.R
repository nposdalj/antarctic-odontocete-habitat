# =============================================================================
# Script  : 02_download_aviso_fsle.R
# Section : 05_environmental_data
# Purpose : Download AVISO+ delayed-time global Finite-Size Lyapunov Exponents
#           (FSLE; backward-in-time, 1/25 degree ~ 0.04 degree, daily;
#           variables fsle_max and theta_max) for the +/- 40 km box around each
#           HARP site (SITES$box40_*; contains the main +/- 20 km box) over
#           the climatology window (deployment start year -
#           CLIM_HALF_WINDOW_YEARS to deployment end year +
#           CLIM_HALF_WINDOW_YEARS, as in 01_download_copernicus.R).
# Inputs  : final_scripts_odontocetes/config.R (SITES boxes and recording
#           times, CLIM_HALF_WINDOW_YEARS, DOWNLOAD_DIR)
#           AVISO+ credentials in the environment variables AVISO_USERNAME and
#           AVISO_PASSWORD (e.g., in ~/.Renviron); internet connection.
# Outputs : DOWNLOAD_DIR/aviso_fsle/
#             {site}_fsle_{YYYY}.nc      (THREDDS NCSS mode: site x year)
#             {site}_fsle_{YYYYMM}.nc    (NCSS monthly fallback for a failed year)
#             {site}_fsle_{YYYYMMDD}.nc  (FTP fallback mode: site x day)
#           Box-cropped NetCDFs with fsle_max and theta_max on (lon, lat,
#           time); read directly by 03_extract_aviso_fsle.R.
#             download_log.txt           (requests, retries, failures)
# Author  : Natalie Posdaljian
# Date    : 2026-11-17
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Requires an AVISO+ account whose subscription includes the
#             "FSLE - Delayed Time - Global" value-added product
#             (https://www.aviso.altimetry.fr). The original analysis used
#             daily global files dt_global_allsat_madt_fsle_YYYYMMDD_20180704.nc
#             downloaded manually.
#           - Credentials are read ONLY from AVISO_USERNAME / AVISO_PASSWORD;
#             they are never written in this script, printed, or logged.
#           - Default method ("ncss"): server-side subsetting with the THREDDS
#             NetCDF Subset Service on the AVISO ODATIS THREDDS server, HTTP
#             basic authentication. The public catalogue (checked 2026-11-17
#             without credentials) lists the aggregated dataset
#             "dataset-duacs-dt-global-allsat-madt-fsle" under
#             L4P > FSLE Data > "FSLE - Delayed Time - Global - FSLE", version
#             label vDT2024, with OPeNDAP and NCSS services; data requests
#             return 401 (basic auth, LDAP realm) without credentials.
#             [CONFIRM] endpoint/dataset path and variable names with a real
#             request (not testable here without credentials).
#             [CONFIRM] the THREDDS dataset is a newer reprocessing (vDT2024)
#             than the original files (production date 20180704); values may
#             differ slightly from the original analysis.
#           - Fallback method ("ftp"): authenticated FTP download of the daily
#             global files, cropped immediately to the site boxes; the global
#             file is deleted after cropping. [CONFIRM] FTP directory.
#           - Resumable: existing output files are skipped. Each file is
#             written to "<name>.part" and renamed when complete, so an
#             interrupted run never leaves a truncated file behind. Requests
#             are retried with back-off and spaced by REQUEST_PAUSE_S seconds.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

library(tidyverse)
library(ncdf4)
library(httr2)
library(curl)

# -----------------------------------------------------------------------------
# 1. Script parameters
# -----------------------------------------------------------------------------
DOWNLOAD_METHOD <- "ncss"     # "ncss" (server-side subset) or "ftp" (fallback)

# [CONFIRM] THREDDS NCSS endpoint for the delayed-time global FSLE aggregation.
NCSS_URL <- paste0("https://tds-odatis.aviso.altimetry.fr/thredds/ncss/grid/",
                   "dataset-duacs-dt-global-allsat-madt-fsle")
# [CONFIRM] FTP folder holding yearly sub-folders of daily global FSLE files.
FTP_BASE_URL <- paste0("ftp://ftp-access.aviso.altimetry.fr/",
                       "value-added/lyapunov/delayed-time/global")
FTP_FILE_REGEX <- "^dt_global_allsat_madt_fsle_(\\d{8})_\\d{8}\\.nc$"

FSLE_VARS       <- c("fsle_max", "theta_max")
MAX_TRIES       <- 5     # attempts per request (with exponential back-off)
REQUEST_PAUSE_S <- 2     # pause between requests, to be polite to the server

FSLE_DOWNLOAD_DIR <- file.path(DOWNLOAD_DIR, "aviso_fsle")
dir.create(FSLE_DOWNLOAD_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(FSLE_DOWNLOAD_DIR, "download_log.txt")

# -----------------------------------------------------------------------------
# 2. Credentials and logging
# -----------------------------------------------------------------------------
aviso_user <- Sys.getenv("AVISO_USERNAME")
aviso_pass <- Sys.getenv("AVISO_PASSWORD")
if (aviso_user == "" || aviso_pass == "") {
  stop("AVISO+ credentials not found: set AVISO_USERNAME and AVISO_PASSWORD ",
       "(e.g., in ~/.Renviron) and restart R.")
}

# Timestamped messages to the console and the log file (never credentials).
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(...))
  message(line)
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

# -----------------------------------------------------------------------------
# 3. Climatology window per site
# -----------------------------------------------------------------------------
# Whole calendar years, identical to 01_download_copernicus.R.
# The DOWNLOAD region is the +/- 40 km box (box40_*), which contains the main
# +/- 20 km analysis box, so one download serves both boxes.
clim_window <- SITES %>%
  transmute(site,
            box_lon_min = box40_lon_min, box_lon_max = box40_lon_max,
            box_lat_min = box40_lat_min, box_lat_max = box40_lat_max,
            start = as.Date(sprintf("%d-01-01",
                                    year(rec_start) - CLIM_HALF_WINDOW_YEARS)),
            end   = as.Date(sprintf("%d-12-31",
                                    year(rec_end) + CLIM_HALF_WINDOW_YEARS)))

# Write to "<file>.part" and rename only once the content is validated.
finalise <- function(tmp, file) {
  nc <- nc_open(tmp)
  ok <- all(FSLE_VARS %in% names(nc$var))
  nc_close(nc)
  if (!ok) {
    unlink(tmp)
    stop("Downloaded file lacks ", paste(FSLE_VARS, collapse = "/"), ": ", basename(file))
  }
  file.rename(tmp, file)
}

# -----------------------------------------------------------------------------
# 4. Method "ncss": one server-side subset per site x year
# -----------------------------------------------------------------------------
# NCSS returns only the box and the requested days, so each request is ~1-2
# MB instead of ~365 global files. Longitudes are sent in 0-360 to match the
# dataset grid. [CONFIRM] that NCSS accepts these bounds as expected.
# One NCSS request for days t0..t1 written to out_file; returns "ok" or the
# error message.
ncss_request <- function(w, t0, t1, out_file) {
  log_msg("NCSS ", w$site, " ", t0, " .. ", t1)
  tmp <- paste0(out_file, ".part")
  req <- request(NCSS_URL) %>%
    req_url_query(var = FSLE_VARS, .multi = "explode") %>%
    req_url_query(north = w$box_lat_max, south = w$box_lat_min,
                  west = w$box_lon_min %% 360, east = w$box_lon_max %% 360,
                  time_start = paste0(t0, "T00:00:00Z"),
                  # 00:00 of the last day: a later end time also returns the
                  # following day (e.g., 1 January of the next year).
                  time_end   = paste0(t1, "T00:00:00Z"),
                  accept = "netcdf") %>%
    req_auth_basic(aviso_user, aviso_pass) %>%
    req_retry(max_tries = MAX_TRIES, backoff = ~ 2^.x) %>%
    req_timeout(600)
  result <- tryCatch({
    req_perform(req, path = tmp)
    finalise(tmp, out_file)
    "ok"
  }, error = function(e) conditionMessage(e))
  if (result != "ok") {
    unlink(tmp)
    log_msg("FAILED ", basename(out_file), ": ", result)
  }
  Sys.sleep(REQUEST_PAUSE_S)
  result
}

# One request per site x year ({site}_fsle_{YYYY}.nc). If a yearly request
# fails (the server has reset the connection on some long requests), the year
# is fetched as monthly files ({site}_fsle_{YYYYMM}.nc) instead, which
# 03_extract_aviso_fsle.R reads the same way.
download_ncss <- function(w) {
  for (yr in year(w$start):year(w$end)) {
    year_file <- file.path(FSLE_DOWNLOAD_DIR, sprintf("%s_fsle_%d.nc", w$site, yr))
    if (file.exists(year_file)) next

    t0 <- max(w$start, as.Date(sprintf("%d-01-01", yr)))
    t1 <- min(w$end, as.Date(sprintf("%d-12-31", yr)))
    months <- seq(floor_date(t0, "month"), t1, by = "month")
    month_files <- file.path(FSLE_DOWNLOAD_DIR,
                             sprintf("%s_fsle_%s.nc", w$site, format(months, "%Y%m")))
    # Resume a year that was already being fetched month by month.
    if (!any(file.exists(month_files)) && ncss_request(w, t0, t1, year_file) == "ok") next

    log_msg("Falling back to monthly requests for ", w$site, " ", yr)
    for (i in seq_along(months)) {
      if (file.exists(month_files[i])) next
      ncss_request(w, max(t0, months[i]), min(t1, ceiling_date(months[i], "month") - 1),
                   month_files[i])
    }
  }
}

# -----------------------------------------------------------------------------
# 5. Method "ftp": daily global files, cropped to every site that needs them
# -----------------------------------------------------------------------------
# Crop one global file to a site box and write a small NetCDF with the same
# variables, coordinates and time axis.
crop_to_site <- function(global_file, w, out_file) {
  nc <- nc_open(global_file)
  on.exit(nc_close(nc))
  lon <- as.vector(ncvar_get(nc, "lon"))
  lat <- as.vector(ncvar_get(nc, "lat"))
  to_file_lon <- function(x) if (max(lon) > 180) x %% 360 else x
  i_lon <- which(lon >= to_file_lon(w$box_lon_min) & lon <= to_file_lon(w$box_lon_max))
  i_lat <- which(lat >= w$box_lat_min & lat <= w$box_lat_max)
  start <- c(min(i_lon), min(i_lat), 1)
  count <- c(length(i_lon), length(i_lat), 1)

  dims <- list(ncdim_def("lon", "degrees_east", lon[i_lon]),
               ncdim_def("lat", "degrees_north", lat[i_lat]),
               ncdim_def("time", ncatt_get(nc, "time", "units")$value,
                         as.vector(ncvar_get(nc, "time")), unlim = TRUE))
  var_units <- function(v) {
    att <- ncatt_get(nc, v, "units")
    if (att$hasatt) att$value else ""
  }
  vars <- map(FSLE_VARS, ~ ncvar_def(.x, var_units(.x), dims, missval = NaN, prec = "float"))
  tmp <- paste0(out_file, ".part")
  out <- nc_create(tmp, vars)
  walk2(vars, FSLE_VARS, ~ ncvar_put(out, .x, ncvar_get(nc, .y, start = start, count = count)))
  nc_close(out)
  finalise(tmp, out_file)
}

download_ftp <- function(windows) {
  ftp_handle <- function(...) {
    new_handle(userpwd = paste0(aviso_user, ":", aviso_pass), use_ssl = 1, ...)
  }
  needed <- windows %>%
    mutate(date = map2(start, end, seq, by = "day")) %>%
    unnest(date) %>%
    mutate(out_file = file.path(FSLE_DOWNLOAD_DIR,
                                sprintf("%s_fsle_%s.nc", site, format(date, "%Y%m%d")))) %>%
    filter(!file.exists(out_file))

  for (yr in sort(unique(year(needed$date)))) {
    # One directory listing per year gives the file name (incl. production
    # date) for every day.
    listing <- curl_fetch_memory(paste0(FTP_BASE_URL, "/", yr, "/"),
                                 handle = ftp_handle(dirlistonly = TRUE))
    files <- tibble(name = str_split_1(rawToChar(listing$content), "\r?\n")) %>%
      filter(str_detect(name, FTP_FILE_REGEX)) %>%
      mutate(date = as.Date(str_match(name, FTP_FILE_REGEX)[, 2], "%Y%m%d"))

    todo <- needed %>% filter(year(date) == yr) %>% left_join(files, by = "date")
    log_msg("FTP ", yr, ": ", n_distinct(todo$date), " day(s) to fetch, ",
            sum(is.na(todo$name)), " site-day(s) without a file on the server")

    for (d in as.list(unique(todo$date[!is.na(todo$name)]))) {
      day <- filter(todo, date == d)
      global_file <- tempfile(fileext = ".nc")
      result <- "failed"
      for (attempt in seq_len(MAX_TRIES)) {
        result <- tryCatch({
          curl_download(paste0(FTP_BASE_URL, "/", yr, "/", day$name[1]), global_file,
                        handle = ftp_handle(), quiet = TRUE)
          for (i in seq_len(nrow(day))) {
            crop_to_site(global_file, filter(windows, site == day$site[i]), day$out_file[i])
          }
          "ok"
        }, error = function(e) conditionMessage(e))
        if (result == "ok") break
        Sys.sleep(2^attempt)
      }
      unlink(global_file)   # never keep the large global file
      if (result != "ok") log_msg("FAILED ", d, ": ", result)
      Sys.sleep(REQUEST_PAUSE_S)
    }
  }
}

# -----------------------------------------------------------------------------
# 6. Run
# -----------------------------------------------------------------------------
log_msg("Start (method = ", DOWNLOAD_METHOD, ")")
if (DOWNLOAD_METHOD == "ncss") {
  walk(split(clim_window, clim_window$site), download_ncss)
} else if (DOWNLOAD_METHOD == "ftp") {
  download_ftp(clim_window)
} else {
  stop("DOWNLOAD_METHOD must be 'ncss' or 'ftp'")
}
log_msg("Done. Files in ", FSLE_DOWNLOAD_DIR, "; failures (if any) are listed in ",
        basename(LOG_FILE), " - re-run the script to retry them.")
