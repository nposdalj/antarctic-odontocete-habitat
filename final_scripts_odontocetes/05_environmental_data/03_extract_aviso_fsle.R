# =============================================================================
# Script  : 03_extract_aviso_fsle.R
# Section : 05_environmental_data
# Purpose : Summarise AVISO+ Finite-Size Lyapunov Exponents (FSLE) over two
#           boxes around each HARP site - the MAIN +/- 20 km analysis box
#           (SITES$box_*) and the +/- 40 km comparison box (SITES$box40_*):
#           daily box mean and SD of FSLE magnitude (fsle_max) and daily axial
#           circular mean and SD of the FSLE orientation (theta_max).
# Inputs  : final_scripts_odontocetes/config.R (SITES boxes and recording
#           times, CLIM_HALF_WINDOW_YEARS, DOWNLOAD_DIR)
#           NetCDFs in FSLE_DIR (default DOWNLOAD_DIR/aviso_fsle/, written by
#           02_download_aviso_fsle.R), in either layout:
#             {site}_fsle_{YYYY|YYYYMM|YYYYMMDD}.nc (box-cropped, year,
#               month or day)
#             dt_global_allsat_madt_fsle_YYYYMMDD_<prod. date>.nc (global
#               daily files as distributed by AVISO; the original analysis
#               used production date 20180704)
#           Both contain fsle_max and theta_max on (lon, lat, time).
# Outputs : F_FSLE_DAILY (outputs/05_environmental_data/fsle_daily.csv; main
#           box) and fsle_daily_box40.csv (+/- 40 km comparison box):
#           date, site, FSLE, FSLE_sd, fsle_orient, fsle_orient_sd,
#           frac_cells_with_fsle
# Author  : Natalie Posdaljian
# Date    : 2026-11-18
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Product: AVISO+ delayed-time global FSLE, computed from DUACS
#             all-satellite absolute geostrophic velocities (MADT), 1/25 degree
#             (~0.04 degree), daily, backward-in-time integration.
#           - FSLE magnitude (fsle_max, day^-1): backward-in-time exponents are
#             stored as NEGATIVE values, more negative values representing more
#             intense and stronger convergence (attracting Lagrangian
#             fronts). FSLE and FSLE_sd are the daily box mean and SD.
#           - Zeros (decided by N. Posdaljian, 2026-11-17): fsle_max is exactly
#             0 in ocean cells where the separation threshold was not reached
#             within the integration time. These are real values, not missing
#             data: the smallest nonzero magnitude is -0.0171 day^-1, the fill
#             value (9.97e36) is unused, and no cell is always 0 (e.g., 17.6%
#             of CI cells in 2011, 1-70% of the box per day). Therefore:
#             FSLE and FSLE_sd INCLUDE zeros (FSLE_INCLUDE_ZEROS = TRUE, as in
#             the original analysis); orientation is summarised ONLY over cells
#             with fsle_max < 0 (FSLE_ORIENT_CELLS = "negative"), because
#             theta_max is noise where FSLE = 0 (it spans the full -180..180
#             there); frac_cells_with_fsle = daily fraction of valid box cells
#             with fsle_max < 0. Days with no cell < 0 get NA orientation and
#             are counted in the console report.
#           - FSLE orientation (theta_max, degrees): orientation of the
#             eigenvector associated with the maximum FSLE, measured
#             anticlockwise from East. Native range in the DT2024 files is
#             -180..180 degrees (confirmed on downloaded files), but it is an
#             AXIAL quantity: theta and theta +/- 180 degrees describe the same
#             line, so axial statistics are used (invariant to +/-180).
#             Default summary (FSLE_ORIENT_METHOD = "axial"): angles are
#             doubled, averaged as unit vectors, and halved back; the result
#             is in (-90, 90]. The SD is the axial circular SD,
#             sqrt(-2 ln R2) / 2 in degrees, where R2 is the mean resultant
#             length of the doubled angles.
#             Why not the arithmetic mean: near-vertical orientations (close
#             to +90 and -90 degrees, i.e., the same line) cancel to ~0 under
#             arithmetic averaging. The original daily orientation means were
#             confined to about -27..16 degrees, consistent with that
#             cancellation. "arithmetic" is kept only as an option that
#             reproduces the original analysis.
#           - Fixes a bug in the original GetAVISO.R, where the orientation SD
#             was a copy of the FSLE SD.
#           - Time window: the climatology window (deployment start year -
#             CLIM_HALF_WINDOW_YEARS to deployment end year +
#             CLIM_HALF_WINDOW_YEARS). Days without data are counted and
#             reported; missing recording-period days trigger a warning.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

library(tidyverse)
library(ncdf4)

# -----------------------------------------------------------------------------
# 1. Script parameters
# -----------------------------------------------------------------------------
# Folder with the FSLE NetCDFs (searched recursively). Defaults to the output
# folder of 02_download_aviso_fsle.R; the environment variable AVISO_FSLE_DIR
# overrides it (e.g., a folder of global daily files such as the original
# D:/FSLE/AVISO_Global_FSLE/).
FSLE_DIR <- Sys.getenv("AVISO_FSLE_DIR", unset = file.path(DOWNLOAD_DIR, "aviso_fsle"))

SITE_FILE_REGEX   <- "^(EI|KGI|CI)_fsle_(\\d{4}|\\d{6}|\\d{8})\\.nc$"

# Boxes to summarise (corners per site) and the suffix of their output file.
BOXES <- list(
  main  = list(suffix = "",
               corners = transmute(SITES, site, lat_min = box_lat_min, lat_max = box_lat_max,
                                   lon_min = box_lon_min, lon_max = box_lon_max)),
  box40 = list(suffix = "_box40",
               corners = transmute(SITES, site, lat_min = box40_lat_min, lat_max = box40_lat_max,
                                   lon_min = box40_lon_min, lon_max = box40_lon_max))
)
with_suffix <- function(file, suffix) sub("\\.csv$", paste0(suffix, ".csv"), file)
GLOBAL_FILE_REGEX <- "^dt_global_allsat_madt_fsle_(\\d{8})_\\d{8}\\.nc$"

# Orientation summary: "axial" (default; circular statistics modulo 180
# degrees) or "arithmetic" (reproduces the original analysis only).
FSLE_ORIENT_METHOD <- "axial"

# Zeros in fsle_max (threshold not reached): include them in the FSLE mean/SD
# (TRUE, original analysis) or use only cells with fsle_max < 0 (FALSE).
FSLE_INCLUDE_ZEROS <- TRUE

# Cells used for the orientation summary: "negative" (only cells with
# fsle_max < 0; default) or "all" (every valid cell).
FSLE_ORIENT_CELLS <- "negative"

if (!dir.exists(FSLE_DIR)) {
  stop("FSLE_DIR does not exist: '", FSLE_DIR, "'. Run 02_download_aviso_fsle.R ",
       "or set the AVISO_FSLE_DIR environment variable.")
}

# -----------------------------------------------------------------------------
# 2. Climatology window per site
# -----------------------------------------------------------------------------
# Whole calendar years around each deployment, as in 01_download_copernicus.R.
clim_window <- SITES %>%
  transmute(site,
            start = as.Date(sprintf("%d-01-01",
                                    year(rec_start) - CLIM_HALF_WINDOW_YEARS)),
            end   = as.Date(sprintf("%d-12-31",
                                    year(rec_end) + CLIM_HALF_WINDOW_YEARS)))

# -----------------------------------------------------------------------------
# 3. Index the files and decide which sites each one serves
# -----------------------------------------------------------------------------
# Box-cropped files serve their own site; global daily files serve every site
# whose window contains that day (global files outside all windows are not
# opened at all).
all_files <- tibble(path = list.files(FSLE_DIR, pattern = "\\.nc$",
                                      recursive = TRUE, full.names = TRUE)) %>%
  mutate(name = basename(path))

# `period` (YYYY, YYYYMM or YYYYMMDD) is the span a box-cropped file is meant to
# cover; days outside it are discarded after reading (yearly NCSS files can
# include 1 January of the following year).
site_files <- all_files %>%
  filter(str_detect(name, SITE_FILE_REGEX)) %>%
  mutate(site   = str_match(name, SITE_FILE_REGEX)[, 2],
         period = str_match(name, SITE_FILE_REGEX)[, 3])

global_files <- all_files %>%
  filter(str_detect(name, GLOBAL_FILE_REGEX)) %>%
  mutate(date = as.Date(str_match(name, GLOBAL_FILE_REGEX)[, 2], "%Y%m%d")) %>%
  inner_join(clim_window, by = join_by(between(date, start, end))) %>%
  select(path, name, site)

jobs <- bind_rows(site_files, global_files) %>%
  group_by(path) %>%
  summarise(sites = list(site), period = first(period), .groups = "drop")
if (nrow(jobs) == 0) stop("No FSLE NetCDFs found in ", FSLE_DIR)
message(sprintf("%d box-cropped and %d global FSLE file(s) to read",
                nrow(site_files), n_distinct(global_files$path)))

# -----------------------------------------------------------------------------
# 4. Box summary helpers
# -----------------------------------------------------------------------------
# Axial (modulo-180 degree) statistics: double the angles, take the mean unit
# vector, halve its direction; SD = sqrt(-2 ln R2) / 2 with R2 the mean
# resultant length of the doubled angles (Mardia & Jupp 2000).
axial_stats <- function(theta_deg) {
  th2 <- 2 * theta_deg * pi / 180
  s <- mean(sin(th2))
  c2 <- mean(cos(th2))
  r2 <- sqrt(s^2 + c2^2)
  c(mean = atan2(s, c2) / 2 * 180 / pi,
    sd   = sqrt(-2 * log(r2)) / 2 * 180 / pi)
}

summarise_box <- function(fsle, theta) {
  ok     <- !is.na(fsle)
  is_neg <- ok & fsle < 0
  # Cells entering the magnitude and orientation summaries (see parameters).
  mag   <- fsle[if (FSLE_INCLUDE_ZEROS) ok else is_neg]
  th    <- theta[(if (FSLE_ORIENT_CELLS == "negative") is_neg else ok) & !is.na(theta)]
  orient <- if (length(th) == 0) {
    c(mean = NA_real_, sd = NA_real_)
  } else if (FSLE_ORIENT_METHOD == "axial") {
    axial_stats(th)
  } else {
    c(mean = mean(th), sd = sd(th))
  }
  tibble(FSLE                 = if (length(mag) > 0) mean(mag) else NA_real_,
         FSLE_sd              = if (length(mag) > 1) sd(mag) else NA_real_,
         fsle_orient          = orient[["mean"]],
         fsle_orient_sd       = orient[["sd"]],   # SD of theta_max (bug fix)
         frac_cells_with_fsle = sum(is_neg) / sum(ok),
         n_cells              = sum(ok))
}

# -----------------------------------------------------------------------------
# 5. Read one file: crop to each site box, summarise every time step
# -----------------------------------------------------------------------------
nc_dates <- function(nc, time_name) {
  units <- ncatt_get(nc, time_name, "units")$value
  step <- c(seconds = 1, hours = 3600, days = 86400)[[word(units, 1)]]
  origin <- as.POSIXct(str_remove(units, "^\\w+ since "), tz = "UTC")
  as.Date(origin + as.vector(ncvar_get(nc, time_name)) * step)
}

read_fsle_file <- function(path, sites) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  dim_name <- function(options) intersect(options, names(nc$dim))[1]
  lon_name <- dim_name(c("lon", "longitude"))
  lat_name <- dim_name(c("lat", "latitude"))
  lon <- as.vector(ncvar_get(nc, lon_name))
  lat <- as.vector(ncvar_get(nc, lat_name))
  dates <- nc_dates(nc, "time")

  # Variables must be stored as (lon, lat, time) for the start/count below.
  for (v in c("fsle_max", "theta_max")) {
    dims <- map_chr(nc$var[[v]]$dim, "name")
    if (!identical(dims, c(lon_name, lat_name, "time"))) {
      stop(v, " dimensions are (", paste(dims, collapse = ", "), ") in ", basename(path))
    }
  }

  # AVISO longitudes are 0-360; config.R boxes are -180..180.
  to_file_lon <- function(x) if (max(lon) > 180) x %% 360 else x

  # Read each variable once (whole file) and crop it to every box.
  fsle_all  <- ncvar_get(nc, "fsle_max",  collapse_degen = FALSE)
  theta_all <- ncvar_get(nc, "theta_max", collapse_degen = FALSE)

  map_dfr(sites, function(s) {
    imap_dfr(BOXES, function(bx, box_name) {
      b <- filter(bx$corners, site == s)
      i_lon <- which(lon >= to_file_lon(b$lon_min) & lon <= to_file_lon(b$lon_max))
      i_lat <- which(lat >= b$lat_min & lat <= b$lat_max)
      if (length(i_lon) == 0 || length(i_lat) == 0) {
        stop("Site ", s, " ", box_name, " box is outside the grid of ", basename(path))
      }
      fsle  <- fsle_all[i_lon, i_lat, , drop = FALSE]
      theta <- theta_all[i_lon, i_lat, , drop = FALSE]
      map_dfr(seq_along(dates), function(k) {
        summarise_box(fsle[, , k], theta[, , k]) %>%
          mutate(theta_min = suppressWarnings(min(theta[, , k], na.rm = TRUE)),
                 theta_max = suppressWarnings(max(theta[, , k], na.rm = TRUE)))
      }) %>%
        mutate(box = box_name, site = s, date = dates, .before = 1)
    })
  })
}

fsle_raw <- pmap_dfr(jobs, function(path, sites, period) {
  out <- read_fsle_file(path, sites)
  # Keep only the days a box-cropped file is named for.
  if (is.na(period)) return(out)
  fmt <- c(`4` = "%Y", `6` = "%Y%m", `8` = "%Y%m%d")[[as.character(nchar(period))]]
  filter(out, format(date, fmt) == period)
}) %>%
  inner_join(clim_window, by = "site") %>%
  filter(date >= start, date <= end) %>%
  select(-start, -end)

# -----------------------------------------------------------------------------
# 6. Checks and coverage report
# -----------------------------------------------------------------------------
# The same day must not come from two files (e.g., cropped and global copies).
dups <- count(fsle_raw, box, site, date) %>% filter(n > 1)
if (nrow(dups) > 0) stop(nrow(dups), " site-days appear in more than one FSLE file")

# Every box must contain ocean cells; none would mean a wrong box/longitude.
if (any(fsle_raw$n_cells == 0)) stop("Site-days with no valid FSLE cells in the box")

# Grid cells per box (0.04 degree grid), for the log.
message("Valid FSLE grid cells per box (median over days):")
print(fsle_raw %>% group_by(box, site) %>% summarise(n_cells = median(n_cells), .groups = "drop") %>%
        pivot_wider(names_from = box, values_from = n_cells))

# Orientation range (degrees). The DT2024 files use -180..180;
# |theta| <= pi everywhere would suggest radians instead.
theta_range <- range(c(fsle_raw$theta_min, fsle_raw$theta_max))
message(sprintf("theta_max observed range: %.1f .. %.1f degrees", theta_range[1], theta_range[2]))
if (max(abs(theta_range)) <= pi + 1e-6) warning("theta_max looks like radians, not degrees")

# Sign convention of the magnitude: backward FSLE should be <= 0.
pct_pos <- 100 * mean(fsle_raw$FSLE > 0, na.rm = TRUE)
if (pct_pos > 0) {
  warning(sprintf("%.1f%% of daily box-mean FSLE values are positive; expected <= 0 ", pct_pos),
          "for backward-in-time FSLE - check the product's sign convention")
}

# Days without any cell reaching the FSLE threshold have NA orientation;
# report them, since they become NA predictors downstream.
message("Days with no box cell reaching the FSLE threshold (orientation NA):")
print(fsle_raw %>% group_by(box, site) %>%
        summarise(days_without_fsle_cells = sum(frac_cells_with_fsle == 0), .groups = "drop"))

# Missing days per site (identical for both boxes: same files).
requests <- clim_window %>%
  mutate(date = map2(start, end, seq, by = "day")) %>%
  select(site, date) %>%
  unnest(date) %>%
  left_join(fsle_raw %>% filter(box == "main") %>% select(site, date, FSLE),
            by = c("site", "date"))
message("FSLE days missing per site (skipped):")
print(requests %>% group_by(site) %>%
        summarise(days_requested = n(), days_missing = sum(is.na(FSLE)), .groups = "drop"))

# Missing recording days matter most: FSLE is a candidate model predictor.
missing_deploy <- requests %>%
  filter(is.na(FSLE)) %>%
  inner_join(select(SITES, site, rec_start, rec_end), by = "site") %>%
  filter(date >= as.Date(rec_start), date <= as.Date(rec_end))
if (nrow(missing_deploy) > 0) {
  warning(nrow(missing_deploy), " recording-period site-days have no FSLE data, e.g.: ",
          paste(head(paste(missing_deploy$site, missing_deploy$date), 20), collapse = ", "))
}

# -----------------------------------------------------------------------------
# 7. Save (main box: data-contract file; comparison box: "_box40" suffix)
# -----------------------------------------------------------------------------
iwalk(BOXES, function(bx, box_name) {
  out <- fsle_raw %>%
    filter(box == box_name) %>%
    select(date, site, FSLE, FSLE_sd, fsle_orient, fsle_orient_sd, frac_cells_with_fsle) %>%
    arrange(factor(site, levels = SITE_ORDER), date)
  write_csv(out, with_suffix(F_FSLE_DAILY, bx$suffix))
  message("Wrote ", nrow(out), " site-days to ", with_suffix(F_FSLE_DAILY, bx$suffix))
})
