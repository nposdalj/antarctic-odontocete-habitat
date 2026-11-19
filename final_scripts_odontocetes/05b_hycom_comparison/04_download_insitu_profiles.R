# =============================================================================
# Script  : 04_download_insitu_profiles.R
# Section : Supplementary (2.5.1) - model validation against in-situ profiles
# Purpose : Download all in-situ temperature, salinity and dissolved-oxygen
#           profiles from the World Ocean Database (WOD) in the wider study
#           region for 2014-2016, so that the Copernicus and HYCOM models can
#           be compared with measurements at depth (05_*.R and 06_*.R).
# Inputs  : config.R (DOWNLOAD_DIR)
#           WOD yearly files per instrument type on the NCEI THREDDS OPeNDAP
#           server, https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/wod/
# Outputs : DOWNLOAD_DIR/wod/wod_<probe>_<year>_region.csv (one per file with
#             casts in the region; cache, so the script can be re-run)
#           outputs/05b_hycom_comparison/insitu_profiles.csv
#             cast_id, probe, platform, time (UTC), lat, lon, depth_m,
#             temperature (deg C), temperature_flag, temperature_profile_flag,
#             salinity, salinity_flag, salinity_profile_flag,
#             oxygen (umol/kg), oxygen_flag, oxygen_profile_flag
# Author  : Natalie Posdaljian
# Date    : 2026-11-19
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Region: INSITU_LAT / INSITU_LON (60-63 S, 52-60 W), which
#             contains all three sites and Bransfield Strait; the same region
#             as the 2026-11-16 availability check (README).
#           - Instrument types (WOD "probe" codes): apb = animal-borne
#             (instrumented seals), ctd, drb = drifting buoy, gld = glider,
#             mrb = moored buoy, osd = bottle / low-resolution CTD, pfl =
#             profiling float (Argo), xbt = expendable bathythermograph
#             (temperature only).
#           - WOD files are ragged arrays: each variable has its own
#             observation dimension and per-cast row size; a variable's
#             values for a cast correspond to that cast's depths (z) when its
#             row size equals the cast's z row size. Casts where this does not
#             hold are skipped for that variable (counted in the log).
#           - Quality flags are kept, not applied here: WOD observation flag
#             0 = accepted, profile flag 0 = accepted. 06_*.R uses only values
#             with both flags = 0.
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
INSITU_LAT   <- c(-63, -60)
INSITU_LON   <- c(-60, -52)
INSITU_YEARS <- 2014:2016
WOD_PROBES   <- c("apb", "ctd", "drb", "gld", "mrb", "osd", "pfl", "xbt")
WOD_OPENDAP  <- "https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/wod"
WOD_VARS     <- c(temperature = "Temperature", salinity = "Salinity",
                  oxygen = "Oxygen")
WOD_DIR      <- file.path(DOWNLOAD_DIR, "wod")
SECTION      <- "05b_hycom_comparison"
dir.create(WOD_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
# Row sizes with missing values (no data for that cast) as 0.
row_sizes <- function(nc, var) {
  if (!var %in% names(nc$var)) return(NULL)
  r <- as.vector(ncvar_get(nc, var))
  r[is.na(r)] <- 0
  r
}

# Values [start, start + count) of a 1-D variable (1-based start).
read_slice <- function(nc, var, start, count) {
  as.vector(ncvar_get(nc, var, start = start, count = count))
}

# Optional per-cast character variable (e.g. platform name)
cast_text <- function(nc, var, i) {
  if (!var %in% names(nc$var)) return(NA_character_)
  out <- tryCatch(ncvar_get(nc, var, start = c(1, i), count = c(-1, 1)),
                  error = function(e) NA_character_)
  trimws(as.character(out))
}

# All casts of one WOD file inside the region, as a long table (one row per
# cast x depth).
read_wod_file <- function(probe, year) {
  url <- sprintf("%s/%d/wod_%s_%d.nc", WOD_OPENDAP, year, probe, year)
  nc <- tryCatch(nc_open(url), error = function(e) NULL)
  if (is.null(nc)) {
    message(sprintf("  %s %d: file not available", probe, year))
    return(NULL)
  }
  on.exit(nc_close(nc))

  lat <- as.vector(ncvar_get(nc, "lat"))
  lon <- as.vector(ncvar_get(nc, "lon"))
  in_region <- which(lat >= INSITU_LAT[1] & lat <= INSITU_LAT[2] &
                       lon >= INSITU_LON[1] & lon <= INSITU_LON[2])
  message(sprintf("  %s %d: %d of %d casts in region", probe, year,
                  length(in_region), length(lat)))
  if (length(in_region) == 0) return(tibble())

  time_units <- ncatt_get(nc, "time", "units")$value
  stopifnot(grepl("days since 1770-01-01", time_units))
  days  <- as.vector(ncvar_get(nc, "time"))
  cast  <- as.vector(ncvar_get(nc, "wod_unique_cast"))

  z_size  <- row_sizes(nc, "z_row_size")
  z_start <- cumsum(c(0, head(z_size, -1))) + 1
  var_info <- map(WOD_VARS, function(v) {
    size <- row_sizes(nc, paste0(v, "_row_size"))
    if (is.null(size)) return(NULL)
    list(size = size, start = cumsum(c(0, head(size, -1))) + 1,
         prof_flag = if (paste0(v, "_WODprofileflag") %in% names(nc$var))
           as.vector(ncvar_get(nc, paste0(v, "_WODprofileflag"))) else NULL)
  })

  skipped <- 0
  casts <- map(in_region, function(i) {
    if (z_size[i] == 0) return(NULL)
    out <- tibble(
      cast_id  = cast[i], probe = probe,
      platform = cast_text(nc, "Platform", i),
      time     = as.POSIXct("1770-01-01", tz = "UTC") + days[i] * 86400,
      lat = lat[i], lon = lon[i],
      depth_m  = read_slice(nc, "z", z_start[i], z_size[i])
    )
    for (name in names(WOD_VARS)) {
      v <- WOD_VARS[[name]]
      info <- var_info[[name]]
      value <- flag <- rep(NA_real_, z_size[i])
      prof_flag <- NA_real_
      if (!is.null(info) && info$size[i] > 0) {
        if (info$size[i] == z_size[i]) {
          value <- read_slice(nc, v, info$start[i], info$size[i])
          flag  <- read_slice(nc, paste0(v, "_WODflag"), info$start[i],
                              info$size[i])
          if (!is.null(info$prof_flag)) prof_flag <- info$prof_flag[i]
        } else {
          skipped <<- skipped + 1
        }
      }
      out[[name]] <- value
      out[[paste0(name, "_flag")]] <- flag
      out[[paste0(name, "_profile_flag")]] <- prof_flag
    }
    out
  }) |>
    list_rbind()
  if (skipped > 0) {
    message(sprintf("  %s %d: %d variable-casts skipped (row size != z)",
                    probe, year, skipped))
  }
  casts
}

# -----------------------------------------------------------------------------
# 3. Download (cached per WOD file)
# -----------------------------------------------------------------------------
for (year in INSITU_YEARS) {
  for (probe in WOD_PROBES) {
    cache <- file.path(WOD_DIR, sprintf("wod_%s_%d_region.csv", probe, year))
    if (file.exists(cache)) next
    message(sprintf("Reading wod_%s_%d.nc", probe, year))
    casts <- read_wod_file(probe, year)
    if (is.null(casts)) next          # file missing on the server: retry later
    write_csv(casts, cache)
  }
}

# -----------------------------------------------------------------------------
# 4. Combine
# -----------------------------------------------------------------------------
profiles <- list.files(WOD_DIR, pattern = "_region\\.csv$", full.names = TRUE) |>
  map(\(f) read_csv(f, show_col_types = FALSE,
                    col_types = cols(platform = col_character(),
                                     .default = col_guess()))) |>
  list_rbind() |>
  arrange(time, cast_id, depth_m)

write_csv(profiles, out_path(SECTION, "insitu_profiles.csv"))

profiles |>
  distinct(cast_id, probe, time) |>
  count(probe, year = lubridate::year(time)) |>
  pivot_wider(names_from = year, values_from = n, values_fill = 0) |>
  print()
