# =============================================================================
# Script  : 01_download_copernicus.R
# Section : 05_environmental_data
# Purpose : Download daily Copernicus Marine model output (GLORYS12 physics
#           reanalysis and PISCES biogeochemistry hindcast) for the +/- 40 km
#           box around each HARP site (SITES$box40_*, which contains the main
#           +/- 20 km analysis box), from the surface to ~1700 m, over the
#           climatology window (deployment start year - CLIM_HALF_WINDOW_YEARS
#           to deployment end year + CLIM_HALF_WINDOW_YEARS).
# Inputs  : final_scripts_odontocetes/config.R (SITES boxes and recording
#           times, COPERNICUS products/layers/variables, DOWNLOAD_DIR,
#           CLIM_HALF_WINDOW_YEARS). Requires a free Copernicus Marine account
#           and an internet connection.
# Outputs : DOWNLOAD_DIR/copernicus/{site}_{phy|bgc}_{variable}_{year}.nc
#           One NetCDF per site x product x variable x calendar year, with
#           dimensions longitude, latitude, [depth,] time (seconds since
#           1970-01-01 UTC), i.e., the same layout as files downloaded from the
#           Copernicus Marine web portal. Read by 04_extract_copernicus.R.
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - The deployment-period files used in the original analysis
#             ("Environmental Data/Copernicus/Physics Reanalysis/{site}40km_*"
#             and ".../Biogeochemistry Hindcast/{site}40km_*") were downloaded
#             manually through the Copernicus Marine web portal
#             (https://data.marine.copernicus.eu/) using the same +/- 40 km box
#             corners and variables. This script reproduces those downloads
#             programmatically and extends them to the climatology window
#             needed for deseasoning and for the lagged chl-a/productivity
#             predictors.
#           - Credentials are NEVER written in this script. Set them once per
#             machine as environment variables (e.g., in ~/.Renviron):
#                 COPERNICUSMARINE_SERVICE_USERNAME=<your user name>
#                 COPERNICUSMARINE_SERVICE_PASSWORD=<your password>
#             or per session with CopernicusMarine::cms_set_username() /
#             cms_set_password(). They are read with cms_get_username() /
#             cms_get_password().
#           - Uses the R package CopernicusMarine 0.3.2 (checked at start; the
#             workarounds below are specific to that version). Its subset
#             download returns a stars object in memory (`destination` is
#             deprecated), so this script writes each subset to NetCDF itself.
#             Reading the ARCO (zarr) service also requires 'blosc'.
#             Three issues of 0.3.2 are handled here:
#             (1) PHY values come back as raw packed int16; scale_factor and
#                 add_offset (read from the zarr metadata) are applied here;
#             (2) longitude coordinates are reported one grid cell east of the
#                 true cell centres; corrected in write_subset_nc();
#             (3) a variable whose unit string cannot be parsed by the 'units'
#                 package (salinity, unit "1e-3") is silently DROPPED from the
#                 result; download_subset() replicates cms_download_subset()
#                 with unparseable units replaced by "1" (the NetCDF units
#                 attribute is taken from the zarr metadata instead).
#             Validated 2026-11-17: EI 2014-03-05..31 downloads match the
#             original portal files cell-by-cell and day-by-day for every
#             variable (see report / validation script in the project log).
#           - Downloads are split by variable and calendar year to keep memory
#             use small and to make the script resumable: files that already
#             exist are skipped, so the script can simply be re-run after an
#             interrupted session.
#           - Requested time ranges are truncated to each layer's temporal
#             extent, which is queried from the Copernicus metadata catalogue
#             at run time.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

library(tidyverse)
library(CopernicusMarine)
library(stars)
library(ncdf4)
library(httr2)

# -----------------------------------------------------------------------------
# 1. Script parameters (not in config.R because only this script uses them)
# -----------------------------------------------------------------------------
# Deepest level requested (m). 1700 m keeps the deepest level used for sperm
# whales (nominal 1665 m: PHY level 1684.3 m, BGC level 1652.6 m) and excludes
# the next model levels (PHY 1941.9 m, BGC 1795.7 m), matching the original
# portal downloads.
DOWNLOAD_MAX_DEPTH_M <- 1700

# Two-dimensional (surface) variables are requested without a vertical range.
SURFACE_VARS <- c("zos", "mlotst", "siconc", "sithick")

# Horizontal grid spacing (degrees) of each layer; used to check that the
# coordinates written to NetCDF are grid-cell centres.
GRID_RES_DEG <- c(phy = 1 / 12, bgc = 0.25)

COP_DOWNLOAD_DIR <- file.path(DOWNLOAD_DIR, "copernicus")
dir.create(COP_DOWNLOAD_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Check credentials and dependencies before starting
# -----------------------------------------------------------------------------
# Fail early with a clear message instead of part-way through the downloads.
if (cms_get_username() == "" || cms_get_password() == "") {
  stop("Copernicus Marine credentials not found. Set the environment variables ",
       "COPERNICUSMARINE_SERVICE_USERNAME and COPERNICUSMARINE_SERVICE_PASSWORD ",
       "(e.g., in ~/.Renviron) or call cms_set_username()/cms_set_password().")
}
if (!requireNamespace("blosc", quietly = TRUE)) {
  stop("Package 'blosc' is required by CopernicusMarine to decompress ARCO ",
       "(zarr) chunks: install.packages('blosc')")
}
# The workarounds in download_subset() and write_subset_nc() rely on the
# internals and behaviour of this exact version; re-validate before changing.
if (packageVersion("CopernicusMarine") != "0.3.2") {
  stop("This script was validated with CopernicusMarine 0.3.2; found ",
       packageVersion("CopernicusMarine"), ". Re-run the validation before use.")
}

# Same steps as CopernicusMarine::cms_download_subset() 0.3.2, with one change:
# variable units that the 'units' package cannot parse (salinity: "1e-3") are
# replaced by "1" before the final formatting step, because 0.3.2 otherwise
# drops such variables from the result without an error.
download_subset <- function(product, layer, variable, region, timerange,
                            verticalrange = NULL) {
  ns <- asNamespace("CopernicusMarine")
  subset_request <- list(product = product, layer = layer, variable = variable,
                         region = region, timerange = timerange,
                         verticalrange = verticalrange)
  service <- ns$.get_best_arco_service_type(subset_request, "", FALSE, "default")
  refsys  <- ns$.get_refsys(attributes(service)$dim_properties)

  var_props <- attr(service, "var_properties")
  for (v in names(var_props)) {
    ok <- tryCatch({units::as_units(1, var_props[[v]]$unit); TRUE},
                   error = function(e) FALSE)
    if (!ok) var_props[[v]]$unit <- "1"
  }
  attr(service, "var_properties") <- var_props

  vars <- names(service$viewVariables)
  dims <- names(attributes(service)$dims)
  serv_props <- ns$.get_service_properties(service, vars)
  if (any(serv_props$zarr_format != 2)) stop("Only zarr format version 2 is supported")
  chunks <- ns$.download_chunks(dims, service, vars, FALSE)
  chunks <- ns$.decompress_data(chunks, serv_props)
  x <- ns$.process_data(chunks, dims, vars, service, serv_props, TRUE, refsys)
  if (!all(variable %in% names(x))) stop("Variable(s) missing from download: ", variable)
  x
}

# -----------------------------------------------------------------------------
# 3. Climatology window per site, truncated to layer availability
# -----------------------------------------------------------------------------
# The window spans whole calendar years so that every month is equally
# represented in the monthly climatology used for deseasoning.
clim_window <- SITES %>%
  transmute(site,
            start = as.Date(sprintf("%d-01-01",
                                    year(rec_start) - CLIM_HALF_WINDOW_YEARS)),
            end   = as.Date(sprintf("%d-12-31",
                                    year(rec_end) + CLIM_HALF_WINDOW_YEARS)))

# Temporal extent of each daily layer, read from the public metadata catalogue
# (no credentials needed). Copernicus extends the multiyear products over time,
# so querying avoids hard-coding an end date that will go stale.
# [CONFIRM] that the full window is served by the multiyear ("my") layer and
# not partly by the interim ("myint") extension; check the extent printed here
# against the product user manual when the download is run.
layer_extent <- function(product, layer) {
  meta <- cms_product_metadata(product)
  i <- which(str_starts(meta$id, fixed(layer)))[1]
  if (is.na(i)) stop("Layer ", layer, " not found in product ", product)
  props <- meta$properties[[i]]
  as.Date(c(props$start_datetime, props$end_datetime))
}
extents <- map(COPERNICUS, ~ layer_extent(.x$product, .x$layer))
print(extents)

# -----------------------------------------------------------------------------
# 4. Variable packing and plausibility
# -----------------------------------------------------------------------------
# The ARCO (zarr) store keeps the PHY variables as packed 16-bit integers
# (physical = raw * scale_factor + add_offset); BGC variables are float32.
# cms_download_subset() 0.3.2 returns the RAW stored values without applying
# scale_factor/add_offset (it only converts the zarr fill value to NA), so the
# packing attributes are read here from the store's public .zattrs/.zarray
# and applied before writing.
read_zarr_json <- function(url) {
  resp_body_json(req_perform(req_retry(request(url), max_tries = 3)), check_type = FALSE)
}
layer_attrs <- function(product, layer, vars) {
  meta <- cms_product_metadata(product)
  i <- which(str_starts(meta$id, fixed(layer)))[1]
  href <- meta$assets[[i]]$timeChunked$href
  set_names(map(vars, function(v) {
    list(zattrs = read_zarr_json(paste0(href, "/", v, "/.zattrs")),
         zarray = read_zarr_json(paste0(href, "/", v, "/.zarray")))
  }), vars)
}
zarr_attrs <- map(COPERNICUS, ~ layer_attrs(.x$product, .x$layer, .x$vars))

# Plausible physical ranges, used to catch packed or otherwise wrong values
# before anything is written (same limits as in 04_extract_copernicus.R).
PLAUSIBLE_RANGE <- list(thetao = c(-3, 15), so = c(30, 36), uo = c(-3, 3), vo = c(-3, 3),
                        zos = c(-3, 3), mlotst = c(0, 2000), siconc = c(0, 1),
                        sithick = c(0, 10), chl = c(0, 50), o2 = c(0, 500),
                        nppv = c(-1, 500))

unpack <- function(raw, var, attrs) {
  za <- attrs$zattrs
  # A remaining float fill value (e.g., 9.97e36 for BGC) is set to NA.
  fill <- attrs$zarray$fill_value
  if (is.numeric(fill)) raw[!is.na(raw) & raw == fill] <- NA
  if (!is.null(za$scale_factor)) {
    # Raw values must still be packed int16 integers; anything else would mean
    # the package already unpacked them. (The stored valid_min/valid_max are
    # not used: real zos values fall slightly outside them, e.g. raw -6252 <
    # valid_min -6228 at EI in 2010. The physical-range check below applies.)
    rng <- range(raw, na.rm = TRUE)
    if (any(raw != round(raw), na.rm = TRUE) || rng[1] < -32768 || rng[2] > 32767) {
      stop(var, ": downloaded values are not packed 16-bit integers")
    }
    raw <- raw * za$scale_factor + (za$add_offset %||% 0)
  }
  lim <- PLAUSIBLE_RANGE[[var]]
  rng <- range(raw, na.rm = TRUE)
  # Allow a small tolerance (0.1% of the range width) so that int16 rounding
  # (e.g., siconc unpacked to 1.00002) does not trigger a false alarm, while
  # still catching packed or mis-scaled values, which are orders of magnitude off.
  tol <- 1e-3 * diff(lim)
  if (rng[1] < lim[1] - tol || rng[2] > lim[2] + tol) {
    stop(sprintf("%s: values %.4g..%.4g outside plausible range %g..%g",
                 var, rng[1], rng[2], lim[1], lim[2]))
  }
  raw
}

# -----------------------------------------------------------------------------
# 5. Write one downloaded stars object to NetCDF
# -----------------------------------------------------------------------------
# The stars object returned by cms_download_subset() has dimensions named after
# the zarr store (longitude/latitude or x/y, elevation, time) in an order that
# depends on the service. This helper writes a NetCDF with the fixed layout
# (longitude, latitude, [depth,] time) used by the portal files, so that
# 04_extract_copernicus.R can read old and new files with the same code.
# Depth is written as positive metres below the surface.
write_subset_nc <- function(x, var, attrs, file, grid_res) {
  dims    <- st_dimensions(x)
  find_dim <- function(pattern) {
    hit <- grep(pattern, names(dims), ignore.case = TRUE, value = TRUE)
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  d_lon   <- find_dim("^(x|lon|longitude)$")
  d_lat   <- find_dim("^(y|lat|latitude)$")
  d_dep   <- find_dim("^(elevation|depth|z)$")
  d_time  <- find_dim("^time$")
  if (anyNA(c(d_lon, d_lat, d_time))) {
    stop("Unexpected dimensions in downloaded object: ",
         paste(names(dims), collapse = ", "))
  }

  # Horizontal coordinates. In CopernicusMarine 0.3.2 the two horizontal
  # dimensions come back differently: latitude is a "point" dimension whose
  # values are the true cell centres, whereas longitude is a regular "cell"
  # dimension whose offset is one full cell EAST of the true first centre
  # (offset = centre + delta/2 instead of centre - delta/2). Verified
  # cell-by-cell against the portal files (see validation note in the header):
  # true longitude centre = reported cell start - delta/2.
  # The result must lie on the product grid, otherwise stop.
  on_grid <- function(v) all(abs(v / grid_res - round(v / grid_res)) < 1e-3)
  get_coord <- function(d) {
    start <- st_get_dimension_values(x, d, center = FALSE)
    is_point <- isTRUE(dims[[d]]$point)
    coord <- if (is_point) start else start - dims[[d]]$delta / 2
    if (!on_grid(coord)) {
      stop("Coordinates of dimension '", d, "' are not on the ", grid_res, " degree grid")
    }
    coord
  }
  lon  <- get_coord(d_lon)
  lat  <- get_coord(d_lat)
  # Daily means are stamped at 00:00 UTC, as in the portal files.
  time <- as.POSIXct(st_get_dimension_values(x, d_time, center = FALSE), tz = "UTC")

  # Unpack to physical units (see unpack()), then reorder the array to
  # (lon, lat, [depth,] time).
  keep <- c(d_lon, d_lat, if (!is.na(d_dep)) d_dep, d_time)
  arr  <- array(as.numeric(x[[var]]), dim = dim(x[[var]]))
  arr  <- aperm(unpack(arr, var, attrs), match(keep, names(dims)))

  nc_dims <- list(ncdim_def("longitude", "degrees_east", lon),
                  ncdim_def("latitude", "degrees_north", lat))
  if (!is.na(d_dep) && var %in% SURFACE_VARS) {
    # Surface (2-D) fields come back with a degenerate one-level vertical
    # dimension; drop it so the file matches the portal layout
    # (longitude, latitude, time).
    if (dim(arr)[3] != 1) stop(var, ": expected a single vertical level")
    arr <- array(arr, dim = dim(arr)[-3])
  } else if (!is.na(d_dep)) {
    # Depth as positive metres, shallowest level first (as in the portal files).
    depth <- abs(as.numeric(st_get_dimension_values(x, d_dep)))
    ord   <- order(depth)
    depth <- depth[ord]
    arr   <- arr[, , ord, , drop = FALSE]
    nc_dims <- c(nc_dims, list(ncdim_def("depth", "m", depth)))
  }
  var_units <- attrs$zattrs$units %||% ""
  nc_dims <- c(nc_dims, list(ncdim_def("time", "seconds since 1970-01-01 00:00:00",
                                       as.numeric(time), unlim = TRUE)))
  nc_var <- ncvar_def(var, units = var_units, dim = nc_dims, missval = NaN,
                      prec = "float")

  # Write to a temporary name first so that an interrupted write never leaves
  # a truncated file that would be skipped on the next run.
  tmp <- paste0(file, ".part")
  nc  <- nc_create(tmp, nc_var)
  ncvar_put(nc, nc_var, arr)
  nc_close(nc)
  file.rename(tmp, file)
}

# -----------------------------------------------------------------------------
# 6. Download loop: site x product x variable x year
# -----------------------------------------------------------------------------
for (s in SITE_ORDER) {
  si <- site_info(s)
  # Download region = the +/- 40 km box (box40_*), which contains the main
  # +/- 20 km analysis box, so one download serves both boxes. Order required
  # by the subset service: xmin, ymin, xmax, ymax.
  region <- c(si$box40_lon_min, si$box40_lat_min, si$box40_lon_max, si$box40_lat_max)
  win <- filter(clim_window, site == s)

  for (p in names(COPERNICUS)) {
    prod <- COPERNICUS[[p]]
    start <- max(win$start, extents[[p]][1])
    end   <- min(win$end, extents[[p]][2])
    if (start > win$start || end < win$end) {
      message(sprintf("%s %s: window truncated to %s .. %s by layer availability",
                      s, p, start, end))
    }

    for (v in prod$vars) {
      for (yr in year(start):year(end)) {
        out_file <- file.path(COP_DOWNLOAD_DIR, sprintf("%s_%s_%s_%d.nc", s, p, v, yr))
        if (file.exists(out_file)) next   # resumable: skip completed years

        t0 <- max(start, as.Date(sprintf("%d-01-01", yr)))
        t1 <- min(end, as.Date(sprintf("%d-12-31", yr)))
        message(sprintf("Downloading %s %s %s %s .. %s", s, p, v, t0, t1))

        args <- list(product   = prod$product,
                     layer     = prod$layer,
                     variable  = v,
                     region    = region,
                     timerange = as.POSIXct(format(c(t0, t1)), tz = "UTC"))
        # Depth is requested only for 3-D variables; values below the sea
        # surface are negative for the subset service.
        if (!v %in% SURFACE_VARS) args$verticalrange <- c(0, -DOWNLOAD_MAX_DEPTH_M)

        x <- do.call(download_subset, args)
        write_subset_nc(x, v, zarr_attrs[[p]][[v]], out_file, GRID_RES_DEG[[p]])
        rm(x)
        gc(verbose = FALSE)
      }
    }
  }
}

message("Copernicus downloads complete: ", COP_DOWNLOAD_DIR)
