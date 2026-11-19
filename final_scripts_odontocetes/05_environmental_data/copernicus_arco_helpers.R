# =============================================================================
# Script  : copernicus_arco_helpers.R
# Section : 05_environmental_data (shared helper functions)
# Purpose : Functions for downloading subsets of Copernicus Marine products
#           from the ARCO (zarr) service with the R package CopernicusMarine
#           0.3.2 and converting them to correctly scaled, correctly located
#           arrays. Sourced by scripts that download Copernicus data outside
#           the main site-box pipeline (e.g., 01_study_area/currents_map.R).
# Inputs  : none (functions only). Callers need a Copernicus Marine account
#           (credentials as environment variables
#           COPERNICUSMARINE_SERVICE_USERNAME / COPERNICUSMARINE_SERVICE_PASSWORD)
#           and the packages CopernicusMarine (0.3.2), stars, httr2, blosc.
# Outputs : none
# Author  : Natalie Posdaljian
# Date    : 2026-11-17
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Same workarounds as 05_environmental_data/01_download_copernicus.R
#           (validated there against the original portal downloads), written
#           as stand-alone functions:
#             (1) PHY values come back as raw packed int16 -> scale_factor and
#                 add_offset are read from the zarr metadata and applied;
#             (2) longitude is reported one grid cell east of the true cell
#                 centres -> true centre = reported cell start - delta/2;
#             (3) variables whose unit string the 'units' package cannot parse
#                 (salinity, "1e-3") are silently dropped by 0.3.2 -> the unit
#                 is replaced by "1" before processing.
#           A fourth issue was found with the monthly PHY layer over a
#           regional (145 x 55 cell) subset: a request spanning several time
#           chunks came back with cells scrambled across time and latitude,
#           while single-time requests were correct. Callers should request
#           one time step per call unless they verify the result (e.g., a
#           land mask that is constant over time).
#           [CONFIRM] 01_download_copernicus.R keeps its own copies of these
#           functions (it was running when this file was written); it can be
#           switched to source this file once its downloads are complete.
# =============================================================================

library(CopernicusMarine)
library(stars)
library(httr2)

# Stop unless the exact package version the workarounds were validated with
# is installed, and the zarr decompression dependency is available.
check_copernicus_setup <- function() {
  if (packageVersion("CopernicusMarine") != "0.3.2") {
    stop("Validated with CopernicusMarine 0.3.2; found ",
         packageVersion("CopernicusMarine"), ". Re-validate before use.")
  }
  if (!requireNamespace("blosc", quietly = TRUE)) {
    stop("Package 'blosc' is required to decompress ARCO (zarr) chunks.")
  }
  if (cms_get_username() == "" || cms_get_password() == "") {
    stop("Copernicus Marine credentials not found. Set the environment ",
         "variables COPERNICUSMARINE_SERVICE_USERNAME and ",
         "COPERNICUSMARINE_SERVICE_PASSWORD (e.g., in ~/.Renviron).")
  }
  invisible(TRUE)
}

# Catalogue entry (metadata row) of a layer within a product. Layer ids in
# the catalogue carry a version suffix (e.g., "_202311"), so prefix matching
# is used.
layer_metadata <- function(product, layer) {
  meta <- cms_product_metadata(product)
  i <- which(startsWith(meta$id, layer))[1]
  if (is.na(i)) stop("Layer ", layer, " not found in product ", product)
  list(properties = meta$properties[[i]],
       href = meta$assets[[i]]$timeChunked$href)
}

# Temporal extent (first and last time stamp) served by a layer.
layer_time_extent <- function(product, layer) {
  props <- layer_metadata(product, layer)$properties
  as.Date(c(props$start_datetime, props$end_datetime))
}

# Packing and fill attributes of each variable, read from the public zarr
# metadata files (.zattrs, .zarray) of the layer.
read_zarr_json <- function(url) {
  resp_body_json(req_perform(req_retry(request(url), max_tries = 3)),
                 check_type = FALSE)
}

layer_zarr_attrs <- function(product, layer, vars) {
  href <- layer_metadata(product, layer)$href
  setNames(lapply(vars, function(v) {
    list(zattrs = read_zarr_json(paste0(href, "/", v, "/.zattrs")),
         zarray = read_zarr_json(paste0(href, "/", v, "/.zarray")))
  }), vars)
}

# Same steps as CopernicusMarine::cms_download_subset() 0.3.2, except that
# unparseable variable units are replaced by "1" so the variable is not
# dropped (workaround 3). Returns a stars object with raw stored values.
download_subset <- function(product, layer, variable, region, timerange,
                            verticalrange = NULL) {
  ns <- asNamespace("CopernicusMarine")
  subset_request <- list(product = product, layer = layer, variable = variable,
                         region = region, timerange = timerange,
                         verticalrange = verticalrange)
  service <- ns$.get_best_arco_service_type(subset_request, "", FALSE,
                                            "default")
  refsys  <- ns$.get_refsys(attributes(service)$dim_properties)

  var_props <- attr(service, "var_properties")
  for (v in names(var_props)) {
    ok <- tryCatch({
      units::as_units(1, var_props[[v]]$unit)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) var_props[[v]]$unit <- "1"
  }
  attr(service, "var_properties") <- var_props

  vars <- names(service$viewVariables)
  dims <- names(attributes(service)$dims)
  serv_props <- ns$.get_service_properties(service, vars)
  if (any(serv_props$zarr_format != 2)) {
    stop("Only zarr format version 2 is supported")
  }
  chunks <- ns$.download_chunks(dims, service, vars, FALSE)
  chunks <- ns$.decompress_data(chunks, serv_props)
  x <- ns$.process_data(chunks, dims, vars, service, serv_props, TRUE, refsys)
  if (!all(variable %in% names(x))) {
    stop("Variable(s) missing from download: ", variable)
  }
  x
}

# Convert raw stored values to physical units (workaround 1) and check them
# against a plausible range (NULL = no range check).
unpack_values <- function(raw, var, attrs, plausible_range = NULL) {
  za <- attrs$zattrs
  fill <- attrs$zarray$fill_value
  if (is.numeric(fill)) raw[!is.na(raw) & raw == fill] <- NA
  if (!is.null(za$scale_factor)) {
    # Values must still be packed integers within the stored valid range;
    # anything else would mean the package already unpacked them.
    rng <- range(raw, na.rm = TRUE)
    if (any(raw != round(raw), na.rm = TRUE) ||
          rng[1] < za$valid_min || rng[2] > za$valid_max) {
      stop(var, ": values are not packed integers within valid_min/valid_max")
    }
    raw <- raw * za$scale_factor + (za$add_offset %||% 0)
  }
  if (!is.null(plausible_range)) {
    rng <- range(raw, na.rm = TRUE)
    if (rng[1] < plausible_range[1] || rng[2] > plausible_range[2]) {
      stop(sprintf("%s: values %.4g..%.4g outside plausible range %g..%g",
                   var, rng[1], rng[2], plausible_range[1],
                   plausible_range[2]))
    }
  }
  raw
}

# Physical array and coordinates of one variable of a downloaded stars
# object, with dimensions reordered to (longitude, latitude, [depth,] time).
# Horizontal coordinates are corrected to true cell centres (workaround 2)
# and must lie on the product grid (`grid_res`, degrees), otherwise stop.
# Depth is returned as positive metres, shallowest first.
subset_to_array <- function(x, var, attrs, grid_res, plausible_range = NULL) {
  dims <- st_dimensions(x)
  find_dim <- function(pattern) {
    hit <- grep(pattern, names(dims), ignore.case = TRUE, value = TRUE)
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  d_lon  <- find_dim("^(x|lon|longitude)$")
  d_lat  <- find_dim("^(y|lat|latitude)$")
  d_dep  <- find_dim("^(elevation|depth|z)$")
  d_time <- find_dim("^time$")
  if (anyNA(c(d_lon, d_lat, d_time))) {
    stop("Unexpected dimensions in downloaded object: ",
         paste(names(dims), collapse = ", "))
  }

  on_grid <- function(v) all(abs(v / grid_res - round(v / grid_res)) < 1e-3)
  get_coord <- function(d) {
    start <- st_get_dimension_values(x, d, center = FALSE)
    coord <- if (isTRUE(dims[[d]]$point)) start else start - dims[[d]]$delta / 2
    if (!on_grid(coord)) {
      stop("Coordinates of dimension '", d, "' are not on the ", grid_res,
           " degree grid")
    }
    coord
  }

  keep <- c(d_lon, d_lat, if (!is.na(d_dep)) d_dep, d_time)
  arr  <- array(as.numeric(x[[var]]), dim = dim(x[[var]]))
  arr  <- aperm(unpack_values(arr, var, attrs, plausible_range),
                match(keep, names(dims)))

  depth <- NULL
  if (!is.na(d_dep)) {
    depth <- abs(as.numeric(st_get_dimension_values(x, d_dep)))
    ord   <- order(depth)
    depth <- depth[ord]
    arr   <- arr[, , ord, , drop = FALSE]
  }

  list(values = arr,
       lon    = get_coord(d_lon),
       lat    = get_coord(d_lat),
       depth  = depth,
       time   = as.POSIXct(st_get_dimension_values(x, d_time, center = FALSE),
                           tz = "UTC"))
}
