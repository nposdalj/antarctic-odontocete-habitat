# =============================================================================
# Script  : 04_extract_copernicus.R
# Section : 05_environmental_data
# Purpose : Summarise Copernicus physics (GLORYS12) and biogeochemistry
#           (PISCES) model output over two boxes around each HARP site - the
#           MAIN +/- 20 km analysis box and the +/- 40 km comparison box:
#           daily box mean and spatial variability at every model level, the
#           fraction of the box that contains water at each level, and daily
#           surface fields (sea surface height, mixed layer, sea ice).
# Inputs  : final_scripts_odontocetes/config.R (SITES boxes: box_* = main
#           +/- 20 km, box40_* = +/- 40 km; COPERNICUS variables,
#           DOWNLOAD_DIR, ENV_DIR)
#           Preferred: DOWNLOAD_DIR/copernicus/{site}_{phy|bgc}_{var}_{year}.nc
#             (climatology window, +/- 40 km download region, from
#             01_download_copernicus.R)
#           Fallback (deployment periods only, original portal downloads):
#             Environmental Data/Copernicus/Physics Reanalysis/
#               {site}40km_cmems_mod_glo_phy_my_0.083deg_P1D-m_*.nc
#             Environmental Data/Copernicus/Biogeochemistry Hindcast/
#               {site}40km_cmems_mod_glo_bgc_my_0.25deg_P1D-m_*.nc
# Outputs : Main box (data-contract names) and the same files with suffix
#           "_box40" for the +/- 40 km comparison box:
#           F_COP_PHY_DEPTH  copernicus_phy_by_depth[_box40].csv
#             date, site, depth_m, temperature, temperature_sd, salinity,
#             salinity_sd, u, v, u_mad, v_mad, coverage
#           F_COP_BGC_DEPTH  copernicus_bgc_by_depth[_box40].csv
#             date, site, depth_m, chla, chla_sd, o2, o2_sd, productivity,
#             productivity_sd, coverage
#           F_COP_SURFACE    copernicus_surface[_box40].csv
#             date, site, ssh, ssh_sd, mixed_layer, mixed_layer_sd, ice_conc,
#             ice_conc_sd, ice_thickness
#           F_COP_COVERAGE   copernicus_depth_coverage[_box40].csv
#             site, product, depth_m, n_cells, coverage
# Author  : Natalie Posdaljian
# Date    : 2026-11-18
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Boxes: a grid cell belongs to a box if its centre lies inside the
#             box corners. The BGC grid is coarse (0.25 degree), so the
#             +/- 20 km box holds only a few BGC cells (printed in the log and
#             in the coverage table, n_cells).
#           - Spatial variability is the standard deviation across box cells,
#             except for velocities (u, v), for which the median absolute
#             deviation (MAD, stats::mad, scaled to be consistent with the SD
#             for normal data) is used, as in the original analysis; u_mad and
#             v_mad feed the EKE_mad predictor in 05_compile_dataset.R.
#           - Coverage = fraction of box grid cells that contain water at a
#             level (from temperature for PHY and oxygen for BGC). Deep levels
#             only partly cover the box, so deep box means represent only the
#             cells deep enough to reach that level. All levels are kept; the
#             coverage tables are supplementary outputs.
#           - Sea ice: in the source files siconc/sithick are NaN (their
#             _FillValue) over open-water ocean cells and never defined over
#             land. NaN over ocean cells is therefore set to 0 (no ice) BEFORE
#             averaging, so that ice_conc is the fraction of the ocean part of
#             the box covered by ice. Ocean cells are those with a sea surface
#             height value. NOTE: the original analysis averaged only the
#             ice-covered cells and set whole-box NaN to 0 afterwards, which
#             overestimates concentration on partially ice-covered days.
#             [CONFIRM] the ocean-cell averaging used here.
#           - Every variable read is checked against a plausible physical range
#             (PLAUSIBLE_RANGE); values outside it (e.g., packed integers that
#             were never scaled) stop the script.
#           - Missing values are checked: within each site and box, the number
#             of valid cells at each level must be the same on every day
#             (structural land/seafloor mask only). Velocities can be masked in
#             a few more bottom cells than tracers; this is reported.
#           - The model levels are kept as provided (PHY and BGC levels
#             differ); 05_compile_dataset.R picks the nearest level of each
#             product to each nominal species depth.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

library(tidyverse)
library(ncdf4)

# -----------------------------------------------------------------------------
# 1. Script parameters
# -----------------------------------------------------------------------------
COP_DOWNLOAD_DIR <- file.path(DOWNLOAD_DIR, "copernicus")
LEGACY_DIRS <- c(phy = file.path(ENV_DIR, "Copernicus", "Physics Reanalysis"),
                 bgc = file.path(ENV_DIR, "Copernicus", "Biogeochemistry Hindcast"))

# Boxes to summarise, with the suffix added to the output file names.
BOXES <- list(
  main  = list(suffix = "",
               corners = transmute(SITES, site, lat_min = box_lat_min, lat_max = box_lat_max,
                                   lon_min = box_lon_min, lon_max = box_lon_max)),
  box40 = list(suffix = "_box40",
               corners = transmute(SITES, site, lat_min = box40_lat_min, lat_max = box40_lat_max,
                                   lon_min = box40_lon_min, lon_max = box40_lon_max))
)

# Source variable -> output name. Velocity spread is a MAD, others an SD.
VAR_NAMES <- c(thetao = "temperature", so = "salinity", uo = "u", vo = "v",
               zos = "ssh", mlotst = "mixed_layer", siconc = "ice_conc",
               sithick = "ice_thickness",
               chl = "chla", o2 = "o2", nppv = "productivity")
MAD_VARS   <- c("uo", "vo")
ICE_VARS   <- c("siconc", "sithick")
OCEAN_MASK_VAR <- "zos"                          # defines ocean cells for ice NA -> 0
COVERAGE_VAR   <- c(phy = "thetao", bgc = "o2")  # defines water cells per level

# Tolerance (degrees) when testing grid-cell centres against box corners.
BOX_TOL_DEG <- 1e-6

# Plausible physical ranges per source variable (same limits as in
# 01_download_copernicus.R). Any value outside stops the script, so packed
# integers or unit errors can never reach the outputs.
PLAUSIBLE_RANGE <- list(thetao = c(-3, 15), so = c(30, 36), uo = c(-3, 3), vo = c(-3, 3),
                        zos = c(-3, 3), mlotst = c(0, 2000), siconc = c(0, 1),
                        sithick = c(0, 10), chl = c(0, 50), o2 = c(0, 500),
                        nppv = c(-1, 500))

with_suffix <- function(file, suffix) sub("\\.csv$", paste0(suffix, ".csv"), file)

# -----------------------------------------------------------------------------
# 2. Locate input files for one site and product
# -----------------------------------------------------------------------------
# Climatology-window downloads are preferred; otherwise the original
# deployment-period portal files are used. The *hfdeg* (0.5 degree box) and
# *_ALLsurf* files in the same folders are deliberately not matched.
find_files <- function(s, p) {
  new <- list.files(COP_DOWNLOAD_DIR, pattern = sprintf("^%s_%s_.*\\.nc$", s, p),
                    full.names = TRUE)
  if (length(new) > 0) {
    message(sprintf("%s %s: %d downloaded file(s) in %s", s, p, length(new), COP_DOWNLOAD_DIR))
    return(new)
  }
  legacy <- list.files(LEGACY_DIRS[[p]],
                       pattern = sprintf("^%s40km_cmems_mod_glo_%s_.*\\.nc$", s, p),
                       full.names = TRUE)
  if (length(legacy) == 0) stop("No Copernicus ", p, " files found for site ", s)
  message(sprintf("%s %s: using deployment-period file(s): %s", s, p,
                  paste(basename(legacy), collapse = ", ")))
  legacy
}

# -----------------------------------------------------------------------------
# 3. Read helpers
# -----------------------------------------------------------------------------
# Convert a CF time axis ("<unit> since <origin>") to Date (daily means are
# stamped at 00:00 UTC of their day).
nc_dates <- function(nc) {
  units <- ncatt_get(nc, "time", "units")$value
  step <- c(seconds = 1, hours = 3600, days = 86400)[[word(units, 1)]]
  origin <- as.POSIXct(str_remove(units, "^\\w+ since "), tz = "UTC")
  as.Date(origin + ncvar_get(nc, "time") * step)
}

# Read one variable over the whole file (the download region). Coordinates are
# rounded so cells match exactly across files (ocean-mask join).
read_var <- function(nc, var) {
  dim_names <- map_chr(nc$var[[var]]$dim, "name")
  has_depth <- "depth" %in% dim_names
  if (!identical(dim_names, c("longitude", "latitude", if (has_depth) "depth", "time"))) {
    stop("Unexpected dimension order for ", var, ": ", paste(dim_names, collapse = ", "))
  }
  list(x     = ncvar_get(nc, var, collapse_degen = FALSE),   # keep degenerate dims
       lon   = round(as.vector(ncvar_get(nc, "longitude")), 5),
       lat   = round(as.vector(ncvar_get(nc, "latitude")), 5),
       has_depth = has_depth,
       # 2 decimals: portal (float32) and downloaded depths can differ in the
       # 3rd decimal (e.g., 163.165 vs 163.164 m) for the same model level.
       depth = if (has_depth) round(as.vector(ncvar_get(nc, "depth")), 2) else NA_real_,
       dates = nc_dates(nc))
}

# Crop a read variable to one box: values as (cell, [depth,] time) plus the
# cell coordinates.
crop_box <- function(r, bx) {
  i_lon <- which(r$lon >= bx$lon_min - BOX_TOL_DEG & r$lon <= bx$lon_max + BOX_TOL_DEG)
  i_lat <- which(r$lat >= bx$lat_min - BOX_TOL_DEG & r$lat <= bx$lat_max + BOX_TOL_DEG)
  if (length(i_lon) == 0 || length(i_lat) == 0) stop("Box ", bx$site, " contains no grid cells")
  x <- if (r$has_depth) r$x[i_lon, i_lat, , , drop = FALSE] else r$x[i_lon, i_lat, , drop = FALSE]
  d <- dim(x)
  list(values = array(x, dim = c(d[1] * d[2], d[-(1:2)])),   # cells first
       cells  = expand_grid(lat = r$lat[i_lat], lon = r$lon[i_lon]) %>% select(lon, lat))
}

# Column-wise (per day x level) box statistics on a cells x columns matrix.
# The SD uses the two-pass formula for numerical accuracy; MAD needs medians,
# so it is computed column by column only for the velocity variables.
box_stats <- function(m, spread = c("sd", "mad")) {
  spread <- match.arg(spread)
  n  <- colSums(!is.na(m))
  mu <- colSums(m, na.rm = TRUE) / n
  sp <- if (spread == "sd") {
    sqrt(colSums((m - rep(mu, each = nrow(m)))^2, na.rm = TRUE) / (n - 1))
  } else {
    apply(m, 2, mad, na.rm = TRUE)
  }
  mu[n == 0] <- NA
  sp[n < 2] <- NA      # spread undefined with fewer than two cells
  tibble(mean = mu, spread = sp, n_valid = n)
}

# -----------------------------------------------------------------------------
# 4. Summarise all files of one site and product, for both boxes
# -----------------------------------------------------------------------------
# Each file is read once and cropped to every box. Returns long-format
# statistics: box, date, depth_m (NA for 2-D fields), var, mean, spread,
# n_valid, n_cells.
summarise_site_product <- function(s, p) {
  files <- find_files(s, p)
  wanted <- COPERNICUS[[p]]$vars
  site_boxes_list <- map(BOXES, ~ filter(.x$corners, site == s))

  # Ocean mask over the download grid (cells with sea surface height on any
  # day) for the ice NA -> 0 rule; built before any ice variable is read.
  ocean <- NULL
  if (p == "phy") {
    zos_files <- keep(files, function(f) {
      nc <- nc_open(f); on.exit(nc_close(nc)); OCEAN_MASK_VAR %in% names(nc$var)
    })
    ocean <- map(zos_files, function(f) {
      nc <- nc_open(f); on.exit(nc_close(nc))
      r <- read_var(nc, OCEAN_MASK_VAR)
      expand_grid(lat = r$lat, lon = r$lon) %>%
        select(lon, lat) %>%
        mutate(is_ocean = as.vector(apply(!is.na(r$x), c(1, 2), any)))
    }) %>%
      bind_rows() %>%
      group_by(lon, lat) %>%
      summarise(is_ocean = any(is_ocean), .groups = "drop")
  }

  map_dfr(files, function(f) {
    nc <- nc_open(f)
    on.exit(nc_close(nc))
    map_dfr(intersect(wanted, names(nc$var)), function(v) {
      r <- read_var(nc, v)
      # Physical-range check: catches packed (unscaled) integers or wrong
      # units before they reach any statistic. Small tolerance (0.1% of the
      # range width) for int16 rounding, e.g. siconc 1.00002.
      rng <- suppressWarnings(range(r$x, na.rm = TRUE))
      lim <- PLAUSIBLE_RANGE[[v]]
      tol <- 1e-3 * diff(lim)
      if (all(is.finite(rng)) && (rng[1] < lim[1] - tol || rng[2] > lim[2] + tol)) {
        stop(sprintf("%s in %s: values %.4g..%.4g outside plausible range %g..%g",
                     v, basename(f), rng[1], rng[2], lim[1], lim[2]))
      }
      imap_dfr(site_boxes_list, function(bx, box_name) {
        b <- crop_box(r, bx)
        vals <- b$values
        if (v %in% ICE_VARS) {
          # NaN over ocean = open water (source _FillValue): set to 0 so the
          # box mean is the ice-covered fraction of the ocean part of the box.
          is_ocean <- left_join(b$cells, ocean, by = c("lon", "lat"))$is_ocean
          if (anyNA(is_ocean)) stop("Ocean mask does not cover all ", v, " cells at ", s)
          vals[is.na(vals) & array(is_ocean, dim = dim(vals))] <- 0
        }
        stats <- box_stats(matrix(vals, ncol = prod(dim(vals)[-1])),
                           if (v %in% MAD_VARS) "mad" else "sd")
        # Columns run over depth fastest, then time (array order).
        n_depth <- length(r$depth)
        stats %>%
          mutate(box     = box_name,
                 date    = rep(r$dates, each = n_depth),
                 depth_m = rep(r$depth, times = length(r$dates)),
                 var     = v,
                 n_cells = nrow(vals),
                 .before = 1)
      })
    })
  }) %>%
    mutate(site = s, product = p, .before = 1)
}

long_all <- map_dfr(SITE_ORDER, function(s) {
  bind_rows(summarise_site_product(s, "phy"), summarise_site_product(s, "bgc"))
})

# -----------------------------------------------------------------------------
# 5. Checks, coverage and contract tables for one box
# -----------------------------------------------------------------------------
to_wide <- function(d, spread_suffix) {
  d <- mutate(d, name = unname(VAR_NAMES[var]))
  bind_rows(transmute(d, site, date, depth_m, name, value = mean),
            transmute(d, site, date, depth_m,
                      name = paste0(name, spread_suffix[name]), value = spread)) %>%
    pivot_wider(names_from = name, values_from = value)
}

# After dropping non-existent levels, box means must be complete. Spreads may
# be NA only where a level has a single valid cell (spread undefined).
check_complete <- function(d, mean_cols, label) {
  bad <- d %>% filter(if_any(all_of(mean_cols), is.na))
  if (nrow(bad) > 0) stop(label, ": ", nrow(bad), " rows with NA box means")
  n_sd_na <- d %>% summarise(across(ends_with(c("_sd", "_mad")), ~ sum(is.na(.x))))
  if (any(n_sd_na > 0)) {
    message(label, ": NA spreads (levels with < 2 valid cells):")
    print(n_sd_na)
  }
}

order_rows <- function(d) {
  d %>% arrange(factor(site, levels = SITE_ORDER), date, across(any_of("depth_m")))
}

process_box <- function(box_name) {
  suffix <- BOXES[[box_name]]$suffix
  message("\n===== Box: ", box_name, " =====")
  long <- filter(long_all, box == box_name)

  # Overlapping files (e.g., a re-download next to an old file) would
  # silently duplicate days, so stop instead.
  dups <- long %>% count(site, product, var, date, depth_m) %>% filter(n > 1)
  if (nrow(dups) > 0) stop("Duplicate site/variable/day/level records in the input files")

  # (a) Within a site, the number of valid cells at a level must not change
  #     from day to day (that would be missing model data, not the fixed
  #     land/seafloor mask). Ice variables are checked after NA -> 0 and must
  #     then match the ocean mask.
  mask_check <- long %>%
    group_by(site, product, var, depth_m) %>%
    summarise(n_min = min(n_valid), n_max = max(n_valid), n_cells = first(n_cells),
              .groups = "drop")
  unstable <- filter(mask_check, n_min != n_max)
  if (nrow(unstable) > 0) {
    print(unstable)
    stop("Number of valid cells varies over time for the variables/levels above")
  }

  # (b) Levels where a variable has fewer valid cells than the coverage
  #     variable (e.g., velocities masked in bottom-adjacent cells).
  ref_counts <- mask_check %>%
    filter(var == COVERAGE_VAR[product]) %>%
    select(site, product, depth_m, n_ref = n_min)
  fewer <- mask_check %>%
    filter(!is.na(depth_m)) %>%
    inner_join(ref_counts, by = c("site", "product", "depth_m")) %>%
    filter(n_min != n_ref)
  if (nrow(fewer) > 0) {
    message("Variables with fewer valid cells than the coverage variable ",
            "(structural, constant in time):")
    print(fewer)
  }

  # Coverage table (with the number of grid cells in the box).
  coverage <- mask_check %>%
    filter(var == COVERAGE_VAR[product]) %>%
    transmute(site, product, depth_m, n_cells, coverage = n_min / n_cells)
  write_csv(coverage, with_suffix(F_COP_COVERAGE, suffix))

  # Contract tables; levels with no water anywhere in the box are dropped
  # (structural, not data gaps).
  cov <- select(coverage, site, product, depth_m, coverage)
  phy_depth <- long %>%
    filter(product == "phy", !is.na(depth_m)) %>%
    to_wide(c(temperature = "_sd", salinity = "_sd", u = "_mad", v = "_mad")) %>%
    left_join(filter(cov, product == "phy") %>% select(-product), by = c("site", "depth_m")) %>%
    filter(coverage > 0) %>%
    select(date, site, depth_m, temperature, temperature_sd, salinity, salinity_sd,
           u, v, u_mad, v_mad, coverage)
  bgc_depth <- long %>%
    filter(product == "bgc") %>%
    to_wide(c(chla = "_sd", o2 = "_sd", productivity = "_sd")) %>%
    left_join(filter(cov, product == "bgc") %>% select(-product), by = c("site", "depth_m")) %>%
    filter(coverage > 0) %>%
    select(date, site, depth_m, chla, chla_sd, o2, o2_sd, productivity,
           productivity_sd, coverage)
  surface <- long %>%
    filter(product == "phy", is.na(depth_m)) %>%
    to_wide(c(ssh = "_sd", mixed_layer = "_sd", ice_conc = "_sd", ice_thickness = "_sd")) %>%
    select(date, site, ssh, ssh_sd, mixed_layer, mixed_layer_sd, ice_conc,
           ice_conc_sd, ice_thickness)

  check_complete(phy_depth, c("temperature", "salinity", "u", "v"), "PHY")
  check_complete(bgc_depth, c("chla", "o2", "productivity"), "BGC")
  check_complete(surface, c("ssh", "mixed_layer", "ice_conc", "ice_thickness"), "Surface")

  write_csv(order_rows(phy_depth), with_suffix(F_COP_PHY_DEPTH, suffix))
  write_csv(order_rows(bgc_depth), with_suffix(F_COP_BGC_DEPTH, suffix))
  write_csv(order_rows(surface), with_suffix(F_COP_SURFACE, suffix))

  # Console summary for the log.
  message("Days per site: ")
  print(surface %>% group_by(site) %>%
          summarise(days = n(), first = min(date), last = max(date)))
  message("Grid cells in the box (surface level):")
  print(coverage %>% group_by(site, product) %>% summarise(n_cells = first(n_cells),
                                                           .groups = "drop") %>%
          pivot_wider(names_from = product, values_from = n_cells))
  message("Coverage at levels deeper than 800 m:")
  print(coverage %>% filter(depth_m > 800) %>% select(-n_cells) %>%
          pivot_wider(names_from = site, values_from = coverage), n = Inf)
}

# -----------------------------------------------------------------------------
# 6. Run for both boxes
# -----------------------------------------------------------------------------
walk(names(BOXES), process_box)
