# =============================================================================
# Script  : 05_compile_dataset.R
# Section : 05_environmental_data
# Purpose : Compile the daily environmental predictor table for each site
#           (one row per site x day, one column per variable x nominal depth),
#           derive eddy kinetic energy (EKE) and lagged chlorophyll-a /
#           productivity, and join daily acoustic presence on recording days
#           to produce the model input table.
# Inputs  : final_scripts_odontocetes/config.R (SPECIES depths, ALL_DEPTHS_M,
#           CHL_PROD_MAX_DEPTH_M, LAG_MONTHS, LAG_WINDOW_DAYS, LAG_VARS)
#           F_COP_PHY_DEPTH, F_COP_BGC_DEPTH, F_COP_SURFACE, F_COP_COVERAGE
#             (main +/- 20 km box; 04_extract_copernicus.R)
#           F_FSLE_DAILY (03_extract_aviso_fsle.R)
#           F_PRESENCE_DAILY (04_acoustic_detection; optional - if absent,
#             the model input table is not written)
# Outputs : F_ENV_DAILY   env_daily_wide.csv
#             date, site, julian_day, ice_conc, ice_thickness, SSH,
#             mixed_layer, FSLE, fsle_orient, temperature_<d>, salinity_<d>,
#             o2_<d>, EKE_<d>, EKE_mad_<d>, chla_<d>, productivity_<d>
#             (d <= CHL_PROD_MAX_DEPTH_M only), chla_lag<N>,
#             productivity_lag<N>; <d> = nominal depth (m) in ALL_DEPTHS_M
#           F_MODEL_INPUT model_input_daily.csv
#             recording days only: F_ENV_DAILY joined to presence
#             (rec_minutes + presence minutes per species)
#           outputs/05_environmental_data/depth_level_map.csv
#             nominal depth -> model level used for PHY and BGC
#           outputs/05_environmental_data/depth_level_fallback.csv
#             nominal depths whose nearest level has no water in the main
#             box, with the deepest available level used instead
# Author  : Natalie Posdaljian
# Date    : 2026-11-18
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Depths: for each nominal depth in ALL_DEPTHS_M (species surface
#             and dive depths), the nearest model level of each product is
#             used (PHY and BGC have different vertical grids). The mapping is
#             saved in depth_level_map.csv. Where the nearest level has no
#             water in the main box at a site (deep levels on the coarse BGC
#             grid), the deepest level with water at that site is used; these
#             substitutions are listed in depth_level_fallback.csv (site,
#             nominal_depth_m, product, level_used_m, n_cells_with_water,
#             reason).
#           - EKE (cm^2 s^-2) keeps the definition of the original
#             GetCopernicus.R: u' = daily box-mean u minus the reference mean
#             u (same for v), EKE = 0.5 (u'^2 + v'^2) x 1e4. Reference mean =
#             the 10-year mean (full extracted climatology window) per site and
#             level.
#           - EKE_mad (cm^2 s^-2) = 0.5 (u_mad^2 + v_mad^2) x 1e4: kinetic
#             energy of the within-box spatial velocity variability (daily
#             spatial MADs of u and v), with no mean subtracted. This departs
#             from the original GetCopernicus.R, which applied the EKE anomaly
#             formula to the MADs, so that both very uneven and very uniform
#             days scored high.
#           - Lagged chl-a and productivity (surface level only): lag<N> = mean
#             of the LAG_WINDOW_DAYS-day window ending N x LAG_WINDOW_DAYS days
#             before each date (N in LAG_MONTHS). A window with any missing day
#             gives NA. Lags need data before each deployment: with
#             deployment-period files only, early recording days have NA lags
#             (reported); the climatology-window downloads remove these NAs.
#           - Box: all predictors come from the MAIN +/- 20 km box.
#           - Variables dropped from the original compiled dataset by design:
#             AAO, deseasoned STL/anomaly columns (deseasoning is done in 06),
#             bathymetry, ice_diff, ice_regime, spatial SD columns, and lags of
#             variables other than chl-a and productivity.
#           - Every modeled variable must be complete on recording days; the
#             script stops otherwise. The only allowed NAs are lags before the
#             start of the input series, which are reported.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

library(tidyverse)

F_DEPTH_MAP <- out_path("05_environmental_data", "depth_level_map.csv")

# -----------------------------------------------------------------------------
# 1. Read the extracted environmental tables
# -----------------------------------------------------------------------------
# Inputs from earlier scripts must exist; a missing file means a step was
# skipped, so stop with the name of the script that makes it.
need <- c(F_COP_PHY_DEPTH = "04_extract_copernicus.R",
          F_COP_BGC_DEPTH = "04_extract_copernicus.R",
          F_COP_SURFACE   = "04_extract_copernicus.R",
          F_COP_COVERAGE  = "04_extract_copernicus.R",
          F_FSLE_DAILY    = "03_extract_aviso_fsle.R")
for (f in names(need)) {
  if (!file.exists(get(f))) stop(get(f), " not found: run ", need[[f]], " first")
}

phy     <- read_csv(F_COP_PHY_DEPTH, show_col_types = FALSE)
bgc     <- read_csv(F_COP_BGC_DEPTH, show_col_types = FALSE)
surface <- read_csv(F_COP_SURFACE,   show_col_types = FALSE)
fsle    <- read_csv(F_FSLE_DAILY,    show_col_types = FALSE)

# -----------------------------------------------------------------------------
# 2. Nominal depth -> nearest model level
# -----------------------------------------------------------------------------
# Model levels are identical at all sites (deep levels only exist where the box
# is deep enough, so use the union).
nearest_level <- function(levels, d) levels[which.min(abs(levels - d))]
phy_levels <- sort(unique(phy$depth_m))
bgc_levels <- sort(unique(bgc$depth_m))

depth_map <- tibble(nominal_depth_m = ALL_DEPTHS_M) %>%
  mutate(phy_level_m = map_dbl(nominal_depth_m, ~ nearest_level(phy_levels, .x)),
         bgc_level_m = map_dbl(nominal_depth_m, ~ nearest_level(bgc_levels, .x)))
write_csv(depth_map, F_DEPTH_MAP)
message("Nominal depth -> model level:")
print(depth_map)

# Site-specific levels. Every site needs a level for every nominal depth, but
# in the +/- 20 km main box some deep levels hold no water at some sites (the
# coarse 0.25 degree BGC grid has only 3-6 cells per box; e.g., no BGC cell at
# KGI reaches below ~800 m). For such site x nominal depth x product
# combinations the DEEPEST level that has water in the main box at that site
# is used instead (decision N. Posdaljian, 2026-11-18), and the substitution is
# written to F_LEVEL_FALLBACK.
F_LEVEL_FALLBACK <- out_path("05_environmental_data", "depth_level_fallback.csv")
coverage <- read_csv(F_COP_COVERAGE, show_col_types = FALSE) %>%
  mutate(product = toupper(product), n_cells_with_water = round(coverage * n_cells))

site_levels <- depth_map %>%
  pivot_longer(-nominal_depth_m, names_to = "product", values_to = "level_nearest_m") %>%
  mutate(product = toupper(str_remove(product, "_level_m"))) %>%
  expand_grid(site = SITE_ORDER, .) %>%
  rowwise() %>%
  mutate(available = list(coverage$depth_m[coverage$site == site & coverage$product == product &
                                             coverage$coverage > 0]),
         has_level = level_nearest_m %in% available,
         level_used_m = if (has_level) level_nearest_m else max(available)) %>%
  ungroup() %>%
  select(-available)

fallback <- site_levels %>%
  filter(!has_level) %>%
  left_join(select(coverage, site, product, depth_m, n_cells_with_water),
            by = c("site", "product", "level_used_m" = "depth_m")) %>%
  transmute(site, nominal_depth_m, product, level_used_m, n_cells_with_water,
            reason = sprintf("no water at nearest level %.2f m in the main box; deepest level with water used",
                             level_nearest_m))
write_csv(fallback, F_LEVEL_FALLBACK)
if (nrow(fallback) > 0) {
  message("Nominal depths without water in the main box (deepest available level used):")
  print(fallback %>% select(-reason))
}
level_map <- function(prod) {
  site_levels %>% filter(product == prod) %>% select(site, nominal_depth_m, depth_m = level_used_m)
}

# -----------------------------------------------------------------------------
# 3. EKE at the selected PHY levels
# -----------------------------------------------------------------------------
# EKE (cm^2 s^-2; 1e4 converts m^2 s^-2): kinetic energy of the daily box-mean
# velocity anomaly from the 10-year mean (full extracted climatology window)
# per site and level.
# EKE_mad: kinetic energy of the within-box spatial velocity variability,
# 0.5 (u_mad^2 + v_mad^2) x 1e4, with u_mad/v_mad the daily spatial MADs; no
# mean is subtracted. (The original GetCopernicus.R squared the deviation of
# the MAD from its mean, so both very uneven and very uniform days scored
# high; redefined 2026-11-18, decision N. Posdaljian.)
eke <- phy %>%
  inner_join(level_map("PHY"), by = c("site", "depth_m"), relationship = "many-to-many") %>%
  group_by(site, nominal_depth_m) %>%
  mutate(EKE     = 0.5 * ((u - mean(u))^2 + (v - mean(v))^2) * 1e4,
         EKE_mad = 0.5 * (u_mad^2 + v_mad^2) * 1e4) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 4. Depth-resolved variables to wide format
# -----------------------------------------------------------------------------
# Column suffix = nominal depth (e.g., temperature_375), so names do not
# depend on the product's model level.
phy_wide <- eke %>%
  select(date, site, nominal_depth_m, temperature, salinity, EKE, EKE_mad) %>%
  pivot_wider(names_from = nominal_depth_m,
              values_from = c(temperature, salinity, EKE, EKE_mad))

bgc_wide <- bgc %>%
  inner_join(level_map("BGC"), by = c("site", "depth_m"), relationship = "many-to-many") %>%
  select(date, site, nominal_depth_m, chla, o2, productivity) %>%
  pivot_longer(c(chla, o2, productivity), names_to = "var") %>%
  # chl-a and productivity are ~0 below the euphotic zone: keep them only at
  # nominal depths <= CHL_PROD_MAX_DEPTH_M.
  filter(var == "o2" | nominal_depth_m <= CHL_PROD_MAX_DEPTH_M) %>%
  pivot_wider(names_from = c(var, nominal_depth_m), values_from = value)

# -----------------------------------------------------------------------------
# 5. Lagged chl-a and productivity (surface level)
# -----------------------------------------------------------------------------
# On a gap-free daily calendar per site, the trailing window mean (all days
# required; stats::filter returns NA otherwise) is shifted by
# N x LAG_WINDOW_DAYS days, so lag<N> on day t averages days
# t - N*W - (W - 1) ... t - N*W.
surface_bgc_level <- depth_map$bgc_level_m[depth_map$nominal_depth_m == 0]

lags <- bgc %>%
  filter(depth_m == surface_bgc_level) %>%
  select(site, date, all_of(LAG_VARS)) %>%
  group_by(site) %>%
  complete(date = seq(min(date), max(date), by = "day")) %>%
  arrange(date, .by_group = TRUE)
for (v in LAG_VARS) {
  win <- paste0(v, "_win")
  lags <- mutate(lags, !!win := as.numeric(
    stats::filter(.data[[v]], rep(1 / LAG_WINDOW_DAYS, LAG_WINDOW_DAYS), sides = 1)))
  for (n in LAG_MONTHS) {
    lags <- mutate(lags, !!paste0(v, "_lag", n) := dplyr::lag(.data[[win]], n * LAG_WINDOW_DAYS))
  }
}
lags <- lags %>%
  ungroup() %>%
  select(site, date, matches("_lag\\d+$"))

# -----------------------------------------------------------------------------
# 6. Assemble the daily environmental table
# -----------------------------------------------------------------------------
# The PHY surface table defines the site-days; BGC and FSLE are joined to it
# and any day missing from them is reported.
env <- surface %>%
  transmute(date, site, ice_conc, ice_thickness, SSH = ssh, mixed_layer) %>%
  left_join(select(fsle, date, site, FSLE, fsle_orient), by = c("date", "site")) %>%
  left_join(phy_wide, by = c("date", "site")) %>%
  left_join(bgc_wide, by = c("date", "site")) %>%
  left_join(lags, by = c("date", "site")) %>%
  mutate(julian_day = yday(date), .after = site)

# Column order of the data contract.
depth_cols <- function(var, depths) paste0(var, "_", depths)
chl_depths <- ALL_DEPTHS_M[ALL_DEPTHS_M <= CHL_PROD_MAX_DEPTH_M]
lag_cols <- names(lags)[str_detect(names(lags), "_lag\\d+$")]
env <- env %>%
  select(date, site, julian_day, ice_conc, ice_thickness, SSH, mixed_layer, FSLE,
         fsle_orient,
         all_of(c(depth_cols("temperature", ALL_DEPTHS_M), depth_cols("salinity", ALL_DEPTHS_M),
                  depth_cols("o2", ALL_DEPTHS_M), depth_cols("EKE", ALL_DEPTHS_M),
                  depth_cols("EKE_mad", ALL_DEPTHS_M), depth_cols("chla", chl_depths),
                  depth_cols("productivity", chl_depths), lag_cols))) %>%
  arrange(factor(site, levels = SITE_ORDER), date)

write_csv(env, F_ENV_DAILY)
message(sprintf("Wrote %s: %d rows x %d columns", F_ENV_DAILY, nrow(env), ncol(env)))

# NA report for the full table (all days, not only recording days).
na_env <- env %>%
  group_by(site) %>%
  summarise(across(-date, ~ sum(is.na(.x)))) %>%
  pivot_longer(-site, names_to = "variable", values_to = "n_na") %>%
  filter(n_na > 0)
if (nrow(na_env) > 0) {
  message("NA counts in ", basename(F_ENV_DAILY), " (all days):")
  print(pivot_wider(na_env, names_from = site, values_from = n_na), n = Inf)
}

# -----------------------------------------------------------------------------
# 7. Join presence on recording days -> model input
# -----------------------------------------------------------------------------
# Presence comes from 04_acoustic_detection; if it is not there yet, stop here
# with a message rather than an error so the environmental table is still made.
if (!file.exists(F_PRESENCE_DAILY)) {
  message("Presence table not found (", F_PRESENCE_DAILY, "): skipping the ",
          "model input table. Run 04_acoustic_detection first, then re-run this script.")
} else {
  presence <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) %>%
    filter(rec_minutes > 0)          # recording days only

  # Recording days without environmental data cannot be modeled; list them.
  no_env <- anti_join(presence, env, by = c("date", "site"))
  if (nrow(no_env) > 0) {
    warning(nrow(no_env), " recording day(s) have no environmental data and are dropped: ",
            paste(sprintf("%s %s (%.1f min)", no_env$site, no_env$date, no_env$rec_minutes),
                  collapse = "; "))
  }

  model_input <- inner_join(env, presence, by = c("date", "site"))

  # Lag NAs are allowed only when the lag window starts before the input
  # series; report them per site and lag.
  series_start <- bgc %>% group_by(site) %>% summarise(first_day = min(date))
  lag_na <- model_input %>%
    select(site, date, all_of(lag_cols)) %>%
    pivot_longer(all_of(lag_cols), names_to = "lag", values_to = "value") %>%
    left_join(series_start, by = "site") %>%
    mutate(n_lag = as.integer(str_extract(lag, "\\d+$")),
           window_start = date - n_lag * LAG_WINDOW_DAYS - (LAG_WINDOW_DAYS - 1))
  unexplained <- filter(lag_na, is.na(value), window_start >= first_day)
  if (nrow(unexplained) > 0) {
    stop(nrow(unexplained), " lag values are NA although their window is inside the input series")
  }
  lag_coverage <- lag_na %>%
    group_by(site, lag) %>%
    summarise(recording_days = n(), n_na = sum(is.na(value)), .groups = "drop") %>%
    mutate(pct_available = round(100 * (1 - n_na / recording_days), 1))
  message("Lag availability on recording days (NA = window before start of input series):")
  print(pivot_wider(select(lag_coverage, site, lag, pct_available),
                    names_from = site, values_from = pct_available), n = Inf)

  # Every other modeled variable must be complete on recording days.
  non_lag <- setdiff(names(env), c("date", "site", lag_cols))
  na_model <- model_input %>%
    summarise(across(all_of(non_lag), ~ sum(is.na(.x)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_na") %>%
    filter(n_na > 0)
  if (nrow(na_model) > 0) {
    print(na_model, n = Inf)
    stop("NA in modeled variables on recording days (see above)")
  }

  write_csv(model_input, F_MODEL_INPUT)
  message(sprintf("Wrote %s: %d rows x %d columns", F_MODEL_INPUT,
                  nrow(model_input), ncol(model_input)))
  print(count(model_input, site, name = "recording_days"))
}
