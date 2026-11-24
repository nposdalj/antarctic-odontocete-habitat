# =============================================================================
# Script  : 00_gam_functions.R
# Section : 2.6 Statistical analyses (also used by 3.2 GAM results)
# Purpose : Functions shared by the three GAM stages and the 08 result
#           scripts: variable naming/labels, deseasoning (10-yr monthly
#           climatology + periodic spline), residual-ACF bin length, binning,
#           ice filter, VIF screening, single-term theme screening, backward
#           elimination, model documentation (diagnostics, per-term deviance,
#           smooth shape classes) and the common on-disk result format.
# Inputs  : none (functions only); requires config.R and 00_gam_parameters.R
#           to be sourced first.
# Outputs : none
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Rewritten from GAMs/DEPTH_BINNED_Automated_GAM_tweedie_variability_
#           deSeasoned_allSpecies.R. Changes relative to that script are noted
#           where the function is defined ("Change:").
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(car)
  library(gratia)
  library(patchwork)
})

# =============================================================================
# 1. Names and labels
# =============================================================================

# Split a candidate column name into its parts, e.g.
#   "temperature_mid_anom" -> base "temperature", layer "mid", anomaly TRUE
#   "temperature_375_anom" -> base "temperature", depth 375, anomaly TRUE
#   "chla_lag3"            -> base "chla", lag 3
# Used to assign themes, depths/layers and labels from the column name alone.
# Layer variables get the theme water_mass_<layer> (temperature, salinity,
# oxygen) or eddy_<layer> (EKE); other variables use THEME_OF_BASE.
parse_variable <- function(variable) {
  is_anom <- grepl("_anom$", variable)
  stem    <- sub("_anom$", "", variable)
  lag     <- suppressWarnings(as.integer(sub("^.*_lag([0-9]+)$", "\\1", stem)))
  lag[!grepl("_lag[0-9]+$", stem)] <- NA_integer_
  stem    <- sub("_lag[0-9]+$", "", stem)
  layer_pat <- paste0("_(", paste(LAYERS$layer, collapse = "|"), ")$")
  layer   <- ifelse(grepl(layer_pat, stem), sub(paste0("^.*", layer_pat), "\\1", stem),
                    NA_character_)
  stem    <- sub(layer_pat, "", stem)
  depth   <- suppressWarnings(as.numeric(sub("^.*_([0-9]+)$", "\\1", stem)))
  depth[!grepl("_[0-9]+$", stem)] <- NA_real_
  base    <- sub("_[0-9]+$", "", stem)
  theme   <- case_when(
    !is.na(layer) & base %in% WATER_MASS_BASES ~ paste0("water_mass_", layer),
    !is.na(layer) & base == "EKE"              ~ paste0("eddy_", layer),
    base == ICE_MODEL_VAR                      ~ "sea ice",
    TRUE                                       ~ unname(THEME_OF_BASE[base]))
  tibble(variable = variable, base = base, layer = layer, depth_m = depth,
         lag_months = lag, is_anom = is_anom, theme = theme)
}

# Human-readable variable label with units, for figures and tables.
# Change: replaces nameVar(); unicode escapes avoid the degree-symbol encoding
# problem of the original scripts.
variable_label <- function(variable) {
  base_labels <- c(
    temperature  = "Temperature (°C)",
    salinity     = "Salinity (psu)",
    o2           = "Oxygen (mmol m-3)",
    EKE          = "EKE (cm2 s-2)",          # [CONFIRM] EKE units
    EKE_mad      = "EKE MAD (cm2 s-2)",
    chla         = if (CHL_LOG10) "log10 chl-a (mg m-3)" else "Chl-a (mg m-3)",
    SSH          = "SSH (m)",
    mixed_layer  = "Mixed layer depth (m)",
    FSLE         = "FSLE (d-1)",
    fsle_orient  = "FSLV orientation (°; cyclic, -90 = 90)",
    ice_conc     = "Sea-ice concentration"
  )
  layer_labels <- setNames(paste0(", ", tolower(sub(" [(].*", "", LAYERS$label)),
                                  " (", LAYERS$top_m, "-", LAYERS$bottom_m, " m)"),
                           LAYERS$layer)
  p <- parse_variable(variable)
  lab <- unname(base_labels[p$base])
  lab[is.na(lab)] <- p$base[is.na(lab)]
  prefix <- ifelse(p$is_anom, "Anomaly: ", "")
  std_chl <- p$is_anom & p$base == "chla" & CHL_STANDARDIZE   # in SD units
  prefix[std_chl] <- "Standardized anomaly: "
  lab[std_chl] <- "log10 chl-a (SD)"
  depth  <- ifelse(!is.na(p$layer), unname(layer_labels[p$layer]),
                   ifelse(is.na(p$depth_m), "",
                          ifelse(p$depth_m == 0, ", surface", paste0(", ", p$depth_m, " m"))))
  lag    <- ifelse(is.na(p$lag_months), "", paste0(", ", p$lag_months, "-month lag"))
  paste0(prefix, lab, depth, lag)
}

species_label <- function(code) vapply(code, function(x) SPECIES[[x]]$name, character(1))
model_id      <- function(species, site) paste0(species, "_", site)

# =============================================================================
# 2. Candidate predictor columns per species
# =============================================================================

# Source columns of the candidates for one species: surface variables,
# layer-averaged water-column variables for the species' layers
# (species_layers(), config.R; oxygen and EKE only below the surface layer)
# and the 3-month chl-a lag. Sea ice is not a habitat-model candidate (it has
# its own model, 04_sea_ice_models.R).
# Change: layer averages replace species-specific depths; one common 3-month
# chl-a lag; same-day chl-a and productivity dropped; SSH and EKE_mad dropped
# (decisions 2026-11-19).
layer_columns <- function(species, bases) {
  layers <- species_layers(species)
  out <- as.vector(outer(bases, layers, paste, sep = "_"))
  out[!(sub("_[a-z]+$", "", out) %in% SUBSURFACE_ONLY & grepl("_surface$", out))]
}

raw_candidates <- function(species) {
  c(SURFACE_VARS, layer_columns(species, LAYER_BASES), LAG_COLUMNS)
}

# Columns that enter as anomalies for one species: everything except the
# mesoscale and front variables (EKE_<layer>, FSLE, fsle_orient).
anomaly_sources <- function(species) {
  c(ANOM_SURFACE_VARS, layer_columns(species, ANOM_LAYER_BASES), LAG_COLUMNS)
}

# =============================================================================
# 3. Deseasoning: 10-yr monthly climatology + periodic cubic spline
# =============================================================================

# Periodic cubic spline through the 12 monthly means, anchored mid-month, so
# the seasonal baseline is smooth across month boundaries (including Dec-Jan)
# rather than a step per calendar month. splinefun(method = "periodic")
# evaluates any day of year modulo the 365-day period.
# Change: the climatology comes from the full climatology window
# (env_daily_wide.csv), not from the deployment's own binned data.
fit_periodic_climatology <- function(month, value) {
  monthly <- tapply(value, factor(month, levels = 1:12), mean, na.rm = TRUE)
  ok <- is.finite(monthly)
  if (sum(ok) < 2) return(NULL)
  x <- MID_MONTH_YDAY[ok]
  y <- as.numeric(monthly[ok])
  stats::splinefun(c(x, x[1] + 365), c(y, y[1]), method = "periodic")
}

# Monthly climatology table (site x variable x month) and coverage check.
# Production runs stop if a variable lacks CLIM_MIN_YEARS distinct years or
# any calendar month inside the site's climatology window; the test override
# only downgrades this to a loud warning.
build_climatology <- function(env_daily, variables) {
  windows <- SITES |>
    transmute(site,
              win_start = as.Date(sprintf("%d-01-01", as.integer(format(rec_start, "%Y")) -
                                            CLIM_HALF_WINDOW_YEARS)),
              win_end   = as.Date(sprintf("%d-12-31", as.integer(format(rec_end, "%Y")) +
                                            CLIM_HALF_WINDOW_YEARS)))
  env_win <- env_daily |>
    inner_join(windows, by = "site") |>
    filter(date >= win_start, date <= win_end) |>
    mutate(month = as.integer(format(date, "%m")),
           year  = as.integer(format(date, "%Y")))

  long <- env_win |>
    select(site, date, month, year, all_of(variables)) |>
    pivot_longer(all_of(variables), names_to = "variable", values_to = "value") |>
    filter(!is.na(value))

  coverage <- long |>
    group_by(site, variable) |>
    summarise(n_years = n_distinct(year), n_months = n_distinct(month),
              first_date = min(date), last_date = max(date), .groups = "drop")
  short <- coverage |> filter(n_years < CLIM_MIN_YEARS | n_months < 12)
  if (nrow(short) > 0) {
    msg <- paste0(nrow(short), " site x variable climatologies have < ", CLIM_MIN_YEARS,
                  " years or < 12 months of data (e.g., ",
                  paste(head(paste0(short$site, ":", short$variable, " (", short$n_years,
                                    " yr, ", short$n_months, " mo)"), 3), collapse = "; "), ").")
    if (!ALLOW_SHORT_CLIMATOLOGY) {
      stop(msg, " The deseasoning requires the 10-yr climatology window (CLIM_HALF_WINDOW_YEARS).")
    }
    warning("\n", strrep("!", 70), "\nTEST MODE: ", msg,
            "\nAnomalies are NOT relative to a 10-yr climatology. Do not report these results.\n",
            strrep("!", 70), call. = FALSE, immediate. = TRUE)
  }

  monthly <- long |>
    group_by(site, variable, month) |>
    summarise(monthly_mean = mean(value), n_days = n(), n_years = n_distinct(year),
              .groups = "drop")

  # Share of daily variance explained by the seasonal cycle over the window:
  # 1 - var(daily - spline climatology) / var(daily). Reported only (every
  # anomaly source is deseasoned; see 00_gam_parameters.R section 5).
  seasonal <- long |>
    group_by(site, variable) |>
    group_modify(function(d, key) {
      m <- monthly |> filter(site == key$site, variable == key$variable)
      f <- fit_periodic_climatology(m$month, m$monthly_mean)
      frac <- if (is.null(f) || var(d$value) == 0) NA_real_ else
        1 - var(d$value - f(as.integer(format(d$date, "%j")))) / var(d$value)
      tibble(seasonal_var_frac = frac)
    }) |>
    ungroup() |>
    mutate(version_used = "anomaly")
  list(monthly = monthly, coverage = coverage, seasonal = seasonal)
}

# Candidate columns of one species (same at every site): raw candidates, with
# every anomaly source replaced by its anomaly (raw OR anomaly, never both).
# `site` and `seasonal` are kept so callers need not change; the choice no
# longer depends on them.
site_candidates <- function(species, site, seasonal) {
  use_anom <- intersect(raw_candidates(species), anomaly_sources(species))
  out <- c(setdiff(raw_candidates(species), use_anom), paste0(use_anom, "_anom"))
  out[!sub("_anom$", "", out) %in% EXCLUDE_VARS]
}

# Add <var>_anom columns: daily value minus the spline climatology evaluated at
# that day of year (per site and variable). Anomalies are computed daily and
# then averaged within bins (the mean of anomalies equals the anomaly of the
# bin mean, but the baseline is evaluated at each day instead of bin start).
add_anomalies <- function(daily, climatology, variables) {
  yday <- as.integer(format(daily$date, "%j"))
  for (v in variables) {
    anom <- rep(NA_real_, nrow(daily))
    for (s in unique(daily$site)) {
      clim <- climatology$monthly |> filter(site == s, variable == v)
      f <- fit_periodic_climatology(clim$month, clim$monthly_mean)
      rows <- which(daily$site == s)
      if (!is.null(f)) anom[rows] <- daily[[v]][rows] - f(yday[rows])
    }
    daily[[paste0(v, "_anom")]] <- anom
  }
  daily
}

# =============================================================================
# 4. Bin length from the block-model residual ACF
# =============================================================================

# Fits the block model to daily presence minutes at one site and returns the
# residual ACF and the first lag >= 1 inside the 95% white-noise bounds.
# Change: the original returned which(acf < CI)[1], an INDEX into a vector
# whose first element is lag 0, i.e. lag + 1 (off by one), and tested only the
# upper bound; here the lag itself is returned and |ACF| is tested.
block_model_acf <- function(daily_site, species) {
  daily_site <- daily_site |> arrange(date)
  form <- as.formula(paste(species, "~", ACF_BLOCK_RHS))
  mod <- gam(form, data = daily_site, family = gam_family(), method = GAM_METHOD)
  r <- residuals(mod, type = "deviance")
  a <- acf(r, lag.max = ACF_LAG_MAX, plot = FALSE)
  bound <- qnorm((1 + ACF_CI_LEVEL) / 2) / sqrt(length(r))
  acf_tbl <- tibble(lag = as.numeric(a$lag), acf = as.numeric(a$acf), ci_bound = bound)
  inside <- acf_tbl |> filter(lag >= 1, abs(acf) < bound)
  list(values = acf_tbl,
       first_lag_inside = if (nrow(inside) > 0) min(inside$lag) else NA_real_)
}

# =============================================================================
# 5. Binning and ice filter
# =============================================================================

# Consecutive bin_days-day blocks counted from the first recording day of the
# site. Response = mean daily presence minutes in the bin; predictors = bin
# means; detection_days counts days with presence.
# Change: the original used lubridate::floor_date(date, "N days"), which
# restarts the blocks on the 1st of every month (short bins at month ends)
# rather than counting consecutive blocks from recording start.
# Axial circular mean of angles in degrees with period 180 (theta == theta +
# 180): double the angles, average the unit vectors, halve the mean direction.
# Result in (-90, 90].
axial_mean_deg <- function(theta) {
  rad <- 2 * theta * pi / 180
  m <- atan2(mean(sin(rad)), mean(cos(rad))) * 180 / pi / 2
  if_else(m <= -90, m + 180, m)
}

# Cyclic (axial) predictors are averaged with the axial circular mean, so that
# e.g. -89 and 89 deg average to +-90 deg, not 0 deg.
bin_site_data <- function(daily_site, species, bin_days, predictors) {
  first_day <- min(daily_site$date)
  linear <- setdiff(predictors, names(CYCLIC_VARS))
  axial  <- intersect(predictors, names(CYCLIC_VARS))
  daily_site |>
    mutate(day_index = as.integer(date - first_day),
           bin_index = day_index %/% bin_days) |>
    group_by(site, bin_index) |>
    summarise(bin_start = first_day + first(bin_index) * bin_days,
              n_days = n(),
              detection_days = sum(.data[[species]] > 0),
              rec_minutes = mean(rec_minutes),
              across(all_of(c(species, linear)), mean),
              across(all_of(axial), axial_mean_deg),
              .groups = "drop") |>
    filter(n_days >= BIN_MIN_DAYS) |>
    relocate(site, bin_index, bin_start)
}

# Per-species ice cutoff: the highest bin-mean ice concentration at which the
# species was detected, pooled across sites. Bins above it are outside the
# conditions in which the species occurred and are removed.
# Change: replaces the per-site 95th percentile with a strict "<" filter,
# which dropped every bin when the threshold was 0 (BW29 at CI).
ice_cutoff <- function(binned_all, species) {
  pres <- binned_all |> filter(.data[[species]] > 0)
  if (nrow(pres) == 0) return(NA_real_)
  max(pres$ice_conc)
}

# =============================================================================
# 6. Predictor screening
# =============================================================================

remove_zero_var <- function(df, vars) {
  vars[vapply(df[vars], function(x) isTRUE(sd(x, na.rm = TRUE) > 0), logical(1))]
}

# Sequential VIF screening (kept as coded): fit a gaussian GLM of the response
# on all predictors, compute car::vif (GVIF^(1/(2*Df)) when reported), drop the
# predictor with the largest value until all are <= threshold. If the GLM or
# VIF fails (e.g., aliased coefficients when predictors outnumber bins), drop
# one member of the most correlated pair: the one with the larger mean |r| to
# all remaining predictors.
vif_stepwise_drop <- function(data, response, predictors, vif_thresh = VIF_THRESH,
                              max_iter = 200) {
  preds <- intersect(predictors, names(data))
  steps <- list()

  get_vif <- function(fit) {
    v <- tryCatch(car::vif(fit), error = function(e) NA)
    if (is.matrix(v)) {
      col <- if ("GVIF^(1/(2*Df))" %in% colnames(v)) "GVIF^(1/(2*Df))" else 1
      v <- setNames(as.numeric(v[, col]), rownames(v))
    }
    v
  }

  fallback_drop <- function(preds) {
    x <- data[, preds, drop = FALSE]
    cm <- suppressWarnings(cor(x, use = "pairwise.complete.obs"))
    upper <- cm
    upper[lower.tri(upper, diag = TRUE)] <- NA
    pair <- which(abs(upper) == max(abs(upper), na.rm = TRUE), arr.ind = TRUE)[1, ]
    a <- preds[pair[1]]
    b <- preds[pair[2]]
    if (mean(abs(cm[a, ]), na.rm = TRUE) >= mean(abs(cm[b, ]), na.rm = TRUE)) a else b
  }

  iter <- 0
  while (iter < max_iter && length(preds) >= 2) {
    iter <- iter + 1
    fml <- as.formula(paste(response, "~", paste(preds, collapse = " + ")))
    fit <- tryCatch(glm(fml, data = data, family = gaussian()), error = function(e) NULL)
    v <- if (is.null(fit)) NA else get_vif(fit)

    if (anyNA(v)) {
      drop_var <- fallback_drop(preds)
      steps[[iter]] <- tibble(step = iter, n_predictors = length(preds),
                              action = "drop (VIF not estimable: correlation fallback)",
                              variable = drop_var, max_vif = NA_real_)
      preds <- setdiff(preds, drop_var)
      next
    }
    worst <- names(v)[which.max(v)]
    if (max(v) <= vif_thresh) {
      steps[[iter]] <- tibble(step = iter, n_predictors = length(preds),
                              action = "stop (all VIF <= threshold)",
                              variable = worst, max_vif = max(v))
      break
    }
    steps[[iter]] <- tibble(step = iter, n_predictors = length(preds),
                            action = "drop (max VIF)", variable = worst, max_vif = max(v))
    preds <- setdiff(preds, worst)
  }
  list(kept = preds, steps = bind_rows(steps))
}

# =============================================================================
# 7. Model fitting
# =============================================================================

# Tweedie GAM with s(x, k = GAM_K) for each term. The formula environment is
# set to the global environment so saved .rds files do not carry copies of
# function frames.
# Cyclic predictors (CYCLIC_VARS, e.g. fsle_orient) get a cyclic cubic
# regression spline whose ends join at the given knots (-90 == 90 deg).
fit_tw_gam <- function(data, response, terms, select = TRUE) {
  is_cyc <- terms %in% names(CYCLIC_VARS)
  smooths <- if_else(is_cyc,
                     paste0("s(", terms, ", bs = \"cc\", k = ", GAM_K, ")"),
                     paste0("s(", terms, ", k = ", GAM_K, ")"))
  rhs <- if (length(terms) == 0) "1" else paste(smooths, collapse = " + ")
  fml <- as.formula(paste(response, "~", rhs), env = globalenv())
  knots <- if (any(is_cyc)) CYCLIC_VARS[terms[is_cyc]] else NULL
  gam(fml, data = data, family = gam_family(), method = GAM_METHOD, select = select,
      knots = knots)
}

term_variable <- function(term_label) gsub("^s\\(|\\)$", "", term_label)

# Single-term screening of every candidate (kept as coded: y ~ s(x, k = 4),
# REML, no select). Results feed theme selection and the Stage 2 table; the
# fitted curve of each single-term model is summarised as a shape class.
single_term_screen <- function(data, response, candidates) {
  fits <- map(candidates, function(v) {
    mod <- tryCatch(fit_tw_gam(data, response, v, select = FALSE),
                    error = function(e) NULL)
    if (is.null(mod)) {
      return(list(result = tibble(variable = v, fit_ok = FALSE), grid = NULL))
    }
    st <- summary(mod)$s.table
    grid <- smooth_grid(mod, v)
    shape <- classify_smooth_shape(grid$.estimate, st[1, "p-value"], grid$x,
                                   cyclic = v %in% names(CYCLIC_VARS))
    list(result = tibble(variable = v, fit_ok = TRUE, AIC = AIC(mod),
                         dev_expl = summary(mod)$dev.expl, edf = st[1, "edf"],
                         F = st[1, "F"], p_value = st[1, "p-value"],
                         shape = shape, direction = shape_direction(shape)),
         grid = mutate(grid, variable = v))
  })
  list(results = bind_rows(map(fits, "result")), grids = bind_rows(map(fits, "grid")))
}

# One winner per theme by lowest single-term AIC among the eligible
# candidates of that theme (`eligible`: usable candidates; VIF is applied to
# the winners afterwards, see 01_stage1_automated_selection.R). The original
# prefer_anoms step is no longer needed: each variable enters either raw or
# as anomaly.
select_theme_winners <- function(screen, eligible) {
  screen |>
    left_join(parse_variable(screen$variable), by = "variable") |>
    mutate(competes = variable %in% eligible & fit_ok) |>
    group_by(theme) |>
    mutate(aic_rank_all = rank(AIC, ties.method = "first", na.last = "keep"),
           aic_rank_theme = if_else(competes,
                                    rank(if_else(competes, AIC, Inf), ties.method = "first"),
                                    NA_real_),
           theme_winner = competes & aic_rank_theme == 1) |>
    ungroup()
}

# Maximum number of smooth terms in a habitat model: floor(n_bins /
# MAX_TERMS_PER_BINS), and at most LOW_SAMPLE_MAX_TERMS in low-sample models.
term_cap <- function(n_bins, low_sample) {
  min(floor(n_bins / MAX_TERMS_PER_BINS), if (low_sample) LOW_SAMPLE_MAX_TERMS else Inf)
}

# Backward elimination (kept as coded): fit all candidate terms with
# select = TRUE, drop the term with the largest p-value until every p <=
# P_THRESH; the last remaining term is always kept. Added: the term with the
# largest p is also dropped while more than max_terms (term_cap()) remain;
# if the last remaining term has p > P_THRESH the model is flagged
# no_significant_terms (fit kept for diagnostics, reported without terms).
# Used by Stage 1 (theme winners) and Stage 3 (harmonized term sets).
backward_eliminate <- function(data, response, predictors, max_terms = Inf) {
  preds <- predictors
  steps <- list()
  repeat {
    mod <- fit_tw_gam(data, response, preds, select = TRUE)
    pv <- setNames(summary(mod)$s.pv, term_variable(rownames(summary(mod)$s.table)))
    worst <- names(pv)[which.max(pv)]
    reason <- case_when(
      length(preds) == 1               ~ "stop (one term left)",
      max(pv) > P_THRESH               ~ "drop (p > threshold)",
      length(preds) > max_terms        ~ "drop (term cap)",
      TRUE                             ~ "stop (all p <= threshold)")
    steps[[length(steps) + 1]] <- tibble(
      step = length(steps) + 1, formula = deparse1(formula(mod)),
      n_terms = length(preds), max_p = max(pv), worst_term = worst, action = reason)
    if (startsWith(reason, "stop")) break
    preds <- setdiff(preds, worst)
  }
  list(model = mod, terms = preds, steps = bind_rows(steps),
       no_significant_terms = max(pv) > P_THRESH)
}

# =============================================================================
# 8. Smooth shapes
# =============================================================================

# Centered partial effect (log scale) of one smooth on an even grid over the
# observed predictor range, with 95% pointwise CI.
# Cyclic predictors are evaluated over the full circle (CYCLIC_VARS knots).
smooth_grid <- function(mod, variable, n = SHAPE_GRID_N) {
  # gratia needs every model covariate in `data`; the others are set to their
  # medians (they do not affect this term's centered smooth).
  grid_data <- if (variable %in% names(CYCLIC_VARS)) {
    g <- tibble(!!variable := seq(CYCLIC_VARS[[variable]][1], CYCLIC_VARS[[variable]][2],
                                  length.out = n))
    for (v in setdiff(names(mod$model)[-1], variable)) g[[v]] <- median(mod$model[[v]])
    g
  }
  smooth_estimates(mod, select = paste0("s(", variable, ")"), n = n, data = grid_data) |>
    add_confint() |>
    rename(x = all_of(variable)) |>
    select(.smooth, x, .estimate, .se, .lower_ci, .upper_ci)
}

# Objective shape class of a centered smooth evaluated on a grid (decision
# 2.6: p-value based, the CI band is not used):
#   not significant - term p > P_THRESH
#   increasing - falls account for <= SHAPE_MONO_TOL of total variation
#   decreasing - rises account for <= SHAPE_MONO_TOL of total variation
#   otherwise the sequence of rising/falling runs is reduced to runs that
#   change the curve by >= SHAPE_MONO_TOL of total variation (slopes smaller
#   than SHAPE_SLOPE_NEGLIG x the steepest slope are ignored):
#   peaked (rise then fall), U-shaped (fall then rise), complex (more turns).
#   A significant smooth penalised to an exactly constant curve is "flat".
# Cyclic predictors (CYCLIC_VARS; x = grid over the full circle) have no
# increasing/decreasing: the curve is treated as a closed loop, the same run
# reduction counts its major peaks, and the shape is reported as
#   "cyclic, unimodal (peak P deg, trough T deg)" or "cyclic, multimodal (...)"
# with P/T = orientations of the global maximum/minimum of the fitted curve.
classify_smooth_shape <- function(estimate, p_value, x = NULL, cyclic = FALSE) {
  if (is.na(p_value) || p_value > P_THRESH) return("not significant")
  if (cyclic) {
    loop <- estimate[-length(estimate)]        # last grid point == first (-90 == 90)
    d <- diff(c(loop, loop[1]))
  } else {
    d <- diff(estimate)
  }
  total <- sum(abs(d))
  if (total == 0) return("flat")
  if (cyclic) {
    keep <- abs(d) >= SHAPE_SLOPE_NEGLIG * max(abs(d))
    runs <- rle(sign(d[keep]))
    run_change <- as.numeric(tapply(d[keep], rep(seq_along(runs$lengths), runs$lengths), sum))
    signs <- rle(sign(run_change[abs(run_change) >= SHAPE_MONO_TOL * total]))$values
    if (length(signs) > 1 && signs[1] == signs[length(signs)]) signs <- signs[-1]
    n_peaks <- sum(signs > 0)
    return(sprintf("cyclic, %s (peak %d°, trough %d°)",
                   if (n_peaks <= 1) "unimodal" else "multimodal",
                   round(x[which.max(estimate)]), round(x[which.min(estimate)])))
  }
  rise <- sum(d[d > 0])
  fall <- -sum(d[d < 0])
  if (fall <= SHAPE_MONO_TOL * total) return("increasing")
  if (rise <= SHAPE_MONO_TOL * total) return("decreasing")
  keep <- abs(d) >= SHAPE_SLOPE_NEGLIG * max(abs(d))
  runs <- rle(sign(d[keep]))
  run_id <- rep(seq_along(runs$lengths), runs$lengths)
  run_change <- as.numeric(tapply(d[keep], run_id, sum))
  signs <- rle(sign(run_change[abs(run_change) >= SHAPE_MONO_TOL * total]))$values
  if (identical(as.numeric(signs), c(1, -1))) return("peaked")
  if (identical(as.numeric(signs), c(-1, 1))) return("U-shaped")
  if (length(signs) == 1) return(if (signs > 0) "increasing" else "decreasing")
  "complex"
}

shape_direction <- function(shape) {
  case_when(startsWith(shape, "cyclic") ~ "cyclic",
            shape == "increasing" ~ "positive",
            shape == "decreasing" ~ "negative",
            shape %in% c("flat", "not significant") ~ "none",
            TRUE                  ~ "mixed")
}

# =============================================================================
# 9. Documentation of a final model (shared by Stage 1 and Stage 3)
# =============================================================================

# Per-term deviance explained: full-model deviance explained minus that of the
# model refitted without the term (same family/settings). Shared deviance is
# attributed to no term, so per-term values need not sum to the model total.
drop_one_deviance <- function(mod, data, response, terms) {
  full <- summary(mod)$dev.expl
  map_dbl(terms, function(v) {
    reduced <- fit_tw_gam(data, response, setdiff(terms, v), select = TRUE)
    full - summary(reduced)$dev.expl
  }) |> setNames(terms)
}

# Diagnostics saved for every final model: gam.check (text + 4-panel plot),
# k-index, concurvity(full = TRUE), gratia::appraise, residual ACF on the
# regular bin grid (ice-removed bins left as gaps so lags stay in bin units).
model_diagnostics <- function(mod, data, id, diag_dir) {
  set.seed(RANDOM_SEED)
  png(file.path(diag_dir, paste0(id, "_gam_check.png")), width = 2400, height = 2400, res = 300)
  check_txt <- capture.output(gam.check(mod))
  dev.off()
  writeLines(check_txt, file.path(diag_dir, paste0(id, "_gam_check.txt")))

  set.seed(RANDOM_SEED)
  kc <- k.check(mod)
  conc <- concurvity(mod, full = TRUE)
  conc_tbl <- as_tibble(t(conc), rownames = "term") |>
    filter(term != "para")

  appraise_plot <- tryCatch(appraise(mod) + plot_annotation(title = paste(id, "- appraise")),
                            error = function(e) NULL)
  if (!is.null(appraise_plot)) {
    ggsave(file.path(diag_dir, paste0(id, "_appraise.png")), appraise_plot,
           width = 10, height = 8, dpi = 300, bg = "white")
  }

  r_grid <- rep(NA_real_, max(data$bin_index) + 1)
  r_grid[data$bin_index + 1] <- residuals(mod, type = "deviance")
  n_obs <- sum(!is.na(r_grid))
  a <- acf(r_grid, lag.max = min(RESID_ACF_MAX_LAG, n_obs - 1), na.action = na.pass,
           plot = FALSE)
  bound <- qnorm((1 + ACF_CI_LEVEL) / 2) / sqrt(n_obs)
  acf_tbl <- tibble(lag = as.numeric(a$lag), acf = as.numeric(a$acf), ci_bound = bound) |>
    mutate(outside = lag >= 1 & abs(acf) > ci_bound)
  acf_plot <- ggplot(acf_tbl, aes(lag, acf)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(-bound, bound), linetype = "dashed", colour = "#2a78d6") +
    geom_segment(aes(xend = lag, yend = 0), linewidth = 0.8) +
    labs(title = paste(id, "- deviance residual ACF"),
         x = "Lag (bins)", y = "ACF") +
    theme_bw()
  ggsave(file.path(diag_dir, paste0(id, "_resid_acf.png")), acf_plot,
         width = 6, height = 4, dpi = 300, bg = "white")

  # Ljung-Box test over lags 1..h, computed from the gap-aware ACF above:
  # Q = n (n + 2) sum_k r_k^2 / (n - k), chi-squared with h df (n = bins
  # with residuals). stats::Box.test() would ignore the gaps.
  h  <- min(LJUNG_BOX_LAGS, n_obs - 2)
  rk <- acf_tbl$acf[match(seq_len(h), acf_tbl$lag)]
  lb_q <- n_obs * (n_obs + 2) * sum(rk^2 / (n_obs - seq_len(h)))
  lb_p <- pchisq(lb_q, df = h, lower.tail = FALSE)

  list(
    summary = tibble(
      k_index_min = min(kc[, "k-index"]),
      k_pvalue_min = min(kc[, "p-value"]),
      concurvity_max = if (nrow(conc_tbl) > 0) max(conc_tbl$worst) else NA_real_,
      resid_acf_lag1 = acf_tbl$acf[acf_tbl$lag == 1],
      resid_acf_lags_outside = paste(acf_tbl$lag[acf_tbl$outside], collapse = ";"),
      ljung_box_lags = h,
      ljung_box_q = lb_q,
      ljung_box_p = lb_p,
      resid_acf_ok = lb_p > ACF_TEST_ALPHA),
    acf = acf_tbl,
    concurvity = conc_tbl,
    k_check = as_tibble(kc, rownames = "term")
  )
}

# Fit-independent documentation of one final model: summary row, term table
# (edf, p, per-term deviance, shape of the final smooth), diagnostics, and the
# .rds file read by 08_results_gams/. Used identically by Stage 1 and Stage 3
# so the 08 scripts can read either stage.
document_final_model <- function(mod, data, species, site, stage, meta, out_dir) {
  id <- model_id(species, site)
  diag_dir <- file.path(out_dir, "diagnostics")
  model_dir <- file.path(out_dir, "models")
  for (d in c(diag_dir, model_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # A Stage 1 model whose last remaining term is not significant is reported
  # as "no significant terms": diagnostics are kept, the term table is empty.
  no_sig <- isTRUE(meta$no_significant_terms)
  summ <- summary(mod)
  diag <- model_diagnostics(mod, data, id, diag_dir)

  term_tbl <- if (no_sig) {
    tibble(species = character(), site = character(), stage = character(),
           term = character(), variable = character(), theme = character(),
           edf = numeric(), ref_df = numeric(), F = numeric(), p_value = numeric(),
           dev_expl_term = numeric(), shape = character(), direction = character())
  } else {
    terms <- term_variable(rownames(summ$s.table))
    dev_terms <- drop_one_deviance(mod, data, species, terms)
    shapes <- map2(terms, summ$s.pv, function(v, p) {
      g <- smooth_grid(mod, v)
      s <- classify_smooth_shape(g$.estimate, p, g$x, cyclic = v %in% names(CYCLIC_VARS))
      tibble(variable = v, shape = s, direction = shape_direction(s))
    }) |> bind_rows()
    tibble(
      species = species, site = site, stage = stage, term = rownames(summ$s.table),
      variable = terms, theme = parse_variable(terms)$theme,
      edf = summ$s.table[, "edf"], ref_df = summ$s.table[, "Ref.df"],
      F = summ$s.table[, "F"], p_value = summ$s.table[, "p-value"],
      dev_expl_term = unname(dev_terms[terms])) |>
      left_join(shapes, by = "variable")
  }

  summary_row <- tibble(
    species = species, site = site, stage = stage,
    formula = deparse1(formula(mod)), n_terms = nrow(term_tbl),
    no_significant_terms = no_sig,
    term_cap = if (is.null(meta$term_cap)) NA_real_ else meta$term_cap,
    n_bins = nrow(data), n_detection_bins = sum(data[[species]] > 0),
    low_sample = meta$low_sample, bin_days = meta$bin_days, ice_cutoff = meta$ice_cutoff,
    dev_expl = summ$dev.expl, adj_r2 = summ$r.sq,
    tweedie_p = mod$family$getTheta(TRUE), scale = summ$scale,
    AIC = AIC(mod), REML = mod$gcv.ubre,
    all_terms_significant = all(summ$s.pv <= P_THRESH)) |>
    bind_cols(diag$summary)

  saveRDS(list(model = mod, data = data, species = species, site = site, stage = stage,
               meta = meta, summary = summary_row, terms = term_tbl),
          file.path(model_dir, paste0(id, ".rds")))

  list(summary = summary_row, terms = term_tbl,
       acf = mutate(diag$acf, species = species, site = site, .before = 1),
       concurvity = mutate(diag$concurvity, species = species, site = site, .before = 1),
       k_check = mutate(diag$k_check, species = species, site = site, .before = 1))
}

# Write the per-model documentation of all models of one stage as CSV tables
# (same file names for Stage 1 and Stage 3).
write_model_tables <- function(docs, table_dir) {
  write_csv(bind_rows(map(docs, "summary")),    file.path(table_dir, "final_model_summary.csv"))
  write_csv(bind_rows(map(docs, "terms")),      file.path(table_dir, "final_model_terms.csv"))
  write_csv(bind_rows(map(docs, "acf")),        file.path(table_dir, "final_model_resid_acf.csv"))
  write_csv(bind_rows(map(docs, "concurvity")), file.path(table_dir, "final_model_concurvity.csv"))
  write_csv(bind_rows(map(docs, "k_check")),    file.path(table_dir, "final_model_k_check.csv"))
}

# Read every final model (.rds) of one stage; used by the 08 scripts.
read_stage_models <- function(stage) {
  files <- list.files(file.path(STATS_DIR, stage, "models"), pattern = "\\.rds$",
                      full.names = TRUE)
  if (length(files) == 0) stop("No models found for stage '", stage, "' in ", STATS_DIR)
  setNames(map(files, readRDS), sub("\\.rds$", "", basename(files)))
}

# Standardized anomalies (used only when CHL_STANDARDIZE is TRUE): each
# <var>_anom is divided by the SD of that anomaly in its calendar month over
# the 10-year record (env_daily), interpolated to daily values with the same
# periodic spline as the climatology. The anomaly variance of chl-a is
# strongly seasonal (e.g. KGI SD 0.04 in September vs 0.21 in May), so raw
# anomalies give large weight to a few months.
standardize_anomalies <- function(daily, env_daily, variables) {
  for (v in variables) {
    col <- paste0(v, "_anom")
    for (s in unique(daily$site)) {
      e <- env_daily[env_daily$site == s, ]
      sd_month <- tapply(e[[col]], as.integer(format(e$date, "%m")), sd, na.rm = TRUE)
      f <- fit_periodic_climatology(as.integer(names(sd_month)), as.numeric(sd_month))
      rows <- which(daily$site == s)
      daily[[col]][rows] <- daily[[col]][rows] / f(as.integer(format(daily$date[rows], "%j")))
    }
  }
  daily
}

# VIF of a set of terms (linear Gaussian fit, as in vif_stepwise_drop)
vif_of <- function(data, response, terms) {
  if (length(terms) < 2) return(setNames(rep(1, length(terms)), terms))
  fit <- glm(reformulate(terms, response), data = data, family = gaussian())
  v <- car::vif(fit)
  if (is.matrix(v)) v <- setNames(as.numeric(v[, ncol(v)]), rownames(v))
  v
}

# Selection with a collinearity CHECK on the final model (decision
# N. Posdaljian 2026-11-20; replaces VIF screening of the candidates, which
# could drop the better-supported of two collinear predictors on a near-tie):
# backward elimination with the term cap from the full candidate set; if the
# final model has a term with VIF > vif_thresh, the least significant of the
# terms above the threshold is removed from the candidates and selection is
# repeated, until every final term has VIF <= vif_thresh.
select_vif_checked <- function(data, response, predictors, max_terms = Inf,
                               vif_thresh = VIF_THRESH) {
  preds <- remove_zero_var(data, predictors)
  checks <- list()
  repeat {
    be <- backward_eliminate(data, response, preds, max_terms = max_terms)
    if (be$no_significant_terms || length(be$terms) < 2) break
    v <- vif_of(data, response, be$terms)
    pv <- setNames(summary(be$model)$s.pv,
                   term_variable(rownames(summary(be$model)$s.table)))
    checks[[length(checks) + 1]] <- tibble(check = length(checks) + 1,
                                           terms = paste(be$terms, collapse = " + "),
                                           max_vif = max(v), worst_vif = names(v)[which.max(v)])
    if (max(v) <= vif_thresh) break
    high <- names(v)[v > vif_thresh]
    drop <- high[which.max(pv[high])]
    checks[[length(checks)]]$dropped <- drop
    preds <- setdiff(preds, drop)
  }
  be$vif_checks <- bind_rows(checks)
  be
}
