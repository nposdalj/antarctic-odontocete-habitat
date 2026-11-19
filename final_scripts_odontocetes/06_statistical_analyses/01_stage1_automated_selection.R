# =============================================================================
# Script  : 01_stage1_automated_selection.R
# Section : 2.6 Statistical analyses (Stage 1 results reported in the supplement)
# Purpose : Stage 1 of the habitat-model workflow: automated predictor
#           selection for one Tweedie GAM per species x site. Steps:
#             1. 10-yr monthly climatology (periodic cubic spline); MLD,
#                temperature, salinity, oxygen (dive depths) and log10 3-month
#                lagged chl-a enter as anomalies; EKE (dive depths), FSLE and FSLV
#                orientation stay raw (decision 2026-11-19)
#             2. Bin length per species from block-model residual ACF
#                (max across sites, capped at BIN_CAP_DAYS)
#             3. Consecutive N-day bins from recording start
#             4. Per-species ice cutoff (max ice_conc at detection, pooled
#                across sites) -> drop bins above it
#             5. Inclusion rule (>= 10 bins with detections; 10-20 = low sample)
#             6. Zero-variance removal -> sequential VIF <= 5 -> one winner per
#                theme by single-term AIC -> backward elimination (p <= 0.05,
#                select = TRUE; at most floor(n_bins / 10) terms, 2 in
#                low-sample models)
#             7. Diagnostics and documentation of every final model
#           fsle_orient (axial, -90 == 90 deg) enters every fit as a cyclic
#           smooth s(fsle_orient, bs = "cc", k = 4), knots at -90/90 deg.
# Inputs  : F_MODEL_INPUT (model_input_daily.csv: recording days, daily env +
#           presence minutes per species + rec_minutes)
#           F_ENV_DAILY (env_daily_wide.csv: full climatology window)
#           [both defined in config.R; see 00_gam_parameters.R for overrides]
# Outputs : <OUTPUT_ROOT>/06_statistical_analyses/stage1/
#             tables/climatology_monthly.csv, climatology_coverage.csv,
#                    seasonal_variance.csv (site, variable, seasonal_var_frac,
#                    version_used)
#             tables/acf_block_model.csv, acf_block_values.csv, bin_lengths.csv
#             tables/ice_cutoffs.csv, model_inclusion.csv
#             tables/candidates.csv, vif_steps.csv, theme_selection.csv,
#                    single_term_smooths.csv, backward_elimination.csv
#             tables/final_model_summary.csv, final_model_terms.csv,
#                    final_model_resid_acf.csv, final_model_concurvity.csv,
#                    final_model_k_check.csv
#             binned_data/<species>_binned.csv (all sites, all candidates,
#                    ice_excluded flag; reused by Stage 3)
#             models/<species>_<site>.rds, diagnostics/<species>_<site>_*
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Consistent with the methods described for baleen whale targets in
#           [COMPANION PAPER CITATION]; differences for odontocetes are
#           species-specific depths, the ice cutoff (sea ice is modeled
#           separately in 04_sea_ice_models.R, not as a habitat candidate),
#           lags 1-6, the inclusion rule and the term cap. No effort offset
#           (decision 2.6: duty-cycle loss is far less than proportional to
#           on-time).
#           Run from the repository root:
#             Rscript final_scripts_odontocetes/06_statistical_analyses/01_stage1_automated_selection.R
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

out_dir    <- stage_dir("stage1")
table_dir  <- stage_dir("stage1", "tables")
binned_dir <- stage_dir("stage1", "binned_data")

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# model_input: one row per site x recording day (response + predictors).
# env_daily: one row per site x day over the climatology window (baselines).
model_input <- read_csv(GAM_MODEL_INPUT, show_col_types = FALSE) |>
  mutate(date = as.Date(date)) |>
  arrange(site, date)
env_daily <- read_csv(GAM_ENV_DAILY, show_col_types = FALSE) |>
  mutate(date = as.Date(date))

# Lagged chl-a on the log10 scale (CHL_LOG10), applied before the
# climatology so that anomalies are departures in log10 chl-a.
if (CHL_LOG10) {
  log_chl <- function(df) {
    lag_cols <- intersect(LAG_COLUMNS, names(df))
    if (any(unlist(df[lag_cols]) <= 0, na.rm = TRUE)) {
      stop("Non-positive lagged chl-a values; cannot log10-transform.")
    }
    mutate(df, across(all_of(lag_cols), log10))
  }
  model_input <- log_chl(model_input)
  env_daily   <- log_chl(env_daily)
}

# -----------------------------------------------------------------------------
# 2. Deseasoning: climatology from the full window, anomalies on recording days
# -----------------------------------------------------------------------------
# Every anomaly source needed by any modeled species, computed once per site.
anom_vars <- unique(unlist(map(MODELED_SPECIES, anomaly_sources)))
missing_anom <- setdiff(anom_vars, names(env_daily))
if (length(missing_anom) > 0) {
  warning("Anomaly sources missing from env_daily (skipped): ",
          paste(missing_anom, collapse = ", "), call. = FALSE)
}
anom_vars <- intersect(anom_vars, intersect(names(env_daily), names(model_input)))

climatology <- build_climatology(env_daily, anom_vars)
write_csv(climatology$monthly,  file.path(table_dir, "climatology_monthly.csv"))
write_csv(climatology$coverage, file.path(table_dir, "climatology_coverage.csv"))
write_csv(climatology$seasonal, file.path(table_dir, "seasonal_variance.csv"))

model_input <- add_anomalies(model_input, climatology, anom_vars) |>
  group_by(site) |>
  mutate(day_index = as.integer(date - min(date))) |>   # block-model time axis
  ungroup()

# -----------------------------------------------------------------------------
# 3. Species loop
# -----------------------------------------------------------------------------
# Results are collected in lists and written as tables at the end.
acf_summary <- list(); acf_values <- list(); bin_lengths <- list()
ice_tbl <- list(); inclusion <- list(); candidates_tbl <- list()
vif_tbl <- list(); theme_tbl <- list(); smooth_tbl <- list()
backward_tbl <- list(); docs <- list()

for (sp in MODELED_SPECIES) {
  message("\n==== ", sp, " (", species_label(sp), ") ====")

  # --- 3a. Residual-ACF bin length per site -----------------------------------
  # Only sites with >= 1 detection day are informative (no residual structure
  # in an all-zero series; Pm at KGI).
  for (st in SITE_ORDER) {
    daily_site <- model_input |> filter(site == st)
    n_det_days <- sum(daily_site[[sp]] > 0)
    res <- if (n_det_days > 0) {
      tryCatch(block_model_acf(daily_site, sp), error = function(e) NULL)
    }
    acf_summary[[length(acf_summary) + 1]] <- tibble(
      species = sp, site = st, n_days = nrow(daily_site), n_detection_days = n_det_days,
      block_model = paste(sp, "~", ACF_BLOCK_RHS),
      acf_lag1 = if (!is.null(res)) res$values$acf[res$values$lag == 1] else NA_real_,
      ci_bound = if (!is.null(res)) res$values$ci_bound[1] else NA_real_,
      first_lag_inside_ci = if (!is.null(res)) res$first_lag_inside else NA_real_,
      status = if (n_det_days == 0) "no detections"
               else if (is.null(res)) "block model failed"
               else if (is.na(res$first_lag_inside)) paste0("no lag <= ", ACF_LAG_MAX,
                                                            " inside CI")
               else "ok")
    if (!is.null(res)) {
      acf_values[[length(acf_values) + 1]] <- mutate(res$values, species = sp, site = st,
                                                    .before = 1)
    }
  }

  # One bin length per species: maximum across sites, capped. A site whose
  # ACF never enters the CI within ACF_LAG_MAX counts as exceeding the cap.
  sp_acf <- bind_rows(acf_summary) |> filter(species == sp, n_detection_days > 0)
  site_lags <- if_else(is.na(sp_acf$first_lag_inside_ci), Inf, sp_acf$first_lag_inside_ci)
  max_lag <- if (length(site_lags) > 0) max(site_lags) else NA_real_
  bin_days <- as.integer(min(max_lag, BIN_CAP_DAYS))
  bin_lengths[[sp]] <- tibble(species = sp, max_site_lag = max_lag, cap_days = BIN_CAP_DAYS,
                              bin_days = bin_days, capped = max_lag > BIN_CAP_DAYS)
  message("Bin length: ", bin_days, " d (max site lag ", max_lag, ")")

  # --- 3b. Binning ----------------------------------------------------------
  # Raw values AND all anomalies of this species are binned so Stage 3 can use
  # any of them without re-running Stage 1.
  sp_candidates <- c(raw_candidates(sp), paste0(intersect(anomaly_sources(sp), anom_vars),
                                                "_anom"))
  missing_cand <- setdiff(sp_candidates, names(model_input))
  if (length(missing_cand) > 0) {
    warning(sp, ": candidate columns missing from model input (skipped): ",
            paste(missing_cand, collapse = ", "), call. = FALSE)
  }
  sp_candidates <- intersect(sp_candidates, names(model_input))

  binned <- model_input |>
    group_split(site) |>
    map(\(d) bin_site_data(d, sp, bin_days, c(sp_candidates, ICE_MODEL_VAR))) |>
    bind_rows()

  # --- 3c. Ice filter ---------------------------------------------------------
  cutoff <- ice_cutoff(binned, sp)
  binned <- binned |> mutate(ice_excluded = ice_conc > cutoff)
  write_csv(binned, file.path(binned_dir, paste0(sp, "_binned.csv")))
  ice_tbl[[sp]] <- binned |>
    group_by(site) |>
    summarise(n_bins_before = n(), n_bins_removed = sum(ice_excluded),
              n_bins_after = sum(!ice_excluded), .groups = "drop") |>
    mutate(species = sp, ice_cutoff = cutoff, .before = 1)
  message("Ice cutoff: ", round(cutoff, 3))

  # --- 3d. Per-site selection -----------------------------------------------
  for (st in SITE_ORDER) {
    id <- model_id(sp, st)
    dat <- binned |> filter(site == st, !ice_excluded)
    n_det <- sum(dat[[sp]] > 0)
    status <- case_when(n_det < MIN_DETECTION_BINS   ~ "descriptive only",
                        n_det <= LOW_SAMPLE_MAX_BINS ~ "modeled (low sample)",
                        TRUE                         ~ "modeled")
    inclusion[[id]] <- tibble(species = sp, site = st, bin_days = bin_days,
                              n_bins = nrow(dat), n_detection_bins = n_det, status = status)
    message(st, ": ", nrow(dat), " bins, ", n_det, " with detections -> ", status)
    if (status == "descriptive only") next

    # Candidates with missing values or no variation at this site cannot be
    # used.
    # Each deseasoning-eligible variable enters raw OR as anomaly at this site
    # (seasonal_variance.csv), never both.
    cand <- intersect(site_candidates(sp, st, climatology$seasonal), sp_candidates)
    has_na <- cand[vapply(dat[cand], anyNA, logical(1))]
    usable <- remove_zero_var(dat, setdiff(cand, has_na))

    # Order (decision 2026-11-19): single-term screening of every usable
    # candidate -> one winner per theme by single-term AIC -> collinearity
    # (VIF) screening among the theme winners only -> backward elimination.
    # Screening the winners (not all candidates) means a strong predictor is
    # never removed only for being collinear with a weaker one of another
    # theme; two layers both survive only if they carry different information.
    screen <- single_term_screen(dat, sp, usable)
    themes <- select_theme_winners(screen$results, usable)
    theme_winners <- themes |> filter(theme_winner) |> pull(variable)
    vif <- vif_stepwise_drop(dat, sp, theme_winners)
    winners <- vif$kept
    theme_tbl[[id]] <- mutate(themes, species = sp, site = st, .before = 1,
                              vif_retained = variable %in% winners)
    smooth_tbl[[id]] <- mutate(screen$grids, species = sp, site = st, .before = 1)
    vif_tbl[[id]] <- mutate(vif$steps, species = sp, site = st, .before = 1)
    candidates_tbl[[id]] <- parse_variable(cand) |>
      mutate(species = sp, site = st, .before = 1,
             screening = case_when(variable %in% has_na ~ "removed (missing values)",
                                   !variable %in% usable ~ "removed (zero variance)",
                                   !variable %in% theme_winners ~ "not theme winner",
                                   variable %in% winners ~ "VIF retained",
                                   TRUE ~ "removed (VIF)"))
    message("Theme winners: ", paste(theme_winners, collapse = ", "),
            "; after VIF: ", paste(winners, collapse = ", "))

    # Backward elimination with the term cap (floor(n_bins / 10); 2 terms in
    # low-sample models).
    low_sample <- status == "modeled (low sample)"
    be <- backward_eliminate(dat, sp, winners,
                             max_terms = term_cap(nrow(dat), low_sample))
    backward_tbl[[id]] <- mutate(be$steps, species = sp, site = st, .before = 1)
    message("Final: ", deparse1(formula(be$model)))

    docs[[id]] <- document_final_model(
      be$model, dat, sp, st, stage = "stage1",
      meta = list(bin_days = bin_days, ice_cutoff = cutoff, low_sample = low_sample,
                  theme_winners = winners, term_cap = term_cap(nrow(dat), low_sample),
                  no_significant_terms = be$no_significant_terms),
      out_dir = out_dir)
  }
}

# -----------------------------------------------------------------------------
# 4. Save tables
# -----------------------------------------------------------------------------
write_csv(bind_rows(acf_summary),    file.path(table_dir, "acf_block_model.csv"))
write_csv(bind_rows(acf_values),     file.path(table_dir, "acf_block_values.csv"))
write_csv(bind_rows(bin_lengths),    file.path(table_dir, "bin_lengths.csv"))
write_csv(bind_rows(ice_tbl),        file.path(table_dir, "ice_cutoffs.csv"))
write_csv(bind_rows(inclusion),      file.path(table_dir, "model_inclusion.csv"))
write_csv(bind_rows(candidates_tbl), file.path(table_dir, "candidates.csv"))
write_csv(bind_rows(vif_tbl),        file.path(table_dir, "vif_steps.csv"))
write_csv(bind_rows(theme_tbl),      file.path(table_dir, "theme_selection.csv"))
write_csv(bind_rows(smooth_tbl),     file.path(table_dir, "single_term_smooths.csv"))
write_csv(bind_rows(backward_tbl),   file.path(table_dir, "backward_elimination.csv"))
write_model_tables(docs, table_dir)

message("\nStage 1 complete: ", length(docs), " models. Outputs in ", out_dir)
