# =============================================================================
# Script  : 06_permutation_test.R
# Section : 2.6 Statistical analyses (supplement) - how often does the model
#           selection find "significant" habitat terms by chance?
# Purpose : Time-shift permutation test of the whole habitat-model selection
#           (single-term screen -> VIF <= 5 -> backward elimination with the
#           term cap, as in Stage 1). For each species x site, the
#           environmental predictors on the recording days are replaced by the
#           same site's predictors from another year (shifted by a whole
#           number of years plus a random 0-60 day offset), and the selection
#           is re-run. Shifting keeps each predictor's autocorrelation, its
#           season and the correlations among predictors, but breaks any real
#           link with the whales. Presence and sea ice (and therefore the
#           ice filter and the bins) are not shifted.
#           The observed result is compared with the null distribution of
#           (a) the number of terms retained and (b) deviance explained.
# Inputs  : GAM_MODEL_INPUT, GAM_ENV_DAILY (10-year daily predictors)
#           <OUTPUT_ROOT>/06_statistical_analyses/stage1/tables/model_inclusion.csv
#           <OUTPUT_ROOT>/06_statistical_analyses/stage1/tables/ice_cutoffs.csv
# Outputs : outputs/06_statistical_analyses/permutation/permutation_runs.csv
#             one row per species x site x run (run 0 = observed data)
#           outputs/06_statistical_analyses/permutation/permutation_summary.csv
#             observed vs null: proportion of shifted runs retaining >= 1
#             term, and permutation p-values for deviance explained and
#             number of terms
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : N_PERM (default 200) via env ODONTOCETE_N_PERM. Runtime ~1-2 h.
#           Tests the Stage 1 selection from all candidates; the Stage 3
#           harmonized sets are a subset of the union of Stage 1 picks, so
#           this is the more complete test of the selection procedure.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

N_PERM     <- as.integer(Sys.getenv("ODONTOCETE_N_PERM", unset = "200"))
MAX_OFFSET <- 60   # days added to the whole-year shift
set.seed(20261120)

out_dir <- file.path(STATS_DIR, "permutation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
s1_tables <- file.path(STATS_DIR, "stage1", "tables")

# -----------------------------------------------------------------------------
# 1. Predictors as in Stage 1 (log10 chl-a, anomalies from the climatology)
# -----------------------------------------------------------------------------
model_input <- read_csv(GAM_MODEL_INPUT, show_col_types = FALSE) |>
  mutate(date = as.Date(date)) |> arrange(site, date)
env_daily <- read_csv(GAM_ENV_DAILY, show_col_types = FALSE) |>
  mutate(date = as.Date(date))
if (CHL_LOG10) {
  model_input <- mutate(model_input, across(any_of(LAG_COLUMNS), log10))
  env_daily   <- mutate(env_daily,   across(any_of(LAG_COLUMNS), log10))
}
anom_vars <- unique(unlist(map(MODELED_SPECIES, anomaly_sources)))
anom_vars <- intersect(anom_vars, intersect(names(env_daily), names(model_input)))
climatology <- build_climatology(env_daily, anom_vars)
model_input <- add_anomalies(model_input, climatology, anom_vars)
env_daily   <- add_anomalies(env_daily,   climatology, anom_vars)
if (CHL_STANDARDIZE) {   # model_input first, using the unstandardized env_daily
  model_input <- standardize_anomalies(model_input, env_daily, LAG_COLUMNS)
  env_daily   <- standardize_anomalies(env_daily, env_daily, LAG_COLUMNS)
}

inclusion <- read_csv(file.path(s1_tables, "model_inclusion.csv"), show_col_types = FALSE) |>
  filter(status != "descriptive only")
cutoffs <- read_csv(file.path(s1_tables, "ice_cutoffs.csv"), show_col_types = FALSE) |>
  distinct(species, ice_cutoff)

# -----------------------------------------------------------------------------
# 2. One selection run (Stage 1 procedure) on one site's daily data
# -----------------------------------------------------------------------------
run_selection <- function(daily_site, sp, st, bin_days, cutoff, low_sample) {
  cand <- intersect(site_candidates(sp, st, NULL), names(daily_site))
  dat <- bin_site_data(daily_site, sp, bin_days, c(cand, ICE_MODEL_VAR)) |>
    filter(ice_conc <= cutoff)
  usable <- remove_zero_var(dat, cand[!vapply(dat[cand], anyNA, logical(1))])
  screen <- single_term_screen(dat, sp, usable)
  winners <- select_theme_winners(screen$results, usable) |>
    filter(theme_winner) |> pull(variable)
  kept <- vif_stepwise_drop(dat, sp, winners)$kept
  be <- backward_eliminate(dat, sp, kept, max_terms = term_cap(nrow(dat), low_sample))
  terms <- if (be$no_significant_terms) character(0) else be$terms
  tibble(n_terms = length(terms),
         dev_expl = if (length(terms) == 0) 0 else 100 * summary(be$model)$dev.expl,
         terms = paste(terms, collapse = " + "))
}

# Replace the predictor columns of the recording days with the same site's
# values `shift_days` later (ice and presence stay on the real dates).
shift_predictors <- function(daily_site, env_site, shift_days, cols) {
  src <- env_site |> select(date, all_of(cols)) |>
    mutate(date = date - shift_days)          # value from date + shift
  daily_site |> select(-all_of(cols)) |> left_join(src, by = "date")
}

# Valid shifts: whole years (not 0) plus an offset, with the shifted window
# fully inside the 10-year record at that site.
draw_shifts <- function(dates, env_dates, n) {
  years <- setdiff(-10:10, 0)
  out <- integer(0)
  while (length(out) < n) {
    s <- round(365.25 * sample(years, 1)) + sample(-MAX_OFFSET:MAX_OFFSET, 1)
    if (all((dates + s) %in% env_dates)) out <- c(out, s)
  }
  out
}

# -----------------------------------------------------------------------------
# 3. Observed run and N_PERM shifted runs per species x site
# -----------------------------------------------------------------------------
runs <- list()
for (i in seq_len(nrow(inclusion))) {
  sp <- inclusion$species[i]; st <- inclusion$site[i]
  bin_days <- inclusion$bin_days[i]
  low_sample <- inclusion$status[i] == "modeled (low sample)"
  cutoff <- cutoffs$ice_cutoff[cutoffs$species == sp]
  daily_site <- filter(model_input, site == st)
  env_site <- filter(env_daily, site == st)
  cols <- setdiff(intersect(site_candidates(sp, st, NULL), names(daily_site)), ICE_MODEL_VAR)

  obs <- run_selection(daily_site, sp, st, bin_days, cutoff, low_sample)
  message(sp, " ", st, ": observed ", obs$n_terms, " terms (", obs$terms, "), ",
          round(obs$dev_expl, 1), "% deviance")
  shifts <- draw_shifts(daily_site$date, env_site$date, N_PERM)
  null <- map_dfr(seq_along(shifts), function(k) {
    r <- tryCatch(
      run_selection(shift_predictors(daily_site, env_site, shifts[k], cols),
                    sp, st, bin_days, cutoff, low_sample),
      error = function(e) tibble(n_terms = NA_integer_, dev_expl = NA_real_,
                                 terms = paste("error:", conditionMessage(e))))
    mutate(r, run = k, shift_days = shifts[k])
  })
  runs[[paste(sp, st)]] <- bind_rows(mutate(obs, run = 0L, shift_days = 0L), null) |>
    mutate(species = sp, site = st, .before = 1)
  write_csv(bind_rows(runs), file.path(out_dir, "permutation_runs.csv"))  # checkpoint
}

runs <- bind_rows(runs)
summary_tbl <- runs |>
  group_by(species, site) |>
  summarise(observed_terms = terms[run == 0],
            observed_n_terms = n_terms[run == 0],
            observed_dev_expl = round(dev_expl[run == 0], 1),
            n_null = sum(run > 0 & !is.na(n_terms)),
            null_pct_with_terms = round(100 * mean(n_terms[run > 0] > 0, na.rm = TRUE), 1),
            null_median_dev_expl = round(median(dev_expl[run > 0], na.rm = TRUE), 1),
            null_95th_dev_expl = round(quantile(dev_expl[run > 0], 0.95, na.rm = TRUE), 1),
            p_dev_expl = (1 + sum(dev_expl[run > 0] >= dev_expl[run == 0], na.rm = TRUE)) /
              (1 + n_null),
            p_n_terms = (1 + sum(n_terms[run > 0] >= n_terms[run == 0], na.rm = TRUE)) /
              (1 + n_null),
            .groups = "drop")
write_csv(summary_tbl, file.path(out_dir, "permutation_summary.csv"))
print(summary_tbl, width = Inf)
