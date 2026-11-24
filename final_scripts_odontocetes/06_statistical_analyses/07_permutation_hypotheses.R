# =============================================================================
# Script  : 07_permutation_hypotheses.R
# Section : 2.6 Statistical analyses (supplement) - time-shift permutation
#           tests of (A) the Stage 3 models and (B) a priori hypotheses
# Purpose : Same null as 06_permutation_test.R: the predictors on each
#           recording day are replaced by the same site's predictors from
#           another year at the same time of year (whole years +/- up to 60
#           days); presence and sea ice (and thus the ice filter and bins) stay
#           on the real dates. The observed deviance explained is compared with
#           its distribution over N_PERM shifted runs.
#   Mode A  : Stage 3 procedure. Backward elimination (with the term cap)
#             from each species' harmonized Stage 3 set, repeated in every
#             shifted run, so the null gets the same selection step.
#   Mode B  : A priori hypotheses (HYPOTHESES below). Each is one fixed
#             single-term model per species x site, with no selection; the
#             statistic is the deviance explained by that term.
# Inputs  : GAM_MODEL_INPUT, GAM_ENV_DAILY, stage1 model_inclusion.csv and
#           ice_cutoffs.csv, stage3 harmonized_terms.csv
# Outputs : outputs/06_statistical_analyses/permutation/
#             A_stage3_runs.csv, A_stage3_summary.csv
#             B_hypotheses_runs.csv, B_hypotheses_summary.csv (with Holm and
#             Benjamini-Hochberg adjusted permutation p-values)
# Usage   : Rscript 07_permutation_hypotheses.R A    (or B)
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : N_PERM via env ODONTOCETE_N_PERM (default 200). The same shifts
#           are used in A and B (fixed seed per species x site).
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

MODE <- toupper(commandArgs(trailingOnly = TRUE)[1])
if (!MODE %in% c("A", "B", "C")) stop("Usage: Rscript 07_permutation_hypotheses.R A|B|C")
OUT_NAME   <- c(A = "stage3", B = "hypotheses", C = "hypothesis_set")[[MODE]]
# Mode C: backward elimination (with the term cap) from each species' a priori
# hypothesis set (HYPOTHESES), repeated in every shifted run.
N_PERM     <- as.integer(Sys.getenv("ODONTOCETE_N_PERM", unset = "200"))
MAX_OFFSET <- 60

# A priori hypotheses (decision N. Posdaljian 2026-11-20). Each term is tested
# alone at every modeled site of the species. Rationale and references are in
# the manuscript (Supplement, permutation tests).
HYPOTHESES <- list(
  Pm   = c("EKE_mid", "EKE_deep", "temperature_mid_anom", "FSLE", "chla_lag3_anom"),
  Gm   = c("temperature_surface_anom", "FSLE", "mixed_layer_anom", "chla_lag3_anom"),
  Oo   = c("temperature_surface_anom", "mixed_layer_anom", "EKE_mid", "chla_lag3_anom"),
  BW29 = c("temperature_deep_anom", "EKE_deep", "FSLE", "chla_lag3_anom"),
  BW37 = c("temperature_deep_anom", "EKE_deep", "FSLE", "chla_lag3_anom")
)

out_dir <- file.path(STATS_DIR, "permutation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Predictors as in Stage 1 (log10 chl-a, anomalies, standardized chl-a)
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

s1_tables <- file.path(STATS_DIR, "stage1", "tables")
inclusion <- read_csv(file.path(s1_tables, "model_inclusion.csv"), show_col_types = FALSE) |>
  filter(status != "descriptive only")
cutoffs <- read_csv(file.path(s1_tables, "ice_cutoffs.csv"), show_col_types = FALSE) |>
  distinct(species, ice_cutoff)
harmonized <- read_csv(file.path(STATS_DIR, "stage3", "tables", "harmonized_terms.csv"),
                       show_col_types = FALSE) |>
  transmute(species, site, terms = str_split(terms, " [+] "))

# -----------------------------------------------------------------------------
# 2. Shifts and statistics
# -----------------------------------------------------------------------------
shift_predictors <- function(daily_site, env_site, shift_days, cols) {
  src <- env_site |> select(date, all_of(cols)) |> mutate(date = date - shift_days)
  daily_site |> select(-all_of(cols)) |> left_join(src, by = "date")
}
draw_shifts <- function(dates, env_dates, n) {
  years <- setdiff(-10:10, 0)
  out <- integer(0)
  while (length(out) < n) {
    s <- round(365.25 * sample(years, 1)) + sample(-MAX_OFFSET:MAX_OFFSET, 1)
    if (all((dates + s) %in% env_dates)) out <- c(out, s)
  }
  out
}
dev_of <- function(dat, sp, terms) 100 * summary(fit_tw_gam(dat, sp, terms, select = TRUE))$dev.expl

# Statistics for one (real or shifted) data set: one row per test
run_stats <- function(daily_site, sp, st, bin_days, cutoff, low_sample, terms_A, terms_B) {
  cols <- unique(c(terms_A, terms_B))
  dat <- bin_site_data(daily_site, sp, bin_days, c(cols, ICE_MODEL_VAR)) |>
    filter(ice_conc <= cutoff)
  if (MODE %in% c("A", "C")) {
    # Same procedure as 03_stage3_harmonized_models.R: backward selection with
    # the term cap and a VIF check on the final model
    be <- select_vif_checked(dat, sp, terms_A, max_terms = term_cap(nrow(dat), low_sample))
    tibble(test = "stage3_selection",
           terms = if (be$no_significant_terms) "" else paste(be$terms, collapse = " + "),
           dev_expl = if (be$no_significant_terms) 0 else
             100 * summary(be$model)$dev.expl)
  } else {
    map_dfr(terms_B, \(v) tibble(test = v, terms = v, dev_expl = dev_of(dat, sp, v)))
  }
}

# -----------------------------------------------------------------------------
# 3. Observed and shifted runs per species x site
# -----------------------------------------------------------------------------
runs <- list()
for (i in seq_len(nrow(inclusion))) {
  sp <- inclusion$species[i]; st <- inclusion$site[i]
  set.seed(20261120 + i)   # same shifts in modes A and B
  daily_site <- filter(model_input, site == st)
  env_site <- filter(env_daily, site == st)
  terms_A <- if (MODE == "C") intersect(HYPOTHESIS_SETS[[sp]], names(daily_site)) else
    harmonized$terms[[which(harmonized$species == sp & harmonized$site == st)]]
  terms_B <- intersect(HYPOTHESES[[sp]], names(daily_site))
  cols <- unique(c(terms_A, terms_B))
  args <- list(sp = sp, st = st, bin_days = inclusion$bin_days[i],
               cutoff = cutoffs$ice_cutoff[cutoffs$species == sp],
               low_sample = inclusion$status[i] == "modeled (low sample)",
               terms_A = terms_A, terms_B = terms_B)
  shifts <- draw_shifts(daily_site$date, env_site$date, N_PERM)
  obs <- do.call(run_stats, c(list(daily_site = daily_site), args)) |> mutate(run = 0L)
  null <- map_dfr(seq_along(shifts), function(k) {
    tryCatch(do.call(run_stats, c(list(daily_site = shift_predictors(daily_site, env_site,
                                                                     shifts[k], cols)), args)),
             error = function(e) tibble(test = NA_character_, terms = conditionMessage(e),
                                        dev_expl = NA_real_)) |>
      mutate(run = k)
  })
  runs[[paste(sp, st)]] <- bind_rows(obs, null) |> mutate(species = sp, site = st, .before = 1)
  message(sp, " ", st, " done")
  write_csv(bind_rows(runs), file.path(out_dir, paste0(MODE, "_", OUT_NAME, "_runs.csv")))
}

runs <- bind_rows(runs)
summary_tbl <- runs |>
  filter(!is.na(test)) |>
  group_by(species, site, test) |>
  summarise(observed_terms = terms[run == 0],
            observed_dev_expl = round(dev_expl[run == 0], 1),
            n_null = sum(run > 0),
            null_median_dev_expl = round(median(dev_expl[run > 0], na.rm = TRUE), 1),
            null_95th_dev_expl = round(quantile(dev_expl[run > 0], 0.95, na.rm = TRUE), 1),
            p_perm = (1 + sum(dev_expl[run > 0] >= dev_expl[run == 0], na.rm = TRUE)) /
              (1 + n_null),
            .groups = "drop") |>
  mutate(p_holm = p.adjust(p_perm, "holm"), p_bh = p.adjust(p_perm, "BH"))
write_csv(summary_tbl, file.path(out_dir, paste0(MODE, "_", OUT_NAME, "_summary.csv")))
print(summary_tbl, n = Inf, width = Inf)
