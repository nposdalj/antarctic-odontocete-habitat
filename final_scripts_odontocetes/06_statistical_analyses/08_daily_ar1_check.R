# =============================================================================
# Script  : 08_daily_ar1_check.R
# Section : 2.6 Statistical analyses (supplement) - daily-resolution check
# Purpose : Refit each final habitat model (Stage 3 terms) on DAILY data
#           instead of multi-day bins, with and without a first-order
#           autoregressive (AR(1)) process on consecutive days, in glmmTMB
#           (Tweedie, log link; smooths as s(x, k = GAM_K); FSLV orientation
#           as sin 2theta + cos 2theta), and test each term by a
#           likelihood-ratio test. Asks whether using all recording days with
#           explicit autocorrelation changes which terms are supported.
# Inputs  : GAM_MODEL_INPUT, GAM_ENV_DAILY (predictors as in Stage 1),
#           stage3 models (final terms), stage1 ice_cutoffs.csv
# Outputs : outputs/06_statistical_analyses/stage3/tables/
#             daily_ar1_terms.csv, daily_ar1_models.csv
# Author  : Natalie Posdaljian
# Date    : 2026-11-23
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Days are filtered by the species' ice cutoff on the daily ice
#           concentration. AR(1) is retained if it lowers AIC by > 2, as in
#           05_ar1_sensitivity.R. Exploratory check; the binned models remain
#           the main analysis unless decided otherwise.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))
suppressPackageStartupMessages(library(glmmTMB))
AR1_DELTA_AIC <- 2

# Predictors as in Stage 1 (log10 chl-a, anomalies, standardized chl-a)
model_input <- read_csv(GAM_MODEL_INPUT, show_col_types = FALSE) |>
  mutate(date = as.Date(date)) |> arrange(site, date)
env_daily <- read_csv(GAM_ENV_DAILY, show_col_types = FALSE) |> mutate(date = as.Date(date))
if (CHL_LOG10) {
  model_input <- mutate(model_input, across(any_of(LAG_COLUMNS), log10))
  env_daily   <- mutate(env_daily,   across(any_of(LAG_COLUMNS), log10))
}
anom_vars <- unique(unlist(map(MODELED_SPECIES, anomaly_sources)))
anom_vars <- intersect(anom_vars, intersect(names(env_daily), names(model_input)))
climatology <- build_climatology(env_daily, anom_vars)
model_input <- add_anomalies(model_input, climatology, anom_vars)
if (CHL_STANDARDIZE) {
  model_input <- standardize_anomalies(model_input, add_anomalies(env_daily, climatology, anom_vars),
                                       LAG_COLUMNS)
}
cutoffs <- read_csv(file.path(STATS_DIR, "stage1", "tables", "ice_cutoffs.csv"),
                    show_col_types = FALSE) |> distinct(species, ice_cutoff)

term_rhs <- function(tm) {
  if (tm %in% names(CYCLIC_VARS)) paste0(tm, "_sin2 + ", tm, "_cos2")
  else paste0("s(", tm, ", k = ", GAM_K, ")")
}
fit_tmb <- function(response, terms, data, ar1) {
  rhs <- if (length(terms) > 0) vapply(terms, term_rhs, character(1)) else "1"
  if (ar1) rhs <- c(rhs, "ar1(day + 0 | series)")
  fml <- as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
  tryCatch(suppressWarnings(glmmTMB(fml, data = data, family = tweedie(link = "log"))),
           error = function(e) e)
}
lrt_p <- function(full, response, terms, data, ar1) {
  map_dbl(terms, function(tm) {
    red <- fit_tmb(response, setdiff(terms, tm), data, ar1)
    if (inherits(red, "error") || !is.finite(as.numeric(logLik(red)))) return(NA_real_)
    df <- attr(logLik(full), "df") - attr(logLik(red), "df")
    pchisq(2 * (as.numeric(logLik(full)) - as.numeric(logLik(red))), df = df, lower.tail = FALSE)
  })
}

models <- read_stage_models("stage3")
res <- imap(models, function(m, id) {
  sp <- m$species; st <- m$site
  terms <- term_variable(rownames(summary(m$model)$s.table))
  if (length(terms) == 0) return(NULL)
  d <- model_input |>
    filter(site == st, ice_conc <= cutoffs$ice_cutoff[cutoffs$species == sp]) |>
    mutate(day = factor(as.integer(date - min(date)),
                        levels = 0:as.integer(max(date) - min(date))),
           series = factor(1))
  for (cv in intersect(terms, names(CYCLIC_VARS))) {
    d[[paste0(cv, "_sin2")]] <- sin(2 * d[[cv]] * pi / 180)
    d[[paste0(cv, "_cos2")]] <- cos(2 * d[[cv]] * pi / 180)
  }
  no_ar <- fit_tmb(sp, terms, d, FALSE)
  ar <- fit_tmb(sp, terms, d, TRUE)
  ar_ok <- !inherits(ar, "error") && is.finite(AIC(ar))
  daic <- if (ar_ok && !inherits(no_ar, "error")) AIC(ar) - AIC(no_ar) else NA_real_
  message(id, ": ", nrow(d), " days; dAIC(AR1) = ", round(daic, 1))
  list(model = tibble(species = sp, site = st, n_days = nrow(d),
                      n_detection_days = sum(d[[sp]] > 0), delta_aic_ar1 = daic,
                      ar1_supported = !is.na(daic) & daic < -AR1_DELTA_AIC),
       terms = tibble(species = sp, site = st, term = terms,
                      p_binned_gam = as.numeric(summary(m$model)$s.pv),
                      p_daily_no_ar1 = if (inherits(no_ar, "error")) NA_real_ else
                        lrt_p(no_ar, sp, terms, d, FALSE),
                      p_daily_ar1 = if (ar_ok) lrt_p(ar, sp, terms, d, TRUE) else NA_real_))
}) |> compact()

model_tbl <- map(res, "model") |> list_rbind()
term_tbl <- map(res, "terms") |> list_rbind() |>
  left_join(select(model_tbl, species, site, ar1_supported), by = c("species", "site")) |>
  mutate(across(starts_with("p_"), function(x) as.numeric(unlist(x))),
         p_daily_reference = if_else(ar1_supported, p_daily_ar1, p_daily_no_ar1),
         holds_daily = p_daily_reference <= P_THRESH,
         across(starts_with("p_"), \(x) signif(x, 3)))
table_dir <- stage_dir("stage3", "tables")
write_csv(model_tbl, file.path(table_dir, "daily_ar1_models.csv"))
write_csv(term_tbl, file.path(table_dir, "daily_ar1_terms.csv"))
print(model_tbl, width = Inf); print(term_tbl, n = Inf, width = Inf)
