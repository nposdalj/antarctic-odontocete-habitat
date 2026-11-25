# =============================================================================
# Script  : 09_chla_sensitivity.R
# Section : 2.6 Statistical analyses (supplement) - chlorophyll-a sensitivity
# Purpose : Chl-a was removed from the habitat-model candidates (decision
#           N. Posdaljian 2026-11-25, PI feedback: no mechanistic basis for a
#           specific lag; the 3-month lag had been chosen from the same data).
#           This script shows how the final habitat models change if chl-a is
#           offered as an extra candidate:
#             none : HYPOTHESIS_SETS (= main Stage 3 models)
#             lag0 : + same-day chl-a (30-day mean ending on the day)
#             lag3 : + chl-a 30-day mean ending 3 months earlier (as before)
#           chl-a enters as a standardized anomaly of log10 chl-a (10-yr
#           climatology, anomaly / monthly SD), exactly as the lag-3 variable
#           did. Same bins, ice filter, term caps and selection
#           (select_vif_checked) as 03_stage3_harmonized_models.R. Where a
#           chl-a term is retained, it is tested for robustness to AR(1)
#           (glmmTMB Tweedie, ar1(bin + 0 | series), LRT; as in
#           05_ar1_sensitivity.R).
#           Descriptive: prewhitened cross-correlation between daily presence
#           minutes and the standardized chl-a anomaly at lags 0-180 days.
#           The chl-a series is prewhitened with an AR(p) model (order by AIC,
#           fitted to the 10-year daily record of the site) and the same
#           filter is applied to the presence series (Box-Jenkins
#           prewhitening), so slow autocorrelation in both series does not
#           produce spurious peaks.
# Inputs  : GAM_MODEL_INPUT, GAM_ENV_DAILY (chla_0, chla_lag3),
#           stage1 binned_data/<species>_binned.csv, bin_lengths.csv,
#           stage3 tables/model_inclusion.csv
# Outputs : outputs/06_statistical_analyses/chla_sensitivity/
#             chla_sensitivity_models.csv (one row per species x site x variant)
#             chla_sensitivity_terms.csv  (final terms of every variant)
#             chla_ccf.csv, chla_ccf.png  (prewhitened cross-correlations)
# Author  : Natalie Posdaljian
# Date    : 2026-11-25
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Run from the repository root after 01 and 03.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))
suppressPackageStartupMessages({ library(glmmTMB); library(patchwork) })

out_dir <- file.path(STATS_DIR, "chla_sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
CHL_VARS <- c(lag0 = "chla_lag0", lag3 = "chla_lag3")
CCF_MAX_LAG <- 180
AR1_DELTA_AIC <- 2

# -----------------------------------------------------------------------------
# 1. Daily chl-a variables: log10 -> anomaly (10-yr climatology) -> standardized
# -----------------------------------------------------------------------------
env_daily <- read_csv(GAM_ENV_DAILY, show_col_types = FALSE) |>
  mutate(date = as.Date(date)) |> arrange(site, date)
model_input <- read_csv(GAM_MODEL_INPUT, show_col_types = FALSE) |>
  mutate(date = as.Date(date)) |> arrange(site, date)

# Same-day chl-a: 30-day mean ending on the day (same window as the lags),
# from the full daily record so the first recording days have a full window.
roll30 <- function(x) as.numeric(stats::filter(x, rep(1 / 30, 30), sides = 1))
env_daily <- env_daily |> group_by(site) |>
  mutate(chla_lag0 = roll30(chla_0)) |> ungroup()
model_input <- model_input |> select(-any_of("chla_lag0")) |>
  left_join(select(env_daily, site, date, chla_lag0), by = c("site", "date"))
stopifnot(!anyNA(model_input$chla_lag0))

env_daily   <- mutate(env_daily,   across(all_of(unname(CHL_VARS)), log10))
model_input <- mutate(model_input, across(all_of(unname(CHL_VARS)), log10))
clim <- build_climatology(filter(env_daily, !is.na(chla_lag0)), unname(CHL_VARS))
env_anom <- add_anomalies(env_daily, clim, unname(CHL_VARS))
model_input <- add_anomalies(model_input, clim, unname(CHL_VARS)) |>
  standardize_anomalies(env_anom, unname(CHL_VARS))
env_anom <- standardize_anomalies(env_anom, env_anom, unname(CHL_VARS))
chl_cols <- paste0(unname(CHL_VARS), "_anom")

# -----------------------------------------------------------------------------
# 2. Refit the habitat models with each chl-a variant
# -----------------------------------------------------------------------------
bin_lengths <- read_csv(file.path(STATS_DIR, "stage1", "tables", "bin_lengths.csv"),
                        show_col_types = FALSE)
inclusion <- read_csv(file.path(STATS_DIR, "stage3", "tables", "model_inclusion.csv"),
                      show_col_types = FALSE) |>
  filter(status != "descriptive only")

term_rhs <- function(tm) {
  if (tm %in% names(CYCLIC_VARS)) paste0(tm, "_sin2 + ", tm, "_cos2")
  else paste0("s(", tm, ", k = ", GAM_K, ")")
}
fit_tmb <- function(response, terms, data, ar1) {
  rhs <- if (length(terms) > 0) vapply(terms, term_rhs, character(1)) else "1"
  if (ar1) rhs <- c(rhs, "ar1(bin + 0 | series)")
  tryCatch(suppressWarnings(glmmTMB(as.formula(paste(response, "~", paste(rhs, collapse = " + "))),
                                    data = data, family = tweedie(link = "log"))),
           error = function(e) e)
}
ok_fit <- function(f) !inherits(f, "error") && is.finite(AIC(f))
# AR(1) robustness of one term (LRT with/without it; AR(1) used if dAIC < -2)
ar1_term_p <- function(dat, sp, terms, tm) {
  d <- dat |> mutate(bin = factor(bin_index, levels = 0:max(bin_index)), series = factor(1))
  for (cv in intersect(terms, names(CYCLIC_VARS))) {
    d[[paste0(cv, "_sin2")]] <- sin(2 * d[[cv]] * pi / 180)
    d[[paste0(cv, "_cos2")]] <- cos(2 * d[[cv]] * pi / 180)
  }
  no_ar <- fit_tmb(sp, terms, d, FALSE); ar <- fit_tmb(sp, terms, d, TRUE)
  use_ar <- ok_fit(ar) && ok_fit(no_ar) && AIC(ar) - AIC(no_ar) < -AR1_DELTA_AIC
  full <- if (use_ar) ar else no_ar
  if (!ok_fit(full)) return(tibble(ar1_used = NA, p_ref = NA_real_))
  red <- fit_tmb(sp, setdiff(terms, tm), d, use_ar)
  if (!ok_fit(red)) return(tibble(ar1_used = use_ar, p_ref = NA_real_))
  df <- attr(logLik(full), "df") - attr(logLik(red), "df")
  tibble(ar1_used = use_ar,
         p_ref = pchisq(2 * (as.numeric(logLik(full)) - as.numeric(logLik(red))), df = df,
                        lower.tail = FALSE))
}

model_rows <- list(); term_rows <- list()
for (i in seq_len(nrow(inclusion))) {
  sp <- inclusion$species[i]; st <- inclusion$site[i]; id <- model_id(sp, st)
  bin_days <- bin_lengths$bin_days[bin_lengths$species == sp]
  chl_binned <- model_input |> filter(site == st) |>
    bin_site_data(sp, bin_days, chl_cols) |>
    select(site, bin_index, all_of(chl_cols))
  dat <- read_csv(file.path(STATS_DIR, "stage1", "binned_data", paste0(sp, "_binned.csv")),
                  show_col_types = FALSE) |>
    filter(site == st, !ice_excluded) |>
    select(-any_of(chl_cols)) |>
    left_join(chl_binned, by = c("site", "bin_index"))
  low_sample <- inclusion$status[i] == "modeled (low sample)"
  cap <- term_cap(nrow(dat), low_sample)

  for (variant in c("none", names(CHL_VARS))) {
    chl <- if (variant == "none") character(0) else paste0(CHL_VARS[[variant]], "_anom")
    be <- select_vif_checked(dat, sp, c(HYPOTHESIS_SETS[[sp]], chl), max_terms = cap)
    s <- summary(be$model)
    terms <- if (be$no_significant_terms) character(0) else be$terms
    chl_kept <- length(chl) > 0 && chl %in% terms
    rob <- if (chl_kept) ar1_term_p(dat, sp, terms, chl) else tibble(ar1_used = NA, p_ref = NA_real_)
    model_rows[[paste(id, variant)]] <- tibble(
      species = sp, site = st, variant = variant, n_bins = nrow(dat),
      n_detection_bins = sum(dat[[sp]] > 0), term_cap = cap,
      final_terms = if (length(terms)) paste(terms, collapse = " + ") else "(none)",
      dev_expl_pct = round(100 * s$dev.expl, 1), aic = round(AIC(be$model), 1),
      chl_retained = chl_kept,
      chl_p = if (chl_kept) s$s.pv[term_variable(rownames(s$s.table)) == chl] else NA_real_,
      chl_ar1_used = rob$ar1_used, chl_p_ar1_ref = rob$p_ref,
      chl_holds_ar1 = chl_kept & !is.na(rob$p_ref) & rob$p_ref <= P_THRESH)
    if (length(terms)) {
      g <- map_chr(terms, \(v) { gr <- smooth_grid(be$model, v)
        classify_smooth_shape(gr$.estimate, s$s.pv[term_variable(rownames(s$s.table)) == v],
                              x = gr$x, cyclic = v %in% names(CYCLIC_VARS)) })
      term_rows[[paste(id, variant)]] <- tibble(
        species = sp, site = st, variant = variant, term = terms,
        edf = round(as.numeric(s$edf), 2), p_value = signif(as.numeric(s$s.pv), 3), shape = g)
    }
    message(id, " [", variant, "]: ", model_rows[[paste(id, variant)]]$final_terms,
            " (", model_rows[[paste(id, variant)]]$dev_expl_pct, "%)")
  }
}
models_tbl <- bind_rows(model_rows) |>
  group_by(species, site) |>
  mutate(same_as_none = final_terms == final_terms[variant == "none"],
         delta_aic_vs_none = aic - aic[variant == "none"]) |> ungroup()
write_csv(models_tbl, file.path(out_dir, "chla_sensitivity_models.csv"))
write_csv(bind_rows(term_rows), file.path(out_dir, "chla_sensitivity_terms.csv"))

# -----------------------------------------------------------------------------
# 3. Prewhitened cross-correlation (descriptive)
# -----------------------------------------------------------------------------
# z(t) = standardized anomaly of the 30-day mean chl-a ending on day t (lag0
# variable); the chl-a value k days before a recording day is z(t - k).
ccf_rows <- list()
for (st in SITE_ORDER) {
  e <- env_anom |> filter(site == st, !is.na(chla_lag0_anom)) |> arrange(date)
  # e must be a continuous daily record for the AR filter
  stopifnot(all(diff(e$date) == 1))
  arf <- ar(e$chla_lag0_anom, order.max = 30, aic = TRUE, method = "yule-walker")
  filt <- function(x) as.numeric(stats::filter(x - mean(x, na.rm = TRUE), c(1, -arf$ar),
                                               sides = 1, method = "convolution"))
  # filter() convolution with sides = 1 uses x[t], x[t-1], ...: coefficients
  # c(1, -phi) give the innovation series e_t = x_t - sum(phi_i x_{t-i}).
  z_f <- tibble(date = e$date, z_w = filt(e$chla_lag0_anom))
  rec <- model_input |> filter(site == st) |> arrange(date)
  for (sp in MODELED_SPECIES) {
    if (sum(rec[[sp]] > 0) < 10) next
    y_w <- filt(rec[[sp]])
    for (k in 0:CCF_MAX_LAG) {
      zk <- z_f$z_w[match(rec$date - k, z_f$date)]
      ok <- is.finite(y_w) & is.finite(zk)
      ccf_rows[[length(ccf_rows) + 1]] <- tibble(
        species = sp, site = st, lag_days = k, n = sum(ok),
        r = if (sum(ok) > 10) cor(y_w[ok], zk[ok]) else NA_real_,
        ar_order = arf$order)
    }
  }
}
ccf_tbl <- bind_rows(ccf_rows) |> mutate(bound = 1.96 / sqrt(n))
write_csv(ccf_tbl, file.path(out_dir, "chla_ccf.csv"))

p <- ccf_tbl |>
  mutate(species = factor(species_label(species), levels = species_label(MODELED_SPECIES)),
         site = factor(site, levels = SITE_ORDER)) |>
  ggplot(aes(lag_days, r)) +
  geom_ribbon(aes(ymin = -bound, ymax = bound), fill = "grey85") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 90, linetype = "dotted", linewidth = 0.3) +
  geom_line(aes(colour = site), linewidth = 0.5) +
  scale_colour_manual(values = SITE_COLORS, guide = "none") +
  facet_grid(species ~ site, labeller = label_wrap_gen(14)) +
  labs(x = "Lag (days; chl-a before presence)",
       y = "Cross-correlation (prewhitened)") +
  theme_bw(base_size = 9)
ggsave(file.path(out_dir, "chla_ccf.png"), p, width = 7, height = 7.5, dpi = 300)

# Lags outside the 95% bounds (expected ~5% by chance, i.e. ~9 of 181 lags)
print(ccf_tbl |> group_by(species, site) |>
        summarise(n_lags_outside = sum(abs(r) > bound, na.rm = TRUE),
                  max_abs_r = round(max(abs(r), na.rm = TRUE), 2),
                  lag_of_max = lag_days[which.max(abs(r))],
                  r_at_90 = round(r[lag_days == 90], 2), .groups = "drop"), n = Inf)
print(select(models_tbl, species, site, variant, final_terms, dev_expl_pct,
             delta_aic_vs_none, chl_p, chl_holds_ar1), n = Inf, width = Inf)
message("Chl-a sensitivity complete. Outputs in ", out_dir)
