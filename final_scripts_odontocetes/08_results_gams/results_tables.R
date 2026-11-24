# =============================================================================
# Script  : results_tables.R
# Section : 3.2 GAM results
# Purpose : Manuscript tables of the final models.
#           Sea-ice models (main text): per species x site n bins, bins with
#           detections, edf, p, deviance explained, shape, direction.
#           Habitat models (main text = Stage 3; Stage 1 = supplementary):
#           selected terms with edf and p, deviance explained by the model and
#           by each term (drop-one refit), n bins, bins with detections,
#           low-sample flag; species x site combinations below the inclusion
#           rule are listed as "descriptive only".
#           Diagnostics tables (supplementary, both model types): bin length,
#           ice cutoff, term cap, estimated Tweedie power, minimum k-index
#           p-value, maximum concurvity, residual ACF check.
# Inputs  : <OUTPUT_ROOT>/06_statistical_analyses/<stage>/
#             models/<species>_<site>.rds, tables/model_inclusion.csv
#           for stage = sea_ice and RESULTS_STAGE (stage1 or stage3)
# Outputs : <OUTPUT_ROOT>/08_results_gams/
#             main_sea_ice_gam_table.csv
#             supplementary_sea_ice_gam_diagnostics.csv
#             <label>_habitat_gam_table.csv        (one row per model)
#             <label>_habitat_gam_table_terms.csv  (one row per model term)
#             supplementary_<stage>_habitat_gam_diagnostics.csv
#           <label> = "main_stage3" or "supplementary_stage1"
# Author  : Natalie Posdaljian
# Date    : 2026-11-16 (updated 2026-11-17)
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Habitat models whose last remaining term had p > 0.05 are listed
#           as "no significant terms" (no terms shown; diagnostics kept).
#           Residual ACF results are reported as a limitation (decision 2.6).
#           Residual ACF OK = no lag in 2..RESID_ACF_MAX_LAG bins outside the
#           95% bounds of the deviance-residual ACF (ice-removed bins kept as
#           gaps); lag 1 is reported separately (resid_acf_lag1).
#           Per-term deviance explained = model deviance explained minus that
#           of the model refitted without the term, so terms need not sum to
#           the model total.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

# --- Parameters ---------------------------------------------------------------
# Habitat-model stage: "stage3" = harmonized models (main text); "stage1" =
# automated selection (supplement). The sea-ice models are always included.
RESULTS_STAGE <- Sys.getenv("ODONTOCETE_RESULTS_STAGE", unset = "stage3")
HABITAT_LABEL <- if (RESULTS_STAGE == "stage3") "main_stage3" else
  paste0("supplementary_", RESULTS_STAGE)

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
out_file <- function(name) file.path(RESULTS_DIR, name)

# p-values formatted as in the manuscript tables.
format_p <- function(p) if_else(p < 0.001, "<0.001", sprintf("%.3f", p))

order_rows <- function(df) {
  df |> arrange(factor(species, levels = MODELED_SPECIES), factor(site, levels = SITE_ORDER))
}

# Model summaries, term tables and inclusion table of one stage.
read_stage <- function(stage) {
  models <- read_stage_models(stage)
  list(summaries = bind_rows(map(models, "summary")),
       terms = bind_rows(map(models, "terms")),
       inclusion = read_csv(file.path(STATS_DIR, stage, "tables", "model_inclusion.csv"),
                            show_col_types = FALSE))
}

# AR(1) sensitivity results of one stage (05_ar1_sensitivity.R), if present:
# per model (AIC change, AR(1) supported) and per term (robust or not).
read_ar1 <- function(stage) {
  d <- file.path(STATS_DIR, stage, "tables")
  if (!file.exists(file.path(d, "ar1_sensitivity_models.csv"))) return(NULL)
  list(models = read_csv(file.path(d, "ar1_sensitivity_models.csv"), show_col_types = FALSE) |>
         select(species, site, delta_aic_ar1, ar1_supported),
       terms = read_csv(file.path(d, "ar1_sensitivity_terms.csv"), show_col_types = FALSE) |>
         transmute(species, site, variable = term_variable(term), robust_ar1 = holds))
}

# Supplementary diagnostics table (same columns for both model types).
# Residual autocorrelation: Ljung-Box test over lags 1-5 (OK = p > 0.05).
diagnostics_table <- function(summaries, ar1 = NULL) {
  out <- summaries |>
    transmute(species, species_name = species_label(species), site, formula,
              bin_days, ice_cutoff = round(ice_cutoff, 3), term_cap,
              tweedie_p = round(tweedie_p, 2),
              k_index_min = round(k_index_min, 2), k_pvalue_min = round(k_pvalue_min, 3),
              concurvity_max = round(concurvity_max, 2),
              resid_acf_lag1 = round(resid_acf_lag1, 2),
              ljung_box_p = signif(ljung_box_p, 2), resid_acf_ok,
              no_significant_terms, all_terms_significant, AIC = round(AIC, 1))
  if (!is.null(ar1)) {
    out <- out |>
      left_join(ar1$models |> mutate(delta_aic_ar1 = round(delta_aic_ar1, 1)),
                by = c("species", "site"))
  }
  order_rows(out)
}

# -----------------------------------------------------------------------------
# 1. Sea-ice models (main text)
# -----------------------------------------------------------------------------
# One row per species x site, including those not modeled.
ice <- read_stage("sea_ice")
ice_ar1 <- read_ar1("sea_ice")
# Shape of every fitted sea-ice smooth, also when not significant (reported
# as a weak, non-significant tendency; decision N. Posdaljian 2026-11-20).
# Models may keep concentration and/or distance to the ice edge: one row per
# term (decision N. Posdaljian 2026-11-20).
ice_fitted_shape <- map_dfr(read_stage_models("sea_ice"), function(obj) {
  map_dfr(obj$terms$variable, function(v) {
    g <- smooth_grid(obj$model, v)
    tibble(species = obj$species, site = obj$site, variable = v,
           fitted_shape = classify_smooth_shape(g$.estimate, p_value = 0, x = g$x))
  })
})
sea_ice_table <- ice$inclusion |>
  left_join(ice$terms |> select(species, site, variable, edf, p_value, shape, direction),
            by = c("species", "site")) |>
  left_join(ice_fitted_shape, by = c("species", "site", "variable")) |>
  left_join(ice$summaries |> select(species, site, dev_expl), by = c("species", "site")) |>
  left_join(if (!is.null(ice_ar1)) ice_ar1$terms else
              tibble(species = character(), site = character(), variable = character(),
                     robust_ar1 = logical()),
            by = c("species", "site", "variable")) |>
  transmute(species, species_name = species_label(species), site, status, variable,
            n_bins, n_detection_bins, low_sample = status == "modeled (low sample)",
            edf = round(edf, 2), p_value = format_p(p_value),
            dev_expl_pct = round(100 * dev_expl, 1), shape, direction, fitted_shape,
            robust_ar1 = if_else(p_value <= P_THRESH, robust_ar1, NA)) |>
  order_rows()
write_csv(sea_ice_table, out_file("main_sea_ice_gam_table.csv"))
write_csv(diagnostics_table(ice$summaries, ice_ar1),
          out_file("supplementary_sea_ice_gam_diagnostics.csv"))

# -----------------------------------------------------------------------------
# 2. Habitat models (Stage 3 main text; Stage 1 supplementary)
# -----------------------------------------------------------------------------
hab <- read_stage(RESULTS_STAGE)
hab_ar1 <- read_ar1(RESULTS_STAGE)

# Long form: one row per term, ready for formatting in the manuscript.
main_terms <- hab$terms |>
  left_join(hab$summaries |>
              select(species, site, dev_expl, n_bins, n_detection_bins, low_sample),
            by = c("species", "site")) |>
  left_join(if (!is.null(hab_ar1)) hab_ar1$terms else
              tibble(species = character(), site = character(), variable = character(),
                     robust_ar1 = logical()),
            by = c("species", "site", "variable")) |>
  transmute(species, species_name = species_label(species), site,
            n_bins, n_detection_bins, low_sample,
            dev_expl_model_pct = round(100 * dev_expl, 1),
            term = variable_label(variable), variable,
            edf = round(edf, 2), p_value = format_p(p_value),
            dev_expl_term_pct = round(100 * dev_expl_term, 1),
            direction, shape, robust_ar1) |>
  order_rows()
write_csv(main_terms, out_file(paste0(HABITAT_LABEL, "_habitat_gam_table_terms.csv")))

# Compact form: one row per species x site, including combinations that were
# not modeled (inclusion rule) or have no significant terms.
main_compact <- main_terms |>
  group_by(species, site, dev_expl_model_pct) |>
  summarise(terms = paste0(term, " (edf ", edf, ", p ", p_value, ", ",
                           dev_expl_term_pct, "%)", collapse = "; "),
            .groups = "drop")
main_table <- hab$inclusion |>
  filter(species %in% MODELED_SPECIES) |>
  left_join(main_compact, by = c("species", "site")) |>
  left_join(hab$summaries |> select(species, site, no_significant_terms, term_cap),
            by = c("species", "site")) |>
  transmute(species, species_name = species_label(species), site,
            status = case_when(
              coalesce(no_significant_terms, FALSE) ~
                paste0(status, ": no significant terms"),
              status != "descriptive only" & is.na(no_significant_terms) ~
                paste0(status, ": not fitted in ", RESULTS_STAGE),
              TRUE ~ status),
            n_bins, n_detection_bins, low_sample = grepl("low sample", status),
            term_cap, dev_expl_model_pct, terms) |>
  order_rows()
write_csv(main_table, out_file(paste0(HABITAT_LABEL, "_habitat_gam_table.csv")))
write_csv(diagnostics_table(hab$summaries, hab_ar1),
          out_file(paste0("supplementary_", RESULTS_STAGE, "_habitat_gam_diagnostics.csv")))

message("Result tables (sea ice + habitat ", RESULTS_STAGE, ") written to ", RESULTS_DIR)
