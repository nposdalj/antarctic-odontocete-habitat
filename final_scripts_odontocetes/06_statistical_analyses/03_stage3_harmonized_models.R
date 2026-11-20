# =============================================================================
# Script  : 03_stage3_harmonized_models.R
# Section : 2.6 Statistical analyses (Stage 3: harmonized habitat models;
#           main-text results with the sea-ice models of 04_sea_ice_models.R)
# Purpose : TEMPLATE. Fit the manually harmonized habitat models: one term set
#           per species (or per species x site), chosen from the Stage 2
#           decision table, then backward-eliminate within it with the same
#           rules as Stage 1 (p <= 0.05, term cap floor(n_bins / 10), 2 terms
#           in low-sample models, "no significant terms" flag), so no final
#           model keeps a non-significant term. Uses the same binned,
#           ice-filtered data, family and settings as Stage 1, saves the same
#           diagnostics and writes the Stage 1 result format, so the
#           08_results_gams/ scripts can read either stage.
# Inputs  : <OUTPUT_ROOT>/06_statistical_analyses/stage1/
#             binned_data/<species>_binned.csv, tables/model_inclusion.csv,
#             tables/ice_cutoffs.csv
#           harmonized_terms (below) [INSERT after reviewing Stage 2]
# Outputs : <OUTPUT_ROOT>/06_statistical_analyses/stage3/
#             tables/model_inclusion.csv, harmonized_terms.csv,
#                    backward_elimination.csv, final_model_summary.csv,
#                    final_model_terms.csv, final_model_resid_acf.csv,
#                    final_model_concurvity.csv, final_model_k_check.csv
#             models/<species>_<site>.rds, diagnostics/<species>_<site>_*
# Author  : Natalie Posdaljian
# Date    : 2026-11-16 (updated 2026-11-17)
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Sea ice is not a habitat candidate (see 04_sea_ice_models.R).
#           Run after 01 and 02.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

# -----------------------------------------------------------------------------
# 1. Harmonized terms (chosen by N. Posdaljian, 2026-11-19)
# -----------------------------------------------------------------------------
# One row per species (site = "all": every modeled site of that species) or
# per species x site (a site-specific row overrides the "all" row). `terms`
# holds candidate column names exactly as in stage2_decision_table.csv; each
# becomes s(term, k = GAM_K). One predictor set per species, used at every
# site: the union of the Stage 1 selections across that species' sites, one
# per theme, with the common 3-month chl-a lag and layer-averaged water-column
# variables, temperature as the only water-mass variable (review page,
# version 7). No species had sites disagreeing within a theme.
# Approved by N. Posdaljian, 2026-11-19.
harmonized_terms <- tribble(
  ~species, ~site, ~terms,
  # Union of the Stage 1 picks per species (re-derived 2026-11-20 after the
  # switch to per-cell EKE and standardized chl-a anomalies)
  "Pm",   "all", c("FSLE", "fsle_orient", "mixed_layer_anom",
                   "temperature_surface_anom", "chla_lag3_anom"),
  "Gm",   "all", c("temperature_surface_anom", "chla_lag3_anom", "mixed_layer_anom"),
  "Oo",   "all", c("EKE_mid", "temperature_surface_anom"),
  "BW29", "all", c("EKE_deep", "FSLE", "chla_lag3_anom", "mixed_layer_anom",
                   "temperature_surface_anom"),
  "BW37", "all", c("FSLE", "fsle_orient", "mixed_layer_anom", "temperature_deep_anom")
)

if (nrow(harmonized_terms) == 0) {
  stop("[INSERT] harmonized_terms is empty: fill in the formulas after reviewing ",
       "the Stage 2 decision table.")
}

# -----------------------------------------------------------------------------
# 2. Inputs from Stage 1
# -----------------------------------------------------------------------------
s1_dir    <- file.path(STATS_DIR, "stage1")
out_dir   <- stage_dir("stage3")
table_dir <- stage_dir("stage3", "tables")

inclusion <- read_csv(file.path(s1_dir, "tables", "model_inclusion.csv"),
                      show_col_types = FALSE)
ice <- read_csv(file.path(s1_dir, "tables", "ice_cutoffs.csv"), show_col_types = FALSE) |>
  distinct(species, ice_cutoff)

# Resolve the term set of every modeled species x site: a site-specific row
# wins over the species-wide "all" row.
model_plan <- inclusion |>
  filter(status != "descriptive only") |>
  left_join(harmonized_terms |> filter(site != "all"), by = c("species", "site")) |>
  left_join(harmonized_terms |> filter(site == "all") |> select(species, terms_all = terms),
            by = "species") |>
  mutate(terms = map2(terms, terms_all, \(x, y) if (!is.null(x)) x else y)) |>
  select(-terms_all) |>
  filter(!map_lgl(terms, is.null))

skipped <- anti_join(inclusion |> filter(status != "descriptive only"), model_plan,
                     by = c("species", "site"))
if (nrow(skipped) > 0) {
  message("No harmonized terms given (not fitted): ",
          paste(skipped$species, skipped$site, collapse = ", "))
}

write_csv(model_plan |> mutate(terms = map_chr(terms, paste, collapse = " + ")),
          file.path(table_dir, "harmonized_terms.csv"))
write_csv(inclusion, file.path(table_dir, "model_inclusion.csv"))

# -----------------------------------------------------------------------------
# 3. Fit, eliminate and document the harmonized models
# -----------------------------------------------------------------------------
# Same binned data (bin length, ice filter, raw/anomaly columns) as Stage 1;
# backward elimination within the harmonized term set with the term caps.
docs <- list()
backward_tbl <- list()
for (i in seq_len(nrow(model_plan))) {
  sp <- model_plan$species[i]
  st <- model_plan$site[i]
  id <- model_id(sp, st)
  terms <- model_plan$terms[[i]]
  dat <- read_csv(file.path(s1_dir, "binned_data", paste0(sp, "_binned.csv")),
                  show_col_types = FALSE) |>
    filter(site == st, !ice_excluded)

  unknown <- setdiff(terms, names(dat))
  if (length(unknown) > 0) {
    stop(sp, " ", st, ": terms not among the binned candidates: ",
         paste(unknown, collapse = ", "))
  }

  low_sample <- model_plan$status[i] == "modeled (low sample)"
  cap <- term_cap(nrow(dat), low_sample)
  be <- backward_eliminate(dat, sp, terms, max_terms = cap)
  backward_tbl[[id]] <- mutate(be$steps, species = sp, site = st, .before = 1)
  message(id, ": ", deparse1(formula(be$model)),
          if (be$no_significant_terms) " (no significant terms)" else "")

  docs[[id]] <- document_final_model(
    be$model, dat, sp, st, stage = "stage3",
    meta = list(bin_days = model_plan$bin_days[i],
                ice_cutoff = ice$ice_cutoff[ice$species == sp],
                low_sample = low_sample, term_cap = cap,
                harmonized_terms = terms,
                no_significant_terms = be$no_significant_terms),
    out_dir = out_dir)
}

write_csv(bind_rows(backward_tbl), file.path(table_dir, "backward_elimination.csv"))
write_model_tables(docs, table_dir)
message("Stage 3 complete: ", length(docs), " models. Outputs in ", out_dir)
