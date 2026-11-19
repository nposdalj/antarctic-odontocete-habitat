# =============================================================================
# Script  : 04_sea_ice_models.R
# Section : 2.6 Statistical analyses (sea-ice models; main-text results)
# Purpose : Sea-ice model per species x site: presence ~ s(ice_conc, k = 4),
#           Tweedie (log link, p in [1.1, 1.9]), REML, select = TRUE, fitted to
#           ALL bins (no ice filter), to describe occurrence across the full
#           sea-ice range. Complements the habitat models (Stages 1-3), which
#           use ice-filtered bins and do not include sea ice as a candidate.
#           Replaces GAMs/seasonalModels/ (concentration only: no ice_diff, no
#           ice_thickness).
# Inputs  : <OUTPUT_ROOT>/06_statistical_analyses/stage1/
#             binned_data/<species>_binned.csv (species bin lengths, all bins,
#             written by 01_stage1_automated_selection.R), tables/bin_lengths.csv
# Outputs : <OUTPUT_ROOT>/06_statistical_analyses/sea_ice/
#             tables/model_inclusion.csv, final_model_summary.csv,
#                    final_model_terms.csv, final_model_resid_acf.csv,
#                    final_model_concurvity.csv, final_model_k_check.csv
#             models/<species>_<site>.rds, diagnostics/<species>_<site>_*
# Author  : Natalie Posdaljian
# Date    : 2026-11-17
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Same binning (species bin length, 4-day cap) and inclusion rule
#           (>= 10 bins with detections, here counted over all bins) as the
#           habitat models. The single term is always reported with its p,
#           edf and deviance explained; its shape is "not significant" when
#           p > 0.05 (classify_smooth_shape()). Kept as a separate script so
#           the habitat workflow (Stages 1-3) stays unchanged and the sea-ice
#           results can be re-run on their own. Run after 01.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

binned_dir <- file.path(STATS_DIR, "stage1", "binned_data")
bin_lengths <- read_csv(file.path(STATS_DIR, "stage1", "tables", "bin_lengths.csv"),
                        show_col_types = FALSE)
out_dir    <- stage_dir("sea_ice")
table_dir  <- stage_dir("sea_ice", "tables")

docs <- list()
inclusion <- list()
for (sp in MODELED_SPECIES) {
  # All bins of all sites (ice_excluded is ignored: no ice filter here).
  binned <- read_csv(file.path(binned_dir, paste0(sp, "_binned.csv")),
                     show_col_types = FALSE)
  bin_days <- bin_lengths$bin_days[bin_lengths$species == sp]

  for (st in SITE_ORDER) {
    id <- model_id(sp, st)
    dat <- binned |> filter(site == st)
    n_det <- sum(dat[[sp]] > 0)
    # Inclusion rule as for the habitat models; a site without any sea-ice
    # variation cannot support the model.
    status <- case_when(n_det < MIN_DETECTION_BINS   ~ "descriptive only",
                        sd(dat[[ICE_MODEL_VAR]]) == 0 ~ "no sea-ice variation",
                        n_det <= LOW_SAMPLE_MAX_BINS ~ "modeled (low sample)",
                        TRUE                         ~ "modeled")
    inclusion[[id]] <- tibble(species = sp, site = st, bin_days = bin_days,
                              n_bins = nrow(dat), n_detection_bins = n_det, status = status)
    message(id, ": ", nrow(dat), " bins, ", n_det, " with detections -> ", status)
    if (!startsWith(status, "modeled")) next

    mod <- fit_tw_gam(dat, sp, ICE_MODEL_VAR, select = TRUE)
    docs[[id]] <- document_final_model(
      mod, dat, sp, st, stage = "sea_ice",
      meta = list(bin_days = bin_days, ice_cutoff = NA_real_,
                  low_sample = status == "modeled (low sample)"),
      out_dir = out_dir)
  }
}

write_csv(bind_rows(inclusion), file.path(table_dir, "model_inclusion.csv"))
write_model_tables(docs, table_dir)
message("Sea-ice models complete: ", length(docs), " models. Outputs in ", out_dir)
