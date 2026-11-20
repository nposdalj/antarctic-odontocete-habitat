# =============================================================================
# Script  : rare_detections.R
# Section : 3.1 / 3.2 Results - qualitative description of species x site
#           combinations that were not modeled (< 10 bins with detections)
# Purpose : For each unmodeled combination (killer whales at EI and KGI, BW37
#           at CI), list every detection day with the environmental
#           conditions on that day, and compare them with the conditions on
#           all recording days at that site (median and range), so the
#           Results can describe when and under what conditions these rare
#           detections occurred.
# Inputs  : F_PRESENCE_DAILY (presence minutes per species and day)
#           F_ENV_DAILY      (daily predictors; anomalies computed here as in
#                             06_statistical_analyses)
#           <OUTPUT_ROOT>/06_statistical_analyses/stage1/tables/model_inclusion.csv
# Outputs : outputs/07_results_occurrence/rare_detection_days.csv
#             one row per detection day: species, site, date, month,
#             presence minutes, ice, temperature anomalies by layer, MLD
#             anomaly, FSLE, EKE by layer, chl-a anomaly
#           outputs/07_results_occurrence/rare_detection_context.csv
#             per species x site x variable: median on detection days vs
#             median and range on all recording days, and the percentile of
#             the detection-day median among recording days
# Author  : Natalie Posdaljian
# Date    : 2026-11-19
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Descriptive only; no tests (3-7 detection bins).
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

VARS <- c("ice_conc", "temperature_surface_anom", "temperature_mid_anom",
          "temperature_deep_anom", "mixed_layer_anom", "FSLE", "EKE_mid",
          "EKE_deep", "chla_lag3_anom")

# Combinations not modeled, excluding those with no detections at all
rare <- read_csv(file.path(STATS_DIR, "stage1", "tables", "model_inclusion.csv"),
                 show_col_types = FALSE) |>
  filter(status == "descriptive only", n_detection_bins > 0) |>
  distinct(species, site)

# Predictors as used in the models
env <- read_csv(F_ENV_DAILY, show_col_types = FALSE) |> mutate(date = as.Date(date))
if (CHL_LOG10) env <- mutate(env, across(any_of(LAG_COLUMNS), log10))
anom_sources <- intersect(sub("_anom$", "", grep("_anom$", VARS, value = TRUE)), names(env))
env <- add_anomalies(env, build_climatology(env, anom_sources), anom_sources)
if (CHL_STANDARDIZE) env <- standardize_anomalies(env, env, intersect(LAG_COLUMNS, anom_sources))

presence <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) |>
  filter(rec_minutes > 0) |>
  mutate(date = as.Date(date))

days <- presence |>
  pivot_longer(any_of(MODELED_SPECIES), names_to = "species", values_to = "presence_min") |>
  semi_join(rare, by = c("species", "site")) |>
  left_join(select(env, site, date, any_of(VARS)), by = c("site", "date"))

detection_days <- days |>
  filter(presence_min > 0) |>
  mutate(month = format(date, "%b")) |>
  select(species, site, date, month, presence_min, any_of(VARS)) |>
  mutate(across(any_of(VARS), ~ signif(.x, 3))) |>
  arrange(species, site, date)
write_csv(detection_days, out_path("07_results_occurrence", "rare_detection_days.csv"))

context <- days |>
  pivot_longer(any_of(VARS), names_to = "variable") |>
  group_by(species, site, variable) |>
  summarise(n_detection_days = sum(presence_min > 0),
            median_detection_days = median(value[presence_min > 0], na.rm = TRUE),
            median_all_days = median(value, na.rm = TRUE),
            min_all_days = min(value, na.rm = TRUE),
            max_all_days = max(value, na.rm = TRUE),
            percentile_of_detection_median =
              round(100 * mean(value <= median(value[presence_min > 0], na.rm = TRUE),
                               na.rm = TRUE)),
            .groups = "drop") |>
  mutate(across(c(median_detection_days, median_all_days, min_all_days, max_all_days),
                ~ signif(.x, 3)))
write_csv(context, out_path("07_results_occurrence", "rare_detection_context.csv"))
print(detection_days, n = Inf, width = Inf)
print(context, n = Inf, width = Inf)
