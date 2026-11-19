# =============================================================================
# Script  : sensitivity_response_metric.R
# Section : 2.6 Statistical analyses (sensitivity analysis; supplement)
# Purpose : Test whether the choice of daily presence metric changes the
#           Stage 1 habitat models. Three daily responses are compared:
#             minutes - minutes per day in positive 5-min bins (current
#                       response; 0-1440)
#             hours   - hours per day with at least one positive 5-min bin
#                       (0-24)
#             days    - daily presence/absence (0/1)
#           Step 1 builds a model-input file for each alternative metric (the
#           current file with the species columns replaced). Step 2 re-runs
#           01_stage1_automated_selection.R unchanged on each file, writing
#           to a separate output root. Step 3 compares the three runs: bin-
#           level agreement between metrics, bin lengths, selected terms,
#           and deviance explained.
# Inputs  : F_PRESENCE_5MIN  (presence_5min.csv)
#           F_MODEL_INPUT    (model_input_daily.csv)
#           outputs/06_statistical_analyses/stage1/ (current run = minutes)
# Outputs : outputs/sensitivity_response_metric/
#             model_input_<metric>.csv
#             <metric>/06_statistical_analyses/stage1/... (full Stage 1 runs)
#             bin_metric_correlation.csv  - Spearman r between metrics per
#                                           species x site (4-day/2-day bins)
#             model_comparison.csv        - per species x site x metric:
#                                           bin_days, terms, dev_expl
#             term_agreement.csv          - per species x site: shared terms
# Author  : Natalie Posdaljian
# Date    : 2026-11-19
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - The ice cutoff (max ice at any detection) and the inclusion rule
#             (>= 10 bins with detections) depend only on presence vs absence,
#             so they are identical across metrics. Bin length (from the
#             residual ACF of the response) can differ.
#           - Every hour contains recording at all three sites (KGI records
#             10 of every 13.33 min), so hourly counts need no effort scaling.
#           - Run from the repository root:
#               Rscript final_scripts_odontocetes/06_statistical_analyses/sensitivity_response_metric.R
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

suppressPackageStartupMessages(library(tidyverse))

out_root <- out_path("sensitivity_response_metric", "")
modeled  <- MODELED_SPECIES
r_exe    <- file.path(R.home("bin"), "Rscript")
stage1_script <- file.path(FINAL_DIR, "06_statistical_analyses",
                           "01_stage1_automated_selection.R")

# -----------------------------------------------------------------------------
# 1. Alternative daily metrics from the 5-min presence grid
# -----------------------------------------------------------------------------
presence_5min <- read_csv(F_PRESENCE_5MIN, show_col_types = FALSE) |>
  filter(recorded) |>
  mutate(date = as.Date(bin_start), hour = lubridate::hour(bin_start))

hourly <- presence_5min |>
  group_by(site, date, hour) |>
  summarise(across(all_of(modeled), any), .groups = "drop")

daily_metrics <- list(
  hours = hourly |>
    group_by(site, date) |>
    summarise(across(all_of(modeled), sum), .groups = "drop"),
  days = hourly |>
    group_by(site, date) |>
    summarise(across(all_of(modeled), ~ as.numeric(any(.x))),
              .groups = "drop")
)

model_input <- read_csv(F_MODEL_INPUT, show_col_types = FALSE)

for (metric in names(daily_metrics)) {
  alt <- model_input |>
    select(-all_of(modeled)) |>
    left_join(daily_metrics[[metric]], by = c("site", "date"))
  # Every model-input day must have a value for the new metric
  stopifnot(!anyNA(alt[modeled]), nrow(alt) == nrow(model_input))
  write_csv(alt, file.path(out_root, paste0("model_input_", metric, ".csv")))
}

# -----------------------------------------------------------------------------
# 2. Re-run Stage 1 with each metric (same code, separate output root)
# -----------------------------------------------------------------------------
# The overrides are set with Sys.setenv() (inherited by the child process):
# system2(env = ...) does not work on Windows.
for (metric in names(daily_metrics)) {
  # Skip a metric whose run already finished (delete its folder to re-run)
  done_file <- file.path(out_root, metric, "06_statistical_analyses", "stage1",
                         "tables", "final_model_summary.csv")
  if (file.exists(done_file)) {
    message("Stage 1 (", metric, ") already done; skipping")
    next
  }
  message("Stage 1 with response = ", metric)
  log_file <- file.path(out_root, paste0("stage1_", metric, ".log"))
  Sys.setenv(
    ODONTOCETE_OUTPUT_ROOT = file.path(out_root, metric),
    ODONTOCETE_MODEL_INPUT = file.path(out_root,
                                       paste0("model_input_", metric, ".csv"))
  )
  status <- system2(r_exe, shQuote(stage1_script),
                    stdout = log_file, stderr = log_file)
  Sys.unsetenv(c("ODONTOCETE_OUTPUT_ROOT", "ODONTOCETE_MODEL_INPUT"))
  if (status != 0) stop("Stage 1 failed for ", metric, "; see ", log_file)
}

# -----------------------------------------------------------------------------
# 3. Compare
# -----------------------------------------------------------------------------
stage1_dir <- function(metric) {
  root <- if (metric == "minutes") OUT_DIR else file.path(out_root, metric)
  file.path(root, "06_statistical_analyses", "stage1")
}
metrics <- c("minutes", "hours", "days")

# 3a. Bin-level agreement: the current bins (minutes run) with each metric
# averaged over the same days.
daily_all <- bind_rows(
  model_input |> select(site, date, all_of(modeled)) |>
    mutate(metric = "minutes"),
  list_rbind(imap(daily_metrics, \(d, m) mutate(d, metric = m)))
) |>
  pivot_longer(all_of(modeled), names_to = "species", values_to = "value")

bin_lengths <- read_csv(file.path(stage1_dir("minutes"), "tables",
                                  "bin_lengths.csv"), show_col_types = FALSE)

bin_correlation <- daily_all |>
  inner_join(model_input |> group_by(site) |>
               summarise(first_day = min(date)), by = "site") |>
  inner_join(bin_lengths |> select(species, bin_days), by = "species") |>
  mutate(bin_index = as.integer(date - first_day) %/% bin_days) |>
  group_by(species, site, metric, bin_index) |>
  summarise(value = mean(value), .groups = "drop") |>
  pivot_wider(names_from = metric, values_from = value) |>
  group_by(species, site) |>
  filter(sum(minutes > 0) > 0) |>
  summarise(n_bins = n(),
            r_minutes_hours = cor(minutes, hours, method = "spearman"),
            r_minutes_days  = cor(minutes, days, method = "spearman"),
            .groups = "drop")
write_csv(bin_correlation, file.path(out_root, "bin_metric_correlation.csv"))

# 3b. Final models from each run
read_run <- function(metric) {
  d <- stage1_dir(metric)
  summary_tbl <- read_csv(file.path(d, "tables", "final_model_summary.csv"),
                          show_col_types = FALSE)
  terms_tbl <- read_csv(file.path(d, "tables", "final_model_terms.csv"),
                        show_col_types = FALSE) |>
    group_by(species, site) |>
    summarise(terms = paste(sort(variable), collapse = " + "),
              .groups = "drop")
  summary_tbl |>
    select(species, site, bin_days, n_bins, n_detection_bins,
           dev_expl, resid_acf_ok) |>
    left_join(terms_tbl, by = c("species", "site")) |>
    mutate(metric = metric, terms = replace_na(terms, "(none)"))
}

comparison <- map(metrics, read_run) |>
  list_rbind() |>
  mutate(metric = factor(metric, levels = metrics),
         dev_expl = round(dev_expl, 3)) |>
  arrange(species, site, metric)
write_csv(comparison, file.path(out_root, "model_comparison.csv"))

# 3c. Term agreement with the current (minutes) model
term_sets <- comparison |>
  select(species, site, metric, terms) |>
  mutate(term_list = str_split(terms, fixed(" + ")))

term_agreement <- term_sets |>
  filter(metric != "minutes") |>
  inner_join(term_sets |> filter(metric == "minutes") |>
               select(species, site, minutes_terms = term_list),
             by = c("species", "site")) |>
  mutate(identical = map2_lgl(term_list, minutes_terms,
                              \(a, b) setequal(a, b)),
         n_shared  = map2_int(term_list, minutes_terms,
                              \(a, b) length(intersect(a, b))),
         n_terms_minutes = lengths(minutes_terms),
         n_terms_alt     = lengths(term_list)) |>
  select(species, site, metric, identical, n_shared, n_terms_minutes,
         n_terms_alt)
write_csv(term_agreement, file.path(out_root, "term_agreement.csv"))

print(bin_correlation, n = Inf)
print(comparison |> select(species, site, metric, bin_days, dev_expl, terms),
      n = Inf, width = Inf)
print(term_agreement, n = Inf)
