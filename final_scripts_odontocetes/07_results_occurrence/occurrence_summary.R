# =============================================================================
# Script  : occurrence_summary.R
# Section : 07_results_occurrence (Results 3.1, Acoustic occurrence)
# Purpose : Summarise odontocete acoustic occurrence per species x site
#           (main table) and per species x site x month (supplementary
#           table): recording effort, days with detections, and presence time.
# Inputs  : F_PRESENCE_DAILY (outputs/04_acoustic_detection/presence_daily.csv;
#           written by 04_acoustic_detection/build_presence.R)
# Outputs : outputs/07_results_occurrence/occurrence_by_site.csv   (main)
#           outputs/07_results_occurrence/occurrence_by_month.csv  (supplement)
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - All six signal types are reported, including BW58 (too rare to
#             model, reported descriptively).
#           - recording_days = calendar days with recording effort (partial
#             first/last days count as recording days).
#           - recording_hours = hours actually recorded (duty-cycle on-time),
#             i.e., sum of rec_minutes / 60.
#           - pct_time_present = presence time / clock time inside the
#             recording period (= proportion of recorded 5-min bins with
#             presence). Presence minutes are counted in 5-min clock bins, so
#             dividing by duty-cycle on-time would inflate the duty-cycled
#             sites by 1/on-fraction (x1.2 EI, x1.33 KGI); the duty-cycle
#             simulation shows the true loss is much smaller than that.
#             [CONFIRM] denominator choice for "% recorded time with presence".
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

if (!file.exists(F_PRESENCE_DAILY)) {
  stop("Missing ", F_PRESENCE_DAILY,
       "; run 04_acoustic_detection/build_presence.R first.")
}

# -----------------------------------------------------------------------------
# 1. Daily presence in long format, with clock-time effort
# -----------------------------------------------------------------------------
# Clock minutes are recovered from rec_minutes by undoing the duty-cycle
# on-fraction (CI is continuous, on-fraction = 1).
on_fraction <- SITES %>%
  transmute(site, on_fraction = if_else(is.na(duty_on_min), 1,
                                        duty_on_min / duty_cycle_min))

daily <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) %>%
  left_join(on_fraction, by = "site") %>%
  mutate(clock_minutes = rec_minutes / on_fraction,
         month = floor_date(date, "month")) %>%
  pivot_longer(all_of(SPECIES_CODES), names_to = "species",
               values_to = "presence_min")

# -----------------------------------------------------------------------------
# 2. Summary function shared by both tables
# -----------------------------------------------------------------------------
summarise_occurrence <- function(data, ...) {
  data %>%
    filter(rec_minutes > 0) %>%
    group_by(species, site, ...) %>%
    summarise(recording_days   = n(),
              recording_hours  = sum(rec_minutes) / 60,
              detection_days   = sum(presence_min > 0),
              pct_days_with_detections = 100 * detection_days / recording_days,
              presence_hours   = sum(presence_min) / 60,
              pct_time_present = 100 * sum(presence_min) / sum(clock_minutes),
              .groups = "drop") %>%
    mutate(species_name = map_chr(species, ~ SPECIES[[.x]]$name),
           species = factor(species, levels = SPECIES_CODES),
           site = factor(site, levels = SITE_ORDER)) %>%
    relocate(species_name, .after = species) %>%
    mutate(across(c(recording_hours, presence_hours), ~ round(.x, 1)),
           across(c(pct_days_with_detections), ~ round(.x, 1)),
           pct_time_present = round(pct_time_present, 3))
}

# -----------------------------------------------------------------------------
# 3. Main table: species x site
# -----------------------------------------------------------------------------
occ_site <- summarise_occurrence(daily) %>%
  arrange(species, site)
write_csv(occ_site, out_path("07_results_occurrence", "occurrence_by_site.csv"))
print(occ_site, n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# 4. Supplementary table: species x site x calendar month
# -----------------------------------------------------------------------------
# Deployments do not repeat a calendar month within a site (KGI runs Feb 2015
# to Jan 2016), so year-month is used directly.
occ_month <- summarise_occurrence(daily, month) %>%
  mutate(month = format(month, "%Y-%m")) %>%
  arrange(species, site, month)
write_csv(occ_month, out_path("07_results_occurrence", "occurrence_by_month.csv"))

message("Wrote ", out_path("07_results_occurrence", "occurrence_by_site.csv"),
        "\nWrote ", out_path("07_results_occurrence", "occurrence_by_month.csv"))
