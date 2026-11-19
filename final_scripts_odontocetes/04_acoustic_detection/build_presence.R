# =============================================================================
# Script  : build_presence.R
# Section : 04_acoustic_detection (Methods 2.4, Acoustic detection)
# Purpose : Convert the analyst-verified odontocete detection logs from the
#           three HARP sites into (1) a 5-min presence grid and (2) daily
#           presence minutes with daily recording effort. These two tables
#           are the response data for all occurrence summaries (07) and
#           habitat models.
# Inputs  : data/Antarc_EI_01_Odontocetes.csv   (EI detection log)
#           data/Antarc_SSI_01_Odontocetes.csv  (KGI detection log)
#           data/Antarc_EIE_01_Odontocetes.csv  (CI detection log)
#           Required columns: Species Code, Call, Start time, End time
#           (times "m/d/Y H:M:S", UTC). File names are set in config.R.
# Outputs : F_PRESENCE_5MIN  (outputs/04_acoustic_detection/presence_5min.csv)
#             bin_start, site, recorded, <species> (logical)
#           F_PRESENCE_DAILY (outputs/04_acoustic_detection/presence_daily.csv)
#             date, site, rec_minutes, <species> (presence minutes per day)
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Presence is defined by echolocation clicks only
#             (PRESENCE_CALLS in config.R); sperm whale slow clicks and
#             killer whale HFM whistles are excluded.
#           - A 5-min bin is "present" for a species if any detection event
#             overlaps it. Daily presence minutes = present bins x 5.
#           - rec_minutes = minutes actually recorded on that day: clock time
#             inside the recording period (partial first/last days) times the
#             duty-cycle on-fraction (duty_on_min / duty_cycle_min; CI = 1).
#           - Bins are aligned to the UTC clock (00:00, 00:05, ...). Because
#             the duty-cycle phase is not stored in the logs, a bin is flagged
#             "recorded" if it overlaps the recording period, regardless of
#             whether the recorder was on during that particular bin.
#           - Replaces step 4 of "Data Visualizations/TimeseriesPlots.R",
#             which used all call types, ignored the duty cycle, and wrote
#             the minutes table to the "_5min_bins" file name by mistake.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# -----------------------------------------------------------------------------
# Script parameters
# -----------------------------------------------------------------------------
# Detection events with identical start and end times (single clicks logged
# at one time stamp) are given this duration so that they occupy a bin.
MIN_EVENT_DURATION_S <- 1

# Calendar days with less than this much clock time inside the recording
# period are dropped. This removes 2016-02-03 at CI, where the configured
# start time (23:59:59) leaves only 1 s of recording on that day.
# [CONFIRM] CI start of 2016-02-03 23:59:59 vs. first full day 2016-02-04.
MIN_DAY_CLOCK_MIN <- 1

BIN_S <- PRESENCE_BIN_MIN * 60

# -----------------------------------------------------------------------------
# 1. Read and clean the detection logs
# -----------------------------------------------------------------------------
# The three logs have different numbers of extra columns (the EI log also has
# click-parameter columns), so columns are selected by name.
read_detection_log <- function(site_code) {
  info <- site_info(site_code)
  read_csv(file.path(RAW_DIR, info$detection_log),
           col_types = cols(.default = col_character()),
           show_col_types = FALSE) %>%
    transmute(
      site    = site_code,
      species = `Species Code`,
      call    = Call,
      start   = mdy_hms(`Start time`, tz = "UTC"),
      end     = mdy_hms(`End time`, tz = "UTC")
    )
}

detections_all <- map_dfr(SITE_ORDER, read_detection_log)

# Stop rather than silently drop events whose times could not be parsed.
n_bad_time <- sum(is.na(detections_all$start) | is.na(detections_all$end))
if (n_bad_time > 0) stop(n_bad_time, " detection events have unparseable times.")

# Keep echolocation clicks only (presence definition), and only species
# codes defined in config.R.
detections <- detections_all %>%
  filter(call %in% PRESENCE_CALLS, species %in% SPECIES_CODES) %>%
  # Zero-length events -> MIN_EVENT_DURATION_S so they overlap one bin.
  mutate(end = pmax(end, start + MIN_EVENT_DURATION_S))

message("Detection events kept (clicks only): ", nrow(detections),
        " of ", nrow(detections_all))
print(count(detections_all, site, species, call) %>%
        mutate(kept = call %in% PRESENCE_CALLS), n = Inf)

# Every event must lie inside its site's recording period; anything outside
# would indicate a time-zone or log error.
outside <- detections %>%
  left_join(select(SITES, site, rec_start, rec_end), by = "site") %>%
  filter(start < rec_start | end > rec_end)
if (nrow(outside) > 0) {
  warning(nrow(outside), " detection events extend outside the recording period.")
}

# -----------------------------------------------------------------------------
# 2. Recording effort per site
# -----------------------------------------------------------------------------
# Fraction of clock time the recorder was on (CI recorded continuously).
duty_on_fraction <- function(info) {
  if (is.na(info$duty_on_min)) 1 else info$duty_on_min / info$duty_cycle_min
}

# Clock minutes of [period_start, period_end) that fall inside [a, b).
overlap_min <- function(a, b, period_start, period_end) {
  pmax(0, as.numeric(difftime(pmin(b, period_end), pmax(a, period_start),
                              units = "mins")))
}

# Daily effort: one row per calendar day touched by the recording period.
daily_effort <- map_dfr(SITE_ORDER, function(site_code) {
  info <- site_info(site_code)
  day_seq <- seq(as.Date(info$rec_start), as.Date(info$rec_end), by = "day")
  day_start <- as.POSIXct(day_seq, tz = "UTC")
  clock_min <- overlap_min(day_start, day_start + days(1),
                           info$rec_start, info$rec_end)
  tibble(site        = site_code,
         date        = day_seq,
         clock_min   = clock_min,
         rec_minutes = clock_min * duty_on_fraction(info))
}) %>%
  filter(clock_min >= MIN_DAY_CLOCK_MIN)

# -----------------------------------------------------------------------------
# 3. 5-min bin grid
# -----------------------------------------------------------------------------
# The grid covers whole calendar days of every retained recording day so the
# file is regular; bins outside the recording period are recorded = FALSE.
bin_grid <- daily_effort %>%
  select(site, date) %>%
  left_join(select(SITES, site, rec_start, rec_end), by = "site") %>%
  mutate(bin_start = map(date, ~ as.POSIXct(.x, tz = "UTC") +
                           seq(0, 86400 - BIN_S, by = BIN_S))) %>%
  unnest(bin_start) %>%
  mutate(recorded = bin_start + BIN_S > rec_start & bin_start < rec_end) %>%
  select(site, date, bin_start, recorded)

# -----------------------------------------------------------------------------
# 4. Assign detections to bins
# -----------------------------------------------------------------------------
# An event spanning [start, end) occupies every bin from the one containing
# its start to the one containing its last second. Working in integer bin
# indices (seconds since epoch / BIN_S) keeps this exact and fast.
present_bins <- detections %>%
  mutate(first_bin = floor(as.numeric(start) / BIN_S),
         last_bin  = floor((as.numeric(end) - 1) / BIN_S)) %>%
  mutate(bin_index = map2(first_bin, last_bin, seq)) %>%
  select(site, species, bin_index) %>%
  unnest(bin_index) %>%
  distinct() %>%
  mutate(bin_start = as.POSIXct(bin_index * BIN_S, origin = "1970-01-01",
                                tz = "UTC"),
         present = TRUE) %>%
  select(site, species, bin_start, present) %>%
  pivot_wider(names_from = species, values_from = present)

# Species never detected anywhere still get a column (all FALSE).
for (sp in setdiff(SPECIES_CODES, names(present_bins))) present_bins[[sp]] <- NA

presence_5min <- bin_grid %>%
  left_join(present_bins, by = c("site", "bin_start")) %>%
  mutate(across(all_of(SPECIES_CODES), ~ replace_na(.x, FALSE))) %>%
  mutate(site = factor(site, levels = SITE_ORDER)) %>%
  arrange(site, bin_start)

# Every detected bin must land on the grid (i.e., inside a recording day).
n_lost <- nrow(semi_join(present_bins, bin_grid, by = c("site", "bin_start")))
if (n_lost != nrow(present_bins)) {
  warning(nrow(present_bins) - n_lost, " present bins fall outside the bin grid.")
}

# -----------------------------------------------------------------------------
# 5. Daily presence minutes
# -----------------------------------------------------------------------------
presence_daily <- presence_5min %>%
  group_by(site, date) %>%
  summarise(across(all_of(SPECIES_CODES), ~ sum(.x) * PRESENCE_BIN_MIN),
            .groups = "drop") %>%
  left_join(mutate(daily_effort, site = factor(site, levels = SITE_ORDER)),
            by = c("site", "date")) %>%
  select(date, site, rec_minutes, all_of(SPECIES_CODES)) %>%
  arrange(site, date)

# -----------------------------------------------------------------------------
# 6. Write outputs (data contract in config.R)
# -----------------------------------------------------------------------------
presence_5min %>%
  transmute(bin_start = format(bin_start, "%Y-%m-%d %H:%M:%S"), site, recorded,
            across(all_of(SPECIES_CODES))) %>%
  write_csv(F_PRESENCE_5MIN)

presence_daily %>%
  mutate(rec_minutes = round(rec_minutes, 2)) %>%
  write_csv(F_PRESENCE_DAILY)

# Console summary: recording days and detection days per species x site.
presence_daily %>%
  group_by(site) %>%
  summarise(rec_days = n(), rec_hours = round(sum(rec_minutes) / 60, 1),
            across(all_of(SPECIES_CODES), ~ sum(.x > 0)), .groups = "drop") %>%
  print()

message("Wrote ", F_PRESENCE_5MIN, "\nWrote ", F_PRESENCE_DAILY)
