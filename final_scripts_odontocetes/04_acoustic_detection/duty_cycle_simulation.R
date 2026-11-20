# =============================================================================
# Script  : duty_cycle_simulation.R
# Section : 04_acoustic_detection (Supplementary analysis, duty-cycle effect)
# Purpose : Quantify how much odontocete presence is lost because EI and KGI
#           recorded on a duty cycle. The continuously recorded CI detection
#           log is subsampled with the EI (5 of 6 min) and KGI (10 of
#           13.33 min) duty cycles at many random phase offsets, and the
#           detection events, 5-min presence bins, and presence days that
#           would still have been detected are counted per species.
# Inputs  : data/Antarc_EIE_01_Odontocetes.csv (CI detection log; continuous
#           recording). Required columns: Species Code, Call, Start time,
#           End time (times "m/d/Y H:M:S", UTC).
# Outputs : F_DUTY_SIM (outputs/04_acoustic_detection/duty_cycle_simulation.csv)
#             one row per duty cycle x phase x species with counts and
#             percentages retained
#           outputs/04_acoustic_detection/duty_cycle_simulation_table.csv
#             formatted supplementary table: mean (range) % retained
#           outputs/04_acoustic_detection/duty_cycle_simulation.png / .pdf
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Presence uses echolocation clicks only (PRESENCE_CALLS), and
#             zero-length events are set to 1 s, exactly as in
#             build_presence.R.
#           - An event is retained if any part of it falls in an "on" window.
#             Bins and days are computed from the retained (clipped) parts of
#             events, so a presence day is retained only if a retained bin
#             falls on that day (not from event start dates).
#           - Results support the decision not to use an effort offset in the
#             habitat models: retention is far higher than the on-fraction.
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
N_PHASES             <- 100   # random phase offsets per duty cycle
RANDOM_SEED          <- 20261116
MIN_EVENT_DURATION_S <- 1     # same as build_presence.R
CONTINUOUS_SITE      <- "CI"  # site recorded without duty cycle
BIN_S                <- PRESENCE_BIN_MIN * 60
SIM_TABLE_FILE       <- out_path("04_acoustic_detection",
                                 "duty_cycle_simulation_table.csv")
SIM_FIG_FILE         <- out_path("04_acoustic_detection",
                                 "duty_cycle_simulation")   # .png / .pdf

# -----------------------------------------------------------------------------
# 1. Continuous detection log (clicks only)
# -----------------------------------------------------------------------------
detections <- read_csv(file.path(RAW_DIR, site_info(CONTINUOUS_SITE)$detection_log),
                       col_types = cols(.default = col_character()),
                       show_col_types = FALSE) %>%
  transmute(species = `Species Code`,
            call    = Call,
            start   = as.numeric(mdy_hms(`Start time`, tz = "UTC")),
            end     = as.numeric(mdy_hms(`End time`, tz = "UTC"))) %>%
  filter(call %in% PRESENCE_CALLS, species %in% SPECIES_CODES) %>%
  mutate(end = pmax(end, start + MIN_EVENT_DURATION_S),
         event_id = row_number())

# Duty cycles to simulate, taken from the duty-cycled sites in config.R.
cycles <- SITES %>%
  filter(!is.na(duty_on_min)) %>%
  transmute(cycle = site,
            cycle_label = sprintf("%s (%g of %g min)", site, duty_on_min,
                                  round(duty_cycle_min, 2)),
            on_s = duty_on_min * 60,
            cycle_s = duty_cycle_min * 60,
            on_fraction = duty_on_min / duty_cycle_min)

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
# Bins (integer index of 5-min bins since the epoch) touched by intervals
# [lo, hi): first bin contains lo, last bin is the last one starting before hi.
bins_of <- function(lo, hi) {
  map2(floor(lo / BIN_S), ceiling(hi / BIN_S) - 1, seq)
}

# Presence of each species summarised as events, bins, and days.
summarise_presence <- function(pieces) {
  pieces %>%
    mutate(bin = bins_of(lo, hi)) %>%
    unnest(bin) %>%
    mutate(day = floor(bin * BIN_S / 86400)) %>%
    group_by(species) %>%
    summarise(n_events = n_distinct(event_id),
              n_bins   = n_distinct(bin),
              n_days   = n_distinct(day),
              .groups  = "drop")
}

# Clip every event to the "on" windows [phase + k*cycle, phase + k*cycle + on)
# and return only the parts that were recorded.
clip_to_duty_cycle <- function(events, on_s, cycle_s, phase_s) {
  events %>%
    mutate(k = map2(floor((start - phase_s) / cycle_s),
                    floor((end - phase_s) / cycle_s), seq)) %>%
    unnest(k) %>%
    mutate(lo = pmax(start, phase_s + k * cycle_s),
           hi = pmin(end, phase_s + k * cycle_s + on_s)) %>%
    filter(hi > lo) %>%
    select(species, event_id, lo, hi)
}

# -----------------------------------------------------------------------------
# 3. Reference (continuous) presence
# -----------------------------------------------------------------------------
full_presence <- detections %>%
  transmute(species, event_id, lo = start, hi = end) %>%
  summarise_presence()

# -----------------------------------------------------------------------------
# 4. Simulation over random phase offsets
# -----------------------------------------------------------------------------
# Phase offsets are drawn uniformly over one full cycle because the true
# alignment of the duty cycle relative to any given event is arbitrary.
set.seed(RANDOM_SEED)
sim_grid <- cycles %>%
  mutate(phase_s = map(cycle_s, ~ runif(N_PHASES, 0, .x)),
         phase_id = list(seq_len(N_PHASES))) %>%
  unnest(c(phase_s, phase_id))

sim_results <- sim_grid %>%
  mutate(kept = pmap(list(on_s, cycle_s, phase_s), function(on_s, cycle_s, phase_s) {
    clip_to_duty_cycle(detections, on_s, cycle_s, phase_s) %>%
      summarise_presence()
  })) %>%
  select(cycle, cycle_label, on_fraction, phase_id, phase_s, kept) %>%
  unnest(kept) %>%
  rename(n_events_kept = n_events, n_bins_kept = n_bins, n_days_kept = n_days) %>%
  # Species with no retained events in a phase still need a row (0 kept).
  right_join(crossing(distinct(sim_grid, cycle, cycle_label, phase_id, phase_s),
                      full_presence) %>%
               left_join(distinct(cycles, cycle, on_fraction), by = "cycle"),
             by = c("cycle", "cycle_label", "on_fraction", "phase_id",
                    "phase_s", "species")) %>%
  mutate(across(c(n_events_kept, n_bins_kept, n_days_kept), ~ replace_na(.x, 0)),
         pct_events_kept = 100 * n_events_kept / n_events,
         pct_bins_kept   = 100 * n_bins_kept / n_bins,
         pct_days_kept   = 100 * n_days_kept / n_days,
         species = factor(species, levels = SPECIES_CODES)) %>%
  select(cycle, cycle_label, on_fraction, phase_id, phase_s, species,
         n_events, n_events_kept, pct_events_kept,
         n_bins, n_bins_kept, pct_bins_kept,
         n_days, n_days_kept, pct_days_kept) %>%
  arrange(cycle, species, phase_id)

write_csv(sim_results %>% mutate(across(where(is.double), ~ round(.x, 3))),
          F_DUTY_SIM)

# -----------------------------------------------------------------------------
# 5. Supplementary table: mean (range) across phases
# -----------------------------------------------------------------------------
fmt_mean_range <- function(x) {
  sprintf("%.1f (%.1f–%.1f)", mean(x), min(x), max(x))
}

sim_table <- sim_results %>%
  group_by(species, cycle_label, on_fraction) %>%
  summarise(`Detection events (CI)`     = first(n_events),
            `Presence bins (CI)`        = first(n_bins),
            `Presence days (CI)`        = first(n_days),
            `% events retained`         = fmt_mean_range(pct_events_kept),
            `% 5-min bins retained`     = fmt_mean_range(pct_bins_kept),
            `% presence days retained`  = fmt_mean_range(pct_days_kept),
            .groups = "drop") %>%
  mutate(Species = map_chr(as.character(species), ~ SPECIES[[.x]]$name),
         `Duty cycle` = cycle_label,
         `Recording on-time (%)` = round(100 * on_fraction, 1)) %>%
  arrange(species, `Duty cycle`) %>%
  select(Species, Code = species, `Duty cycle`, `Recording on-time (%)`,
         `Detection events (CI)`, `Presence bins (CI)`, `Presence days (CI)`,
         `% events retained`, `% 5-min bins retained`, `% presence days retained`)

write_excel_csv(sim_table, SIM_TABLE_FILE)
print(sim_table, width = Inf)

# -----------------------------------------------------------------------------
# 6. Supplementary figure: retention by species and duty cycle
# -----------------------------------------------------------------------------
# Points = mean across phases, bars = range; dashed lines mark the on-time
# fraction, i.e., the retention expected if loss were proportional to effort.
fig_data <- sim_results %>%
  pivot_longer(c(pct_events_kept, pct_bins_kept, pct_days_kept),
               names_to = "metric", values_to = "pct") %>%
  mutate(metric = factor(metric,
                         levels = c("pct_events_kept", "pct_bins_kept", "pct_days_kept"),
                         labels = c("Detection events", "5-min presence bins",
                                    "Presence days"))) %>%
  group_by(cycle, cycle_label, on_fraction, species, metric) %>%
  summarise(mean = mean(pct), lo = min(pct), hi = max(pct), .groups = "drop")

cycle_colors <- SITE_COLORS[cycles$cycle] %>% set_names(cycles$cycle_label)

# Flag species that were never detected at the duty-cycled site itself (e.g.,
# sperm whales at KGI): the simulation still shows how that duty cycle WOULD
# perform on CI data, but readers should not mistake it for a site result.
# Open symbols + a short label mark these combinations.
detected_at_site <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) %>%
  pivot_longer(any_of(SPECIES_CODES), names_to = "species", values_to = "minutes") %>%
  group_by(site, species) %>%
  summarise(detected = any(minutes > 0), .groups = "drop")
fig_data <- fig_data %>%
  left_join(detected_at_site, by = c("cycle" = "site", "species")) %>%
  mutate(species = factor(species, levels = intersect(SPECIES_CODES, unique(species))),
         detected = coalesce(detected, FALSE),
         status = if_else(detected, "Detected at this site",
                          "Not detected at this site"))
not_detected <- fig_data %>%
  filter(!detected) %>%
  mutate(label = paste("not detected\nat", cycle))
dodge <- position_dodge(width = 0.5)

p_sim <- ggplot(fig_data, aes(x = species, colour = cycle_label)) +
  geom_hline(data = distinct(fig_data, cycle_label, on_fraction),
             aes(yintercept = 100 * on_fraction, colour = cycle_label),
             linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  geom_linerange(aes(ymin = lo, ymax = hi), position = dodge, linewidth = 0.6) +
  geom_point(aes(y = mean, shape = status), position = dodge, size = 1.8,
             fill = "white", stroke = 0.8) +
  facet_wrap(~ metric, nrow = 1) +
  scale_colour_manual(values = cycle_colors, name = "Simulated duty cycle") +
  scale_shape_manual(values = c("Detected at this site" = 16,
                                "Not detected at this site" = 21),
                     name = NULL) +
  scale_y_continuous(limits = c(60, 100), breaks = seq(60, 100, 10)) +
  # No in-figure caption: the explanation is in the manuscript figure legend
  # (04_acoustic_detection/methods_draft.md).
  labs(x = NULL, y = "Retained relative to continuous recording (%)") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95"))

ggsave(paste0(SIM_FIG_FILE, ".png"), p_sim, width = 7, height = 3.3, dpi = 600)
ggsave(paste0(SIM_FIG_FILE, ".pdf"), p_sim, width = 7, height = 3.3,
       device = cairo_pdf)

message("Wrote ", F_DUTY_SIM, "\nWrote ", SIM_TABLE_FILE,
        "\nWrote ", SIM_FIG_FILE, ".png/.pdf")
