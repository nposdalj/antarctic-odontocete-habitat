# =============================================================================
# Script  : occurrence_figures.R
# Section : 07_results_occurrence (Results 3.1, Acoustic occurrence)
# Purpose : Publication figures of odontocete acoustic occurrence:
#           (i)   MAIN TEXT  - weekly presence (days with detections per week)
#                 per species x site, with recording effort per week and
#                 periods without recording shaded; BW58 excluded.
#           (ii)  SUPPLEMENT - weekly presence of the rare BW58 signal.
#           (iii) SUPPLEMENT - per site, stacked daily time series of the
#                 habitat-model predictors (as used in the models) and sea ice
#                 series over the full calendar year(s) of the deployment,
#                 recording period highlighted, daily presence at the bottom.
# Inputs  : F_PRESENCE_DAILY (outputs/04_acoustic_detection/presence_daily.csv;
#             from 04_acoustic_detection/build_presence.R) - required
#           F_ENV_DAILY (outputs/05_environmental_data/env_daily_wide.csv;
#             from 05_environmental_data) - optional; figure (iii) is skipped
#             with a message if it does not exist yet
# Outputs : outputs/07_results_occurrence/
#             fig_weekly_presence.png/.pdf          (i)
#             figS_weekly_presence_BW58.png/.pdf    (ii)
#             figS_env_presence_<site>.png/.pdf     (iii)
# Author  : Natalie Posdaljian
# Date    : 2026-11-16 (revised 2026-11-19: figure (iii) shows the final
#           habitat-model predictors as anomalies by depth layer)
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Weeks start on Monday (UTC). Grey bars show the number of days
#             with recording in each week, so partial-effort weeks (start and
#             end of each deployment) are visible as short grey bars.
#           - A day "with detections" has >= 1 present 5-min bin (clicks only).
#           - PNG at 600 dpi and vector PDF are saved for every figure.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
# Deseasoning settings and functions, so figure (iii) shows the predictors
# exactly as they enter the habitat models (anomalies from the 10-yr
# climatology, log10 lagged chl-a)
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# Script parameters
# -----------------------------------------------------------------------------
FIG_WIDTH_IN <- 7.5   # double-column width
FIG_DPI      <- 600
MAIN_SPECIES <- setdiff(SPECIES_CODES, "BW58")   # BW58 shown in supplement only
RARE_SPECIES <- "BW58"

# Environmental variables for figure (iii), in panel order, with axis labels:
# the habitat-model predictors as they enter the models (revised 2026-11-19),
# plus sea-ice concentration (sea-ice models, raw). "anom." = anomaly from the
# 10-yr climatology; layers: surface 0-100 m, mid-water 100-500 m, deep
# 500-1000 m. Mid-water variables are used by sperm, pilot and killer whales; chl-a is not used.
ENV_VARS <- c(
  ice_conc                 = "Sea ice\n(fraction)",
  temperature_surface_anom = "Temp. anom.,\n0–100 m (°C)",
  temperature_mid_anom     = "Temp. anom.,\n100–500 m (°C)",
  temperature_deep_anom    = "Temp. anom.,\n500–1000 m (°C)",
  mixed_layer_anom         = "MLD anom.\n(m)",
  FSLE                     = "FSLE\n(day⁻¹)",
  fsle_orient              = "FSLV orient.\n(°)",
  EKE_mid                  = "EKE, 100–\n500 m (cm² s⁻²)",
  EKE_deep                 = "EKE, 500–\n1000 m (cm² s⁻²)"
)
ANOM_PANELS <- grep("_anom$", names(ENV_VARS), value = TRUE)

SITE_LABELS <- setNames(sprintf("%s (%s)", SITES$site_name, SITES$site), SITES$site)

save_figure <- function(plot, file_stem, height_in, width_in = FIG_WIDTH_IN) {
  ggsave(out_path("07_results_occurrence", paste0(file_stem, ".png")), plot,
         width = width_in, height = height_in, dpi = FIG_DPI, bg = "white")
  ggsave(out_path("07_results_occurrence", paste0(file_stem, ".pdf")), plot,
         width = width_in, height = height_in, device = cairo_pdf)
  message("Saved ", file_stem, ".png/.pdf")
}

# Wrapped common name + code for panel strips, e.g. "Sperm whale\n(Pm)".
species_label <- function(code) {
  map_chr(code, ~ paste0(str_wrap(SPECIES[[.x]]$name, 14), "\n(", .x, ")"))
}

theme_occurrence <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        strip.background = element_rect(fill = "grey95", colour = "grey60"),
        strip.text.y.left = element_text(angle = 0, hjust = 0),
        strip.text.y.right = element_text(angle = 0, hjust = 0),
        legend.position = "top",
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_blank())

# -----------------------------------------------------------------------------
# 1. Daily presence and weekly summaries
# -----------------------------------------------------------------------------
if (!file.exists(F_PRESENCE_DAILY)) {
  stop("Missing ", F_PRESENCE_DAILY,
       "; run 04_acoustic_detection/build_presence.R first.")
}
daily <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) %>%
  filter(rec_minutes > 0) %>%
  mutate(site = factor(site, levels = SITE_ORDER))

daily_long <- daily %>%
  pivot_longer(all_of(SPECIES_CODES), names_to = "species",
               values_to = "presence_min")

# Weekly counts: days recorded (effort) and days with detections per species.
weekly <- daily_long %>%
  mutate(week_start = floor_date(date, "week", week_start = 1)) %>%
  group_by(site, species, week_start) %>%
  summarise(rec_days = n(), det_days = sum(presence_min > 0),
            # Same unit as the habitat-model response: mean presence minutes
            # per recorded day (positive 5-min bins), here over each week
            min_per_day = sum(presence_min) / rec_days, .groups = "drop")

# Time without recording, shaded in the weekly plots: before the first
# deployment, between deployments, and after the last one. Edges are snapped
# outward to whole weeks, so a partial-effort week sits on the shading and
# the unrecorded part of that week reads as "no recording".
week_floor <- function(x) floor_date(as.Date(x), "week", week_start = 1)
x_range <- c(floor_date(week_floor(min(SITES$rec_start)), "month"),
             ceiling_date(week_floor(max(SITES$rec_end)) + 7, "month"))
no_recording <- SITES %>%
  arrange(rec_start) %>%
  transmute(xmin = c(x_range[1], week_floor(head(rec_end, -1))),
            xmax = week_floor(rec_start) + 7) %>%
  bind_rows(tibble(xmin = week_floor(max(SITES$rec_end)), xmax = x_range[2]))

# -----------------------------------------------------------------------------
# 2. Weekly presence plot (used for figures i and ii)
# -----------------------------------------------------------------------------
# Bars are centred on each week (week_start + 3.5 d). Grey = days recorded,
# coloured = days with detections, so partial-effort weeks and zero-detection
# weeks are both distinguishable from no-recording periods (dark shading).
plot_weekly <- function(species_codes) {
  wk <- weekly %>%
    filter(species %in% species_codes) %>%
    mutate(species = factor(species, levels = species_codes,
                            labels = species_label(species_codes)),
           week_mid = week_start + 3.5)

  gap_labels <- no_recording %>%
    filter(xmax - xmin > 60) %>%
    transmute(x = xmin + (xmax - xmin) / 2,
              species = factor(levels(wk$species)[1], levels = levels(wk$species)))

  ggplot(wk, aes(x = week_mid)) +
    geom_rect(data = no_recording, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey72") +
    geom_col(aes(y = rec_days, fill = "Days recorded"), width = 7) +
    geom_col(data = filter(wk, det_days > 0),
             aes(y = det_days, fill = as.character(site)), width = 6,
             show.legend = TRUE) +
    # Label long gaps once, in the top panel only.
    geom_text(data = gap_labels, inherit.aes = FALSE,
              aes(x = x, y = 3.5, label = "No recording"),
              size = 2.4, colour = "white", fontface = "bold") +
    facet_grid(species ~ ., switch = "y") +
    scale_fill_manual(values = c(`Days recorded` = "grey90", SITE_COLORS),
                      limits = c(SITE_ORDER, "Days recorded"),
                      labels = c(SITE_LABELS[SITE_ORDER], "Days recorded")) +
    scale_x_date(limits = x_range, expand = c(0, 0),
                 date_breaks = "3 months", date_labels = "%b\n%Y",
                 date_minor_breaks = "1 month") +
    scale_y_continuous(limits = c(0, 7), breaks = c(0, 7), expand = c(0, 0)) +
    labs(x = NULL, y = "Days with detections per week") +
    theme_occurrence +
    theme(strip.placement = "outside",
          panel.spacing.y = unit(0.3, "cm"))
}

# Weekly mean presence minutes per recorded day (decision N. Posdaljian
# 2026-11-25: show presence in the units of the habitat-model response instead
# of days with detections). Free y axis per species; weeks recorded without
# detections show a grey baseline tick, unrecorded periods are shaded.
plot_weekly_minutes <- function(species_codes) {
  wk <- weekly %>%
    filter(species %in% species_codes) %>%
    mutate(species = factor(species, levels = species_codes,
                            labels = species_label(species_codes)),
           week_mid = week_start + 3.5)
  ggplot(wk, aes(x = week_mid)) +
    geom_rect(data = no_recording, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey72") +
    geom_segment(aes(xend = week_mid, y = 0, yend = 0), colour = "grey55",
                 linewidth = 1.2) +
    geom_col(data = filter(wk, min_per_day > 0),
             aes(y = min_per_day, fill = as.character(site)), width = 6) +
    facet_grid(species ~ ., switch = "y", scales = "free_y") +
    scale_fill_manual(values = SITE_COLORS, limits = SITE_ORDER,
                      labels = SITE_LABELS[SITE_ORDER], name = NULL) +
    scale_x_date(limits = x_range, expand = c(0, 0),
                 date_breaks = "3 months", date_labels = "%b
%Y",
                 date_minor_breaks = "1 month") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Presence (minutes per day, weekly mean)") +
    theme_occurrence +
    theme(strip.placement = "outside",
          panel.spacing.y = unit(0.3, "cm"))
}

# (i) Main text: all species, BW58 included as the last panel
#     (decision N. Posdaljian 2026-11-20; the separate BW58 supplementary
#     figure is no longer produced).
all_species <- c(MAIN_SPECIES, RARE_SPECIES)
save_figure(plot_weekly_minutes(all_species), "fig_weekly_presence",
            height_in = 1.0 + 0.95 * length(all_species))

# -----------------------------------------------------------------------------
# 3. (iii) Supplement: environmental time series + presence per site
# -----------------------------------------------------------------------------
# Window: 1 Jan of the deployment start year to 31 Dec of that year, or to the
# end of the last recording month if the deployment runs into the next year
# (KGI: Jan 2015 - Jan 2016).
plot_env_site <- function(site_code, env) {
  info <- site_info(site_code)
  win_start <- floor_date(as.Date(info$rec_start), "year")
  win_end   <- max(ceiling_date(win_start, "year"),
                   ceiling_date(as.Date(info$rec_end), "month")) - 1
  x_lim <- c(win_start, win_end + 1)
  rec_period <- tibble(xmin = as.Date(info$rec_start), xmax = as.Date(info$rec_end))

  env_long <- env %>%
    filter(site == site_code, date >= win_start, date <= win_end) %>%
    select(date, any_of(names(ENV_VARS))) %>%
    pivot_longer(-date, names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = names(ENV_VARS),
                             labels = ENV_VARS))

  if (nrow(env_long) == 0) {
    message("No environmental data for ", site_code, " in ", win_start, " to ",
            win_end, "; skipping.")
    return(NULL)
  }

  zero_lines <- tibble(variable = factor(ENV_VARS[ANOM_PANELS], levels = ENV_VARS))

  p_env <- ggplot(env_long, aes(date, value)) +
    geom_rect(data = rec_period, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = SITE_COLORS[[site_code]], alpha = 0.10) +
    geom_hline(data = zero_lines, aes(yintercept = 0), colour = "grey60",
               linewidth = 0.3) +
    geom_line(linewidth = 0.35, colour = "grey15", na.rm = TRUE) +
    facet_grid(variable ~ ., scales = "free_y", switch = "y") +
    scale_x_date(limits = x_lim, expand = c(0, 0), date_breaks = "1 month",
                 date_labels = "%b") +
    labs(x = NULL, y = NULL,
         title = paste0(SITE_LABELS[[site_code]], ": recording period ",
                        format(as.Date(info$rec_start), "%d %b %Y"), " – ",
                        format(as.Date(info$rec_end), "%d %b %Y"), " (shaded)")) +
    theme_occurrence +
    theme(strip.placement = "outside",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.spacing.y = unit(0.1, "cm"),
          plot.title = element_text(size = 8, face = "bold"))

  # Presence strip: one tile per detection day, species detected at the site.
  pres <- daily_long %>%
    filter(site == site_code, presence_min > 0) %>%
    mutate(species = factor(species, levels = rev(SPECIES_CODES)))

  p_pres <- ggplot(pres, aes(date, species)) +
    geom_rect(data = rec_period, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = SITE_COLORS[[site_code]], alpha = 0.10) +
    geom_tile(fill = SITE_COLORS[[site_code]], height = 0.7, width = 1) +
    scale_y_discrete(drop = FALSE) +
    scale_x_date(limits = x_lim, expand = c(0, 0), date_breaks = "1 month",
                 date_labels = "%b\n%Y") +
    labs(x = NULL, y = "Detection\ndays") +
    theme_occurrence +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

  n_vars <- n_distinct(env_long$variable)
  p_env / p_pres + plot_layout(heights = c(n_vars, 0.3 * length(SPECIES_CODES)))
}

if (!file.exists(F_ENV_DAILY)) {
  message("Skipping figure (iii): ", F_ENV_DAILY, " not found. Run the ",
          "05_environmental_data scripts first, then re-run this script.")
} else {
  env <- read_csv(F_ENV_DAILY, show_col_types = FALSE) |>
    mutate(date = as.Date(date))
  # Same transformation as 01_stage1_automated_selection.R: log10 lagged chl-a,
  # then anomalies from the 10-yr periodic-spline climatology.
  if (CHL_LOG10) env <- mutate(env, across(any_of(LAG_COLUMNS), log10))
  anom_sources <- intersect(sub("_anom$", "", ANOM_PANELS), names(env))
  env <- add_anomalies(env, build_climatology(env, anom_sources), anom_sources)
  if (CHL_STANDARDIZE) env <- standardize_anomalies(env, env, intersect(LAG_COLUMNS, anom_sources))
  missing_vars <- setdiff(names(ENV_VARS), names(env))
  if (length(missing_vars) > 0) {
    message("Figure (iii): variables not in ", basename(F_ENV_DAILY),
            " and not plotted: ", paste(missing_vars, collapse = ", "))
  }
  for (site_code in SITE_ORDER) {
    p <- plot_env_site(site_code, env)
    if (!is.null(p)) {
      save_figure(p, paste0("figS_env_presence_", site_code),
                  height_in = 1.2 + 0.75 * sum(names(ENV_VARS) %in% names(env)) +
                    0.15 * length(SPECIES_CODES))
    }
  }
}
