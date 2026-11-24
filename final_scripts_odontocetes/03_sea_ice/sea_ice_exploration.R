# =============================================================================
# Script  : sea_ice_exploration.R
# Section : 2.3 / 3.2 Sea ice - exploratory analysis (descriptive)
# Purpose : Describe sea-ice conditions at EI, KGI and CI over the full
#           10-year record (+/-20 km box, Copernicus reanalysis) to put the
#           acoustic recording years and the sea-ice results in context:
#           (1) seasonal cycle by site; (2) interannual variability and how
#           unusual each recording year was; (3) ice phenology (advance,
#           retreat, duration, ice days, maximum, mean, cumulative exposure);
#           (4) timing of acoustic detections relative to the phase of the
#           seasonal ice cycle, by species and site.
#           Descriptive only: the habitat and sea-ice models are not changed,
#           and no causal claims are made.
# Inputs  : F_ENV_DAILY (daily ice_conc from AMSR2 ASI, 2013-2021; fraction of the box covered by ice;
#           open water counted as 0), F_PRESENCE_DAILY
# Outputs : outputs/03_sea_ice/
#             ice_phenology_by_year.csv   one row per site x ice year x threshold
#             ice_recording_year_rank.csv recording year vs the site's other years
#             ice_site_year_anova.csv     site vs year variance of annual metrics
#             ice_site_synchrony.csv      Spearman r of annual metrics between sites
#             ice_phase_occurrence.csv    detection days by ice phase, species, site
#             ice_seasonal_climatology.png/.pdf   (figure a)
#             ice_year_doy_heatmap.png/.pdf       (figure b)
#             ice_phenology_metrics.png/.pdf      (figure c)
#             ice_occurrence_timeline.png/.pdf    (figure d)
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Ice year = 1 March to the end of February (starts near the summer
#           ice minimum), labelled by its starting calendar year; only complete
#           ice years are used. Phenology follows Stammerjohn et al. (2008,
#           JGR 113, C03S90) [CONFIRM]: advance = first day of the ice year
#           starting a run of >= 5 consecutive days with concentration >=
#           threshold; retreat = last day of the last such run; duration =
#           retreat - advance + 1; ice days = days >= threshold. Because the
#           box mean includes open water, a 15% box mean is a moderately
#           ice-covered box; 5% is also reported as a lower threshold.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })

THRESHOLDS <- c(0.15, 0.05)
RUN_DAYS <- 5
PHASE_WINDOW <- 30   # days after advance / before retreat counted as transition phases
REC_YEAR <- c(EI = 2014, KGI = 2015, CI = 2016)
SECTION <- "03_sea_ice"
save_fig <- function(p, stem, w, h) {
  for (ext in c("png", "pdf")) ggsave(out_path(SECTION, paste0(stem, ".", ext)), p,
                                      width = w, height = h, dpi = 300)
}
site_f <- \(x) factor(x, SITE_ORDER)

# Short AMSR2 data gaps (missing swaths; 1-4 days in May 2013, Nov 2018 and
# Sep 2019) are filled by linear interpolation so that whole ice years are
# not lost; longer gaps (and the period before AMSR2) stay missing.
MAX_GAP_DAYS <- 5
fill_short_gaps <- function(x, max_gap) {
  na <- is.na(x)
  if (!any(na) || all(na)) return(x)
  filled <- approx(which(!na), x[!na], xout = seq_along(x), rule = 1)$y
  r <- rle(na); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  for (k in which(r$values & r$lengths > max_gap)) filled[starts[k]:ends[k]] <- NA
  filled
}

# -----------------------------------------------------------------------------
# 1. Daily ice in ice years
# -----------------------------------------------------------------------------
ice <- read_csv(F_ENV_DAILY, show_col_types = FALSE) |>
  transmute(site, date = as.Date(date), sic = ice_conc) |>
  group_by(site) |> arrange(date, .by_group = TRUE) |>
  mutate(sic = fill_short_gaps(sic, MAX_GAP_DAYS)) |> ungroup() |>
  filter(!is.na(sic)) |>   # AMSR2 ice starts in 2013 (July 2012)
  mutate(ice_year = as.integer(format(date, "%Y")) - (as.integer(format(date, "%m")) < 3),
         day = as.integer(date - as.Date(paste0(ice_year, "-03-01"))) + 1) |>   # 1 = 1 March
  group_by(site, ice_year) |> filter(n() >= 365) |> ungroup() |>
  filter(day <= 365) |>   # drop 29 Feb overflow so every ice year has 365 days
  mutate(site = site_f(site), recording = ice_year == REC_YEAR[as.character(site)])
day_breaks <- c(1, 32, 93, 154, 215, 276, 337)   # 1 Mar, Apr, Jun, Aug, Oct, Dec, Feb
day_labels <- c("Mar", "Apr", "Jun", "Aug", "Oct", "Dec", "Feb")

# -----------------------------------------------------------------------------
# 2. Phenology per site x ice year x threshold
# -----------------------------------------------------------------------------
runs_above <- function(x, thr) {   # start/end of runs of >= RUN_DAYS days >= thr
  r <- rle(x >= thr)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  keep <- r$values & r$lengths >= RUN_DAYS
  list(start = starts[keep], end = ends[keep])
}
phenology <- ice |>
  group_by(site, ice_year, recording) |>
  reframe(map_dfr(THRESHOLDS, function(thr) {
    r <- runs_above(sic, thr)
    tibble(threshold_pct = 100 * thr,
           advance_day = if (length(r$start)) min(r$start) else NA_integer_,
           retreat_day = if (length(r$end)) max(r$end) else NA_integer_,
           ice_days = sum(sic >= thr),
           max_pct = 100 * max(sic), mean_pct = 100 * mean(sic),
           cumulative_ice_days = sum(sic))     # sum of daily fraction = full-ice-day equivalents
  })) |>
  mutate(duration_days = retreat_day - advance_day + 1,
         advance_date = format(as.Date("2001-03-01") + advance_day - 1, "%d %b"),
         retreat_date = format(as.Date("2001-03-01") + retreat_day - 1, "%d %b"))
write_csv(mutate(phenology, across(where(is.double), \(x) round(x, 2))),
          out_path(SECTION, "ice_phenology_by_year.csv"))

metrics <- c(advance_day = "Advance (day of ice year)", retreat_day = "Retreat (day of ice year)",
             duration_days = "Ice-season duration (days)", ice_days = "Ice days (>= threshold)",
             max_pct = "Maximum concentration (%)", mean_pct = "Mean concentration (%)",
             cumulative_ice_days = "Cumulative ice exposure (days)")
ph15 <- filter(phenology, threshold_pct == 15)
ph_long <- phenology |>
  pivot_longer(all_of(names(metrics)), names_to = "metric") |>
  filter(threshold_pct == 15 | metric %in% c("advance_day", "retreat_day", "duration_days", "ice_days"))

# How unusual was each recording year? Rank among the site's ice years.
rank_tbl <- ph_long |>
  group_by(site, threshold_pct, metric) |>
  summarise(n_years = sum(!is.na(value)),
            n_years_no_ice_season = sum(is.na(value)),
            recording_value = round(value[recording], 1),
            median_other = round(median(value[!recording], na.rm = TRUE), 1),
            min_other = round(suppressWarnings(min(value[!recording], na.rm = TRUE)), 1),
            max_other = round(suppressWarnings(max(value[!recording], na.rm = TRUE)), 1),
            percentile = round(100 * mean(value[!recording] <= value[recording], na.rm = TRUE)),
            .groups = "drop")
write_csv(rank_tbl, out_path(SECTION, "ice_recording_year_rank.csv"))

# Spatial vs interannual: two-way ANOVA (site + year) on annual metrics over
# the ice years common to all three sites; share of the sum of squares.
common_years <- Reduce(intersect, split(ph15$ice_year, ph15$site))
anova_tbl <- map_dfr(c("ice_days", "max_pct", "mean_pct", "cumulative_ice_days"), function(m) {
  d <- filter(ph15, ice_year %in% common_years) |> mutate(year = factor(ice_year))
  a <- anova(lm(reformulate(c("site", "year"), m), data = d))
  ss <- setNames(a$`Sum Sq`, rownames(a))
  tibble(metric = m, n_years = length(common_years),
         pct_site = round(100 * ss[["site"]] / sum(ss), 1),
         pct_year = round(100 * ss[["year"]] / sum(ss), 1),
         pct_residual = round(100 * ss[["Residuals"]] / sum(ss), 1),
         p_site = signif(a["site", "Pr(>F)"], 2), p_year = signif(a["year", "Pr(>F)"], 2))
})
write_csv(anova_tbl, out_path(SECTION, "ice_site_year_anova.csv"))

# Synchrony: do the same years have more ice at all sites?
synchrony <- map_dfr(c("cumulative_ice_days", "max_pct"), function(m) {
  w <- ph15 |> filter(ice_year %in% common_years) |> select(site, ice_year, all_of(m)) |>
    pivot_wider(names_from = site, values_from = all_of(m))
  combn(SITE_ORDER, 2, simplify = FALSE) |> map_dfr(\(p) tibble(
    metric = m, pair = paste(p, collapse = "-"),
    spearman_r = round(cor(w[[p[1]]], w[[p[2]]], method = "spearman"), 2)))
})
write_csv(synchrony, out_path(SECTION, "ice_site_synchrony.csv"))

# -----------------------------------------------------------------------------
# 3. Acoustic detections by phase of the ice cycle (recording years; 15%)
# -----------------------------------------------------------------------------
phase_of <- function(day, adv, ret) {
  case_when(is.na(adv) ~ "no ice season",
            day < adv ~ "before advance",
            day < adv + PHASE_WINDOW ~ "advance (first 30 d)",
            day <= ret - PHASE_WINDOW ~ "ice season",
            day <= ret ~ "retreat (last 30 d)",
            TRUE ~ "after retreat")
}
phase_levels <- c("before advance", "advance (first 30 d)", "ice season",
                  "retreat (last 30 d)", "after retreat", "no ice season")
presence <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) |>
  filter(rec_minutes > 0) |> mutate(date = as.Date(date), site = site_f(site))
rec_days <- presence |>
  inner_join(ice, by = c("site", "date")) |>
  left_join(select(ph15, site, ice_year, advance_day, retreat_day), by = c("site", "ice_year")) |>
  mutate(phase = factor(phase_of(day, advance_day, retreat_day), phase_levels))
phase_tbl <- rec_days |>
  pivot_longer(any_of(c(MODELED_SPECIES, "BW58")), names_to = "species", values_to = "minutes") |>
  group_by(site, species, phase) |>
  summarise(recording_days = n(), detection_days = sum(minutes > 0),
            pct_days_detected = round(100 * mean(minutes > 0), 1),
            mean_ice_pct = round(100 * mean(sic), 1), .groups = "drop")
write_csv(phase_tbl, out_path(SECTION, "ice_phase_occurrence.csv"))

# -----------------------------------------------------------------------------
# 4. Figures
# -----------------------------------------------------------------------------
theme_ice <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "bottom")

# (b) Year x day heatmap, recording year outlined
p_b <- ggplot(ice, aes(day, factor(ice_year), fill = 100 * sic)) +
  geom_raster() +
  geom_tile(data = distinct(filter(ice, recording), site, ice_year), inherit.aes = FALSE,
            aes(x = 183, y = factor(ice_year)), width = 365, height = 1,
            fill = NA, colour = "black", linewidth = 0.5) +
  facet_wrap(~ site, nrow = 1, scales = "free_y") +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Sea-ice concentration (%)",
                       limits = c(0, 100)) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels, expand = c(0, 0)) +
  labs(x = NULL, y = "Ice year (starting 1 March)") +
  theme_ice + theme(legend.key.width = unit(1.2, "cm"))
save_fig(p_b, "ice_year_doy_heatmap", 8, 4.5)

# (c) Phenology metrics by site and year (15% threshold), recording year highlighted
p_c <- ph_long |>
  filter(threshold_pct == 15) |>
  mutate(metric = factor(metrics[metric], metrics)) |>
  ggplot(aes(site, value, colour = site)) +
  geom_jitter(data = \(d) filter(d, !recording), width = 0.12, height = 0, size = 1.3,
              alpha = 0.6, na.rm = TRUE) +
  geom_point(data = \(d) filter(d, recording), shape = 23, size = 3, fill = "white",
             stroke = 1.1, na.rm = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = SITE_COLORS, guide = "none") +
  labs(x = NULL, y = NULL,
       subtitle = "Threshold 15%. Circles: other ice years; diamonds: recording year. Missing = no ice season") +
  theme_ice
save_fig(p_c, "ice_phenology_metrics", 8, 4.6)

# (a+d) Combined figure (decision N. Posdaljian 2026-11-20): per site, top =
#     daily sea ice over the ice year (10-90% range and median of all ice
#     years, recording year as a solid line, days without acoustic recording
#     in grey, advance/retreat of the recording year as vertical lines);
#     bottom = days with detections of each species on the same axis.
#     The axis starts on 1 February so that February recordings before the
#     ice year (KGI 2015, CI 2016) are shown; day 1 = 1 March.
FEB_START <- -27
rec_start <- \(st) as.Date(paste0(REC_YEAR[[as.character(st)]], "-03-01"))
to_day <- function(st, date) as.integer(date - rec_start(st)) + 1
feb_prev <- ice |> filter(day >= 365 + FEB_START) |> mutate(day = day - 365)   # Feb of previous ice year
clim2 <- bind_rows(ice, feb_prev) |>
  group_by(site, day) |>
  summarise(med = median(sic), lo = quantile(sic, 0.1), hi = quantile(sic, 0.9), .groups = "drop")
rec_year_ice <- map_dfr(SITE_ORDER, function(st) {
  d <- ice |> filter(site == st, ice_year %in% (REC_YEAR[[st]] + c(-1, 0)))
  d |> mutate(day = to_day(st, date)) |> filter(day >= FEB_START, day <= 365)
})
rec_days2 <- rec_days |> mutate(day = map2_int(site, date, to_day))
not_rec2 <- map_dfr(SITE_ORDER, function(st) {
  recd <- rec_days2$day[rec_days2$site == st]
  tibble(day = FEB_START:365, rec = (FEB_START:365) %in% recd) |>
    mutate(run = cumsum(rec != lag(rec, default = first(rec)))) |>
    filter(!rec) |> group_by(run) |>
    summarise(xmin = min(day) - 0.5, xmax = max(day) + 0.5, .groups = "drop") |>
    mutate(site = st)
}) |> mutate(site = site_f(site))
marks2 <- ph15 |> filter(recording) |>
  transmute(site, advance = advance_day, retreat = retreat_day) |>
  pivot_longer(c(advance, retreat), names_to = "event", values_to = "day") |>
  filter(!is.na(day))
det2 <- rec_days2 |>
  pivot_longer(any_of(c(MODELED_SPECIES, "BW58")), names_to = "species", values_to = "minutes") |>
  filter(minutes > 0) |>
  mutate(species = factor(species, rev(c(MODELED_SPECIES, "BW58"))))
x_axis <- scale_x_continuous(limits = c(FEB_START - 0.5, 365.5), expand = c(0, 0),
                             breaks = c(FEB_START, day_breaks), labels = c("Feb", day_labels))
p_top <- ggplot(clim2, aes(day)) +
  geom_rect(data = not_rec2, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey88") +
  geom_ribbon(aes(ymin = 100 * lo, ymax = 100 * hi, fill = site), alpha = 0.25) +
  geom_line(aes(y = 100 * med, colour = site), linewidth = 0.4, linetype = "dashed") +
  geom_line(data = rec_year_ice, aes(y = 100 * sic, colour = site), linewidth = 0.55) +
  geom_vline(data = marks2, aes(xintercept = day, linetype = event), linewidth = 0.4) +
  facet_wrap(~ site, nrow = 1) +
  scale_colour_manual(values = SITE_COLORS, guide = "none") +
  scale_fill_manual(values = SITE_COLORS, guide = "none") +
  scale_linetype_manual(values = c(advance = "solid", retreat = "dashed"), name = "Ice (15%)") +
  x_axis + labs(x = NULL, y = "Sea-ice concentration (%)") + theme_ice +
  theme(axis.text.x = element_blank(), legend.position = "top")
p_bot <- ggplot(det2, aes(day, species)) +
  geom_rect(data = not_rec2, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey88") +
  geom_tile(aes(fill = site), height = 0.7, width = 1) +
  geom_vline(data = marks2, aes(xintercept = day, linetype = event), linewidth = 0.4,
             show.legend = FALSE) +
  facet_wrap(~ site, nrow = 1) +
  scale_fill_manual(values = SITE_COLORS, guide = "none") +
  scale_linetype_manual(values = c(advance = "solid", retreat = "dashed")) +
  scale_y_discrete(drop = FALSE) +
  x_axis + labs(x = NULL, y = "Detection days") + theme_ice +
  theme(strip.text = element_blank())
save_fig(p_top / p_bot + plot_layout(heights = c(1, 1.2)), "ice_climatology_detections", 10, 5)

print(anova_tbl, width = Inf); print(synchrony)
print(filter(rank_tbl, threshold_pct == 15) |> select(-threshold_pct), n = Inf, width = Inf)
print(filter(ph15, recording) |> select(site, ice_year, advance_date, retreat_date,
                                        duration_days, ice_days, max_pct, mean_pct), width = Inf)
print(phase_tbl |> filter(species != "BW58") |> select(-mean_ice_pct) |>
        pivot_wider(names_from = phase, values_from = c(detection_days, recording_days, pct_days_detected),
                    id_cols = c(site, species)) |> select(site, species, starts_with("pct")), width = Inf)
