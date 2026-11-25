# =============================================================================
# Script  : 08_site_year_tests.R
# Section : Supplement S7 - are the sites and the recording years different?
# Purpose : Statistical tests that complement 07_site_year_comparison.R.
#           Years are the replicates: one mean per site per year, which avoids
#           the day-to-day autocorrelation that would make tests on daily
#           values far too optimistic.
#           (A) Do the sites differ? Annual (calendar-year) means; two-way
#               ANOVA value ~ year + site (years = blocks), F test for site;
#               pairwise paired t-tests between sites (Holm-adjusted within
#               each variable) and the number of years in which the pairwise
#               ordering held. Repeated for the window covered by all three
#               recordings (5 Mar - 17 Jul) as a sensitivity check.
#           (B) Was the recording period unusual for its site? Mean over the
#               recording window compared with the same calendar window in
#               every other year at that site: rank, exact two-sided rank
#               p-value and z-score.
#           (C) How much of the difference between two recordings is site and
#               how much is year? Within the common window, the observed
#               difference between recordings (e.g. EI 2014 - KGI 2015) is
#               split into the long-term site difference (mean over all
#               years) and the remainder (the recording years' anomalies).
#           (D) Figure of the EI 2014 mid-water/deep temperature anomaly.
# Inputs  : F_ENV_DAILY (daily predictors per site, 2009-2021; ice 2013-2021)
# Outputs : outputs/05_environmental_data/site_year_tests_sites.csv
#           outputs/05_environmental_data/site_year_tests_pairs.csv
#           outputs/05_environmental_data/site_year_tests_recording_rank.csv
#           outputs/05_environmental_data/site_year_tests_recording_split.csv
#           outputs/05_environmental_data/site_year_tests.png/.pdf
#           outputs/05_environmental_data/ei2014_temperature.png/.pdf
# Author  : Natalie Posdaljian
# Date    : 2026-11-25
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Means of consecutive years are treated as independent; the
#           largest lag-1 correlation of the yearly means is reported as a
#           check. Raw values (not anomalies) are used. A year counts only if
#           >= 90% of its days have data, and only years with all three sites
#           are used in (A).
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages(library(tidyverse))

VARS <- c(ice_conc = "Sea-ice concentration (%)",
          temperature_surface = "Temperature, surface (°C)",
          temperature_mid = "Temperature, mid-water (°C)",
          temperature_deep = "Temperature, deep (°C)",
          mixed_layer = "Mixed layer depth (m)",
          FSLE = "FSLE (per day)",
          EKE_mid = "EKE, mid-water (cm2 s-2)",
          EKE_deep = "EKE, deep (cm2 s-2)")
COMMON <- c(start = "03-05", end = "07-17")       # covered by all three recordings
RECORDING <- tibble(site = c("EI", "KGI", "CI"),
                    start = as.Date(c("2014-03-05", "2015-02-10", "2016-02-04")),
                    end   = as.Date(c("2014-07-17", "2016-01-29", "2016-12-02")))
RECORDING_YEAR <- c(EI = 2014, KGI = 2015, CI = 2016)
PAIRS <- list(c("EI", "KGI"), c("EI", "CI"), c("KGI", "CI"))

env <- read_csv(F_ENV_DAILY, show_col_types = FALSE) |>
  mutate(date = as.Date(date), ice_conc = 100 * ice_conc,
         year = as.integer(format(date, "%Y")), md = format(date, "%m-%d"))
long <- env |> select(site, date, year, md, all_of(names(VARS))) |>
  pivot_longer(all_of(names(VARS)), names_to = "variable")

# Mean per site and year within a calendar window (>= 90% of days with data)
yearly_means <- function(start_md, end_md) {
  n_days <- as.integer(as.Date(paste0("2015-", end_md)) - as.Date(paste0("2015-", start_md))) + 1
  long |>
    filter(md >= start_md, md <= end_md) |>
    group_by(variable, site, year) |>
    summarise(n = sum(!is.na(value)), value = mean(value, na.rm = TRUE), .groups = "drop") |>
    filter(n >= 0.9 * n_days) |> select(-n) |>
    mutate(site = factor(site, SITE_ORDER)) |>
    group_by(variable, year) |> filter(n_distinct(site) == 3) |> ungroup()
}
yearly <- list(annual = yearly_means("01-01", "12-31"),
               common_window = yearly_means(COMMON[["start"]], COMMON[["end"]]))

# -----------------------------------------------------------------------------
# A. Do the sites differ?
# -----------------------------------------------------------------------------
site_test <- function(d) {
  if (var(d$value) == 0) return(tibble())
  a <- anova(lm(value ~ factor(year) + site, data = d))
  ss <- setNames(a$`Sum Sq`, rownames(a))
  lag1 <- d |> arrange(site, year) |> group_by(site) |>
    summarise(r = suppressWarnings(cor(head(value, -1), tail(value, -1))), .groups = "drop")
  tibble(n_years = n_distinct(d$year),
         years = paste(range(d$year), collapse = "-"),
         F_site = a["site", "F value"], df_resid = a["Residuals", "Df"],
         p_site = a["site", "Pr(>F)"],
         pct_var_site = 100 * ss[["site"]] / sum(ss),
         pct_var_year = 100 * ss[["factor(year)"]] / sum(ss),
         max_abs_lag1_r = max(abs(lag1$r), na.rm = TRUE))
}
pair_test <- function(d) {
  w <- d |> pivot_wider(names_from = site, values_from = value)
  map_dfr(PAIRS, function(p) {
    diff <- w[[p[1]]] - w[[p[2]]]
    if (all(diff == 0)) return(tibble())
    tt <- t.test(diff)
    tibble(pair = paste(p, collapse = " - "), mean_diff = mean(diff),
           ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], p = tt$p.value,
           n_same_sign = max(sum(diff > 0), sum(diff < 0)), n_years = length(diff))
  }) |> mutate(p_holm = p.adjust(p, "holm"))
}
site_tests <- imap_dfr(yearly, \(y, w) y |> group_by(variable) |>
                         group_modify(\(d, k) site_test(d)) |> ungroup() |> mutate(window = w, .before = 1))
pair_tests <- imap_dfr(yearly, \(y, w) y |> group_by(variable) |>
                         group_modify(\(d, k) pair_test(d)) |> ungroup() |> mutate(window = w, .before = 1))

# -----------------------------------------------------------------------------
# B. Was the recording period unusual for its site? (own recording window)
# -----------------------------------------------------------------------------
shift_mean <- function(d, start, end, k) {
  s <- seq(start, by = paste(k, "years"), length.out = 2)[2]
  e <- seq(end,   by = paste(k, "years"), length.out = 2)[2]
  x <- d$value[d$date >= s & d$date <= e]
  if (sum(!is.na(x)) < 0.9 * as.integer(end - start + 1)) NA_real_ else mean(x, na.rm = TRUE)
}
rec_rank <- pmap_dfr(RECORDING, function(site, start, end) {
  map_dfr(names(VARS), function(v) {
    d <- filter(long, .data$site == !!site, variable == v)
    shifts <- -10:10
    m <- map_dbl(shifts, \(k) shift_mean(d, start, end, k))
    keep <- !is.na(m); m <- m[keep]; is_rec <- shifts[keep] == 0
    if (!any(is_rec) || var(m) == 0) return(tibble())
    n <- length(m); r <- rank(m)[is_rec]; other <- m[!is_rec]
    tibble(site = site, variable = v, n_years = n, recording_value = m[is_rec],
           other_mean = mean(other), other_sd = sd(other),
           z = (m[is_rec] - mean(other)) / sd(other),
           rank_low_to_high = r,
           p_rank_two_sided = min(1, 2 * min(r, n + 1 - r) / n))
  })
})

# -----------------------------------------------------------------------------
# C. Difference between two recordings = site part + year part (common window)
# -----------------------------------------------------------------------------
rec_split <- yearly$common_window |> group_by(variable) |> group_modify(function(d, key) {
  site_mean <- d |> group_by(site) |> summarise(m = mean(value), .groups = "drop") |> deframe()
  rec <- map_dbl(names(RECORDING_YEAR), \(s) {
    v <- d$value[d$site == s & d$year == RECORDING_YEAR[[s]]]
    if (length(v)) v else NA_real_
  }) |> setNames(names(RECORDING_YEAR))
  map_dfr(PAIRS, function(p) {
    obs <- rec[[p[1]]] - rec[[p[2]]]
    site_part <- site_mean[[p[1]]] - site_mean[[p[2]]]
    tibble(pair = sprintf("%s %d - %s %d", p[1], RECORDING_YEAR[[p[1]]], p[2], RECORDING_YEAR[[p[2]]]),
           observed_diff = obs, site_part = site_part, year_part = obs - site_part)
  })
}) |> ungroup()

sig3 <- \(d) mutate(d, across(where(is.double), \(x) signif(x, 3)))
write_csv(sig3(site_tests), out_path("05_environmental_data", "site_year_tests_sites.csv"))
write_csv(sig3(pair_tests), out_path("05_environmental_data", "site_year_tests_pairs.csv"))
write_csv(sig3(rec_rank),   out_path("05_environmental_data", "site_year_tests_recording_rank.csv"))
write_csv(sig3(rec_split),  out_path("05_environmental_data", "site_year_tests_recording_split.csv"))

# -----------------------------------------------------------------------------
# Figure: annual mean per year at each site; recording years circled
# -----------------------------------------------------------------------------
ann <- yearly$annual
rec_pts <- ann |> filter(year == RECORDING_YEAR[as.character(site)])
p <- ggplot(ann, aes(site, value, colour = site)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey60") +
  geom_point(position = position_jitter(width = 0.08, seed = 1), size = 1, alpha = 0.6) +
  geom_point(data = rec_pts, shape = 21, size = 3, stroke = 1, fill = NA, colour = "black") +
  facet_wrap(~ factor(VARS[variable], VARS), scales = "free_y", ncol = 4) +
  scale_colour_manual(values = SITE_COLORS, guide = "none") +
  labs(x = NULL, y = "Annual mean") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank())
for (ext in c("png", "pdf")) ggsave(out_path("05_environmental_data", paste0("site_year_tests.", ext)), p, width = 8, height = 4.5, dpi = 300)

# -----------------------------------------------------------------------------
# D. EI 2014 temperature anomaly: daily mid-water and deep temperature at each
#    site by day of year, 2014 highlighted against all other years
# -----------------------------------------------------------------------------
ts <- env |>
  select(site, date, year, temperature_mid, temperature_deep) |>
  pivot_longer(c(temperature_mid, temperature_deep), names_to = "layer") |>
  mutate(layer = factor(c(temperature_mid = "Mid-water (100-500 m)",
                          temperature_deep = "Deep (500-1000 m)")[layer],
                        c("Mid-water (100-500 m)", "Deep (500-1000 m)")),
         site = factor(site, SITE_ORDER),
         doy = as.Date(paste0("2000-", format(date, "%m-%d"))))
rec_shade <- RECORDING |> filter(site == "EI") |>
  mutate(site = factor(site, SITE_ORDER),
         s = as.Date(paste0("2000-", format(start, "%m-%d"))),
         e = as.Date(paste0("2000-", format(end, "%m-%d"))))
p2 <- ggplot(ts, aes(doy, value, group = year)) +
  geom_rect(data = rec_shade, aes(xmin = s, xmax = e, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85") +
  geom_line(data = filter(ts, year != 2014), colour = "grey55", linewidth = 0.25) +
  geom_line(data = filter(ts, year == 2014), colour = "black", linewidth = 0.6) +
  facet_grid(layer ~ site) +
  scale_x_date(date_labels = "%b", date_breaks = "3 months") +
  labs(x = NULL, y = "Temperature (°C)") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank())
for (ext in c("png", "pdf")) ggsave(out_path("05_environmental_data", paste0("ei2014_temperature.", ext)), p2, width = 7.5, height = 4, dpi = 300)

print(site_tests, n = Inf, width = Inf)
print(filter(pair_tests, window == "annual"), n = Inf, width = Inf)
print(rec_rank, n = Inf, width = Inf)
