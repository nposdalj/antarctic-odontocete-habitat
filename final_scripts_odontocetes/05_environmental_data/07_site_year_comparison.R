# =============================================================================
# Script  : 07_site_year_comparison.R
# Section : Supplement - separating site from year in the environmental data
# Purpose : Each site was recorded in a different year (EI 2014, KGI 2015,
#           CI 2016), so site and year are confounded in the acoustic data.
#           The environmental records, however, cover all three sites in all
#           three years. This script asks, for each model predictor and sea
#           ice:
#           (1) How much of the variation among monthly means is due to site,
#               year (2014-2016) and their interaction, after removing the
#               seasonal cycle (month)? Large site and small year/interaction
#               shares mean that conditions differ mainly between places, so
#               differences between recordings are more likely site effects.
#           (2) Was each site's recording period typical for that site?
#               The mean over the recording window is compared with the same
#               calendar window in every other year of the 10-year record.
# Inputs  : F_ENV_DAILY (10-year daily predictors per site)
# Outputs : outputs/05_environmental_data/site_year_variance.csv
#           outputs/05_environmental_data/site_year_recording_percentile.csv
#           outputs/05_environmental_data/site_year_comparison.png/.pdf
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Raw values (not anomalies) are used, with month as a factor, so
#           the seasonal cycle is removed in the same way at every site.
#           FSLV orientation is axial and is summarized by its doubled-angle
#           components only in the variance table (not plotted).
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages(library(tidyverse))

YEARS <- 2014:2016
VARS <- c(ice_conc = "Sea-ice concentration (%)",
          temperature_surface = "Temperature, surface layer (°C)",
          temperature_mid = "Temperature, mid-water layer (°C)",
          temperature_deep = "Temperature, deep layer (°C)",
          mixed_layer = "Mixed layer depth (m)",
          FSLE = "FSLE (per day)",
          EKE_mid = "EKE, mid-water layer (cm2 s-2)",
          EKE_deep = "EKE, deep layer (cm2 s-2)",
          chla_lag3 = "log10 chl-a, 3-month lag (mg m-3)")
RECORDING <- tibble(site = c("EI", "KGI", "CI"),
                    start = as.Date(c("2014-03-05", "2015-02-10", "2016-02-04")),
                    end   = as.Date(c("2014-07-17", "2016-01-29", "2016-12-02")))

env <- read_csv(F_ENV_DAILY, show_col_types = FALSE) |>
  mutate(date = as.Date(date), ice_conc = 100 * ice_conc, chla_lag3 = log10(chla_lag3),
         fsle_orient_cos = cos(2 * fsle_orient * pi / 180),
         fsle_orient_sin = sin(2 * fsle_orient * pi / 180))
vars_all <- c(names(VARS), "fsle_orient_cos", "fsle_orient_sin")

# -----------------------------------------------------------------------------
# 1. Variance partition of monthly means, 2014-2016, all three sites
# -----------------------------------------------------------------------------
monthly <- env |>
  filter(as.integer(format(date, "%Y")) %in% YEARS) |>
  mutate(year = factor(format(date, "%Y")), month = factor(format(date, "%m"))) |>
  group_by(site, year, month) |>
  summarise(across(all_of(vars_all), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  mutate(site = factor(site, SITE_ORDER))

partition <- map_dfr(vars_all, function(v) {
  d <- monthly |> filter(!is.na(.data[[v]]))
  if (var(d[[v]]) == 0) return(NULL)       # e.g. no ice in any month
  a <- anova(lm(reformulate(c("month", "site", "year", "site:year"), v), data = d))
  ss <- setNames(a$`Sum Sq`, rownames(a))
  non_seasonal <- sum(ss) - ss[["month"]]
  tibble(variable = v,
         pct_site = 100 * ss[["site"]] / non_seasonal,
         pct_year = 100 * ss[["year"]] / non_seasonal,
         pct_site_x_year = 100 * ss[["site:year"]] / non_seasonal,
         pct_residual = 100 * ss[["Residuals"]] / non_seasonal,
         pct_total_seasonal = 100 * ss[["month"]] / sum(ss))
}) |> mutate(across(where(is.numeric), \(x) round(x, 1)))
write_csv(partition, out_path("05_environmental_data", "site_year_variance.csv"))

# -----------------------------------------------------------------------------
# 2. Was the recording period typical for its site? (same calendar window,
#    every year of the 10-year record)
# -----------------------------------------------------------------------------
window_mean <- function(d, start, end, shift_years, v) {
  s <- seq(start, by = paste(shift_years, "years"), length.out = 2)[2]
  e <- seq(end,   by = paste(shift_years, "years"), length.out = 2)[2]
  x <- d[[v]][d$date >= s & d$date <= e]
  if (length(x) < 0.9 * as.integer(end - start + 1)) NA_real_ else mean(x, na.rm = TRUE)
}
percentiles <- pmap_dfr(RECORDING, function(site, start, end) {
  d <- filter(env, .data$site == !!site)
  map_dfr(names(VARS), function(v) {
    shifts <- -10:10
    means <- map_dbl(shifts, \(k) window_mean(d, start, end, k, v))
    obs <- means[shifts == 0]; other <- means[shifts != 0 & !is.na(means)]
    tibble(site = site, variable = v, recording_mean = signif(obs, 3),
           other_years_median = signif(median(other), 3),
           other_years_min = signif(min(other), 3), other_years_max = signif(max(other), 3),
           n_other_years = length(other),
           percentile = round(100 * mean(other <= obs)))
  })
})
write_csv(percentiles, out_path("05_environmental_data", "site_year_recording_percentile.csv"))

# -----------------------------------------------------------------------------
# 3. Figure: monthly means at each site (colour) in each year (line type),
#    with each site's recording window highlighted
# -----------------------------------------------------------------------------
plot_dat <- monthly |>
  select(site, year, month, all_of(names(VARS))) |>
  pivot_longer(all_of(names(VARS)), names_to = "variable") |>
  mutate(variable = factor(VARS[variable], VARS), month = as.integer(month))
rec_dat <- plot_dat |>
  mutate(date = as.Date(paste(year, month, 15, sep = "-"))) |>
  inner_join(RECORDING, by = "site", relationship = "many-to-many") |>
  filter(date >= start - 15, date <= end + 15)

p <- ggplot(plot_dat, aes(month, value, colour = site, linetype = year)) +
  geom_line(linewidth = 0.4, alpha = 0.6) +
  geom_point(data = rec_dat, size = 1.3) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = SITE_COLORS, name = "Site") +
  scale_linetype_manual(values = c("2014" = "solid", "2015" = "dashed", "2016" = "dotted"),
                        name = "Year") +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        panel.grid.minor = element_blank())
for (ext in c("png", "pdf")) ggsave(out_path("05_environmental_data", paste0("site_year_comparison.", ext)), p, width = 7.5, height = 7, dpi = 300)

print(partition, width = Inf)
print(percentiles, n = Inf, width = Inf)
