# =============================================================================
# Script  : 03_compare_hycom_copernicus.R
# Section : Supplementary Fig. S3 (2.5.1) — HYCOM vs Copernicus comparison
# Purpose : Quantify how closely the Copernicus reanalysis used in the habitat
#           models agrees with (a) satellite SST and (b) an independent ocean
#           model (HYCOM) at the species-relevant depths, over the same
#           +/- 40 km boxes and the same days.
#             (a) surface temperature: HYCOM vs Aqua MODIS, Copernicus vs MODIS
#             (b) species depths: HYCOM vs Copernicus, temperature and salinity
#           Metrics: bias (mean difference), RMSE, Pearson r and n, on raw
#           daily values and on deseasoned anomalies.
# Inputs  : config.R
#           F_COP_PHY_DEPTH (05_environmental_data/copernicus_phy_by_depth.csv)
#             Copernicus daily box means by model level
#           DOWNLOAD_DIR/hycom/<site>/hycom_<site>_<date>.nc  (01_download_hycom.R)
#           DOWNLOAD_DIR/modis/modis_sst_<site>.nc            (02_download_modis_sst.R)
# Outputs : out_path("05b_hycom_comparison", ...)
#           hycom_daily_box_means.csv   date, site, depth_m (Copernicus level),
#                                       temperature, salinity, coverage
#           modis_daily_box_means.csv   date, site, sst, coverage
#           comparison_metrics.csv      comparison, variable, site,
#                                       nominal_depth_m, depth_m, values,
#                                       n_days, bias, rmse, r
#           hycom_copernicus_profile.csv  raw-value bias/RMSE/r at every
#                                       Copernicus level (figure panel b)
#           FigS3_hycom_copernicus.png (600 dpi) / .pdf
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   :
#   * All three products are reduced to daily means over the SAME main box
#     (+/- 20 km, SITES$box_*; downloads cover the +/- 40 km box)
#     (grid cells whose centres fall inside the box) and joined BY DATE. The
#     original script (Environmental Data/Cop_vs_HYCOM.R) combined the series
#     with cbind() by row position, which misaligns dates whenever a product
#     has a missing day; that is fixed here.
#   * Vertical matching: each nominal species depth (ALL_DEPTHS_M) is mapped to
#     the nearest Copernicus level, as in the main analysis. HYCOM is linearly
#     interpolated, cell by cell, to that same Copernicus level depth before
#     box-averaging, so both models are compared at the same depth; cells
#     whose HYCOM bottom lies above the target depth are excluded, matching
#     Copernicus (cells deeper than the level only).
#     [CONFIRM interpolation of HYCOM to Copernicus levels; the alternative is
#     the nearest HYCOM level, which can be up to 25 m away (e.g. 350/400 m
#     for 375 m)].
#   * "Surface" is the top Copernicus level (0.494 m); HYCOM is interpolated
#     to it from its 0 and 2 m levels.
#   * MODIS days are used only if at least MODIS_MIN_COVERAGE of the box's
#     ocean pixels are cloud- and ice-free, because a mean of a few pixels is
#     not comparable to a full-box model mean. [CONFIRM threshold]
#   * Both models assimilate satellite SST, so (a) is not a fully independent
#     validation, and (b) measures agreement between two models, not accuracy.
#   * Deseasoned anomalies: each product's value minus that product's own
#     mean for the same site and calendar month, computed over the matched
#     comparison days (months with < MIN_DAYS_PER_MONTH matched days are
#     dropped). The bias of anomalies is zero by construction and is
#     reported as NA.
#     [CONFIRM the main analysis deseasons against a 10-yr monthly
#     climatology smoothed with a periodic spline; no 10-yr HYCOM record is
#     downloaded here, so the anomaly comparison uses within-period monthly
#     means instead.]
#   * No hypothesis tests (Shapiro-Wilk, paired t, Wilcoxon, Diebold-Mariano)
#     are applied: daily values are strongly autocorrelated and the purpose
#     is to describe agreement, not to test a null of equality.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

suppressPackageStartupMessages({
  library(ncdf4)
  library(tidyverse)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
SECTION            <- "05b_hycom_comparison"
HYCOM_DIR          <- file.path(DOWNLOAD_DIR, "hycom")
MODIS_DIR          <- file.path(DOWNLOAD_DIR, "modis")
MODIS_MIN_COVERAGE <- 0.2  # min fraction of box ocean pixels with valid SST
MIN_DAYS_PER_MONTH <- 5    # min matched days in a month to compute anomalies

PRODUCT_COLORS <- c(Copernicus = "#0072B2", HYCOM = "#D55E00", MODIS = "#333333")

if (!file.exists(F_COP_PHY_DEPTH)) {
  stop("Copernicus input not found: ", F_COP_PHY_DEPTH,
       "\nRun the 05_environmental_data extraction first.")
}

# -----------------------------------------------------------------------------
# 2. Copernicus daily box means (already computed in 05_environmental_data)
# -----------------------------------------------------------------------------
cop <- read_csv(F_COP_PHY_DEPTH, show_col_types = FALSE) %>%
  transmute(date = as.Date(date), site, depth_m, temperature, salinity, coverage)

cop_levels <- sort(unique(cop$depth_m))
surface_m  <- min(cop_levels)

# Nominal species depths -> nearest Copernicus level (as in the main analysis).
depth_map <- tibble(nominal_depth_m = ALL_DEPTHS_M) %>%
  mutate(depth_m = map_dbl(nominal_depth_m,
                           ~ cop_levels[which.min(abs(cop_levels - .x))]))

# -----------------------------------------------------------------------------
# 3. HYCOM daily box means at the Copernicus levels
# -----------------------------------------------------------------------------
# Linear interpolation of a [lon, lat, depth] array to depth z, cell by cell.
# A cell is NA if either bracketing level is below the HYCOM sea floor.
interp_to_depth <- function(a, depths, z) {
  k <- findInterval(z, depths)
  if (depths[k] == z) return(a[, , k])
  w <- (z - depths[k]) / (depths[k + 1] - depths[k])
  (1 - w) * a[, , k] + w * a[, , k + 1]
}

# Downloads cover the +/- 40 km box (SITES$box40_*); every product is cropped
# here to the MAIN +/- 20 km box (SITES$box_*, cell centres inside the box),
# the box used for the Copernicus predictors. Longitudes in 0..360 are
# converted to -180..180.
in_main_box <- function(site, lon, lat) {
  s   <- SITES[SITES$site == site, ]
  lon <- ifelse(lon > 180, lon - 360, lon)
  list(i_lon = which(lon >= s$box_lon_min & lon <= s$box_lon_max),
       i_lat = which(lat >= s$box_lat_min & lat <= s$box_lat_max))
}

read_hycom_day <- function(file, site, targets) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  depths <- as.vector(ncvar_get(nc, "depth"))
  date   <- as.Date(as.POSIXct("2000-01-01", tz = "UTC") + as.vector(ncvar_get(nc, "time")) * 3600)
  box    <- in_main_box(site, as.vector(ncvar_get(nc, "lon")),
                        as.vector(ncvar_get(nc, "lat")))
  temp   <- ncvar_get(nc, "water_temp", collapse_degen = FALSE)[box$i_lon, box$i_lat, , 1, drop = FALSE]
  salt   <- ncvar_get(nc, "salinity",   collapse_degen = FALSE)[box$i_lon, box$i_lat, , 1, drop = FALSE]
  temp   <- array(temp, dim(temp)[1:3])
  salt   <- array(salt, dim(salt)[1:3])

  map_dfr(targets, function(z) {
    t_z <- interp_to_depth(temp, depths, z)
    s_z <- interp_to_depth(salt, depths, z)
    tibble(date = date, site = site, depth_m = z,
           temperature = mean(t_z, na.rm = TRUE),
           salinity    = mean(s_z, na.rm = TRUE),
           coverage    = mean(!is.na(t_z)))
  })
}

hycom_files <- tibble(site = SITE_ORDER) %>%
  mutate(file = map(site, ~ list.files(file.path(HYCOM_DIR, .x),
                                       pattern = "\\.nc$", full.names = TRUE))) %>%
  unnest(file)
if (nrow(hycom_files) == 0) stop("No HYCOM files found in ", HYCOM_DIR)

# Only levels within the downloaded HYCOM depth range can be interpolated.
hycom_max_depth <- {
  nc <- nc_open(hycom_files$file[1])
  d  <- max(ncvar_get(nc, "depth"))
  nc_close(nc)
  d
}
targets <- cop_levels[cop_levels <= hycom_max_depth]

hycom <- map2_dfr(hycom_files$file, hycom_files$site, read_hycom_day,
                  targets = targets) %>%
  mutate(across(c(temperature, salinity), ~ if_else(is.nan(.x), NA_real_, .x)))
write_csv(hycom, out_path(SECTION, "hycom_daily_box_means.csv"))

# -----------------------------------------------------------------------------
# 4. MODIS daily box means
# -----------------------------------------------------------------------------
# Coverage is relative to the box's ocean pixels, defined as pixels with a
# valid SST on at least one day of the record (land is never valid).
read_modis_site <- function(site) {
  nc <- nc_open(file.path(MODIS_DIR, sprintf("modis_sst_%s.nc", site)))
  on.exit(nc_close(nc))
  box  <- in_main_box(site, as.vector(ncvar_get(nc, "longitude")),
                      as.vector(ncvar_get(nc, "latitude")))
  grid <- ncvar_get(nc, "sstMasked", collapse_degen = FALSE)[box$i_lon, box$i_lat, , drop = FALSE]  # lon, lat, time
  date <- as.Date(as.POSIXct(as.vector(ncvar_get(nc, "time")), origin = "1970-01-01", tz = "UTC"))
  ocean <- apply(!is.na(grid), c(1, 2), any)
  tibble(date = date, site = site,
         sst      = apply(grid, 3, mean, na.rm = TRUE),
         coverage = apply(grid, 3, function(x) sum(!is.na(x[ocean]))) / sum(ocean)) %>%
    mutate(sst = if_else(is.nan(sst), NA_real_, sst))
}
modis <- map_dfr(SITE_ORDER, read_modis_site)
write_csv(modis, out_path(SECTION, "modis_daily_box_means.csv"))

# -----------------------------------------------------------------------------
# 5. Metrics
# -----------------------------------------------------------------------------
# Deseasoned anomalies: value minus the product's own site x calendar-month
# mean over the matched days (see Notes).
add_anomalies <- function(pairs) {
  pairs %>%
    group_by(site, month = format(date, "%m")) %>%
    filter(n() >= MIN_DAYS_PER_MONTH) %>%
    mutate(x = x - mean(x), y = y - mean(y)) %>%
    ungroup() %>%
    select(-month)
}

# x = product being evaluated, y = reference; bias = mean(x - y).
agreement <- function(x, y) {
  tibble(n_days = length(x),
         bias   = mean(x - y),
         rmse   = sqrt(mean((x - y)^2)),
         r      = if (length(x) > 2) cor(x, y) else NA_real_)
}

# pairs: site, date, x, y plus grouping columns -> metrics on raw values and
# on anomalies.
compute_metrics <- function(pairs, ...) {
  raw  <- pairs %>% group_by(site, ...) %>%
    reframe(agreement(x, y)) %>% mutate(values = "raw")
  anom <- pairs %>% group_by(...) %>% group_modify(~ add_anomalies(.x)) %>%
    group_by(site, ...) %>%
    reframe(agreement(x, y)) %>% mutate(values = "anomaly", bias = NA_real_)
  bind_rows(raw, anom)
}

# (a) Surface temperature vs MODIS, on days where all three products exist,
#     so that HYCOM and Copernicus are scored on the same days.
surface <- cop %>%
  filter(depth_m == surface_m) %>%
  select(site, date, Copernicus = temperature) %>%
  inner_join(hycom %>% filter(depth_m == surface_m) %>%
               select(site, date, HYCOM = temperature), by = c("site", "date")) %>%
  inner_join(modis %>% filter(coverage >= MODIS_MIN_COVERAGE) %>%
               select(site, date, MODIS = sst), by = c("site", "date")) %>%
  drop_na()

metrics_surface <- surface %>%
  pivot_longer(c(Copernicus, HYCOM), names_to = "model", values_to = "x") %>%
  rename(y = MODIS) %>%
  compute_metrics(model) %>%
  transmute(comparison = paste(model, "vs MODIS"), variable = "temperature",
            site, nominal_depth_m = 0, depth_m = surface_m,
            values, n_days, bias, rmse, r)

# (b) HYCOM vs Copernicus at every Copernicus level within the HYCOM range.
model_pairs <- inner_join(
  hycom %>% pivot_longer(c(temperature, salinity), names_to = "variable", values_to = "x"),
  cop   %>% pivot_longer(c(temperature, salinity), names_to = "variable", values_to = "y"),
  by = c("site", "date", "depth_m", "variable"), suffix = c("_hycom", "_cop")
) %>%
  drop_na(x, y)

metrics_depth <- model_pairs %>%
  semi_join(depth_map, by = "depth_m") %>%
  compute_metrics(variable, depth_m) %>%
  inner_join(depth_map, by = "depth_m", relationship = "many-to-many") %>%
  transmute(comparison = "HYCOM vs Copernicus", variable, site,
            nominal_depth_m, depth_m, values, n_days, bias, rmse, r)

metrics <- bind_rows(metrics_surface, metrics_depth) %>%
  arrange(comparison, variable, factor(site, SITE_ORDER), nominal_depth_m,
          desc(values))
write_csv(metrics, out_path(SECTION, "comparison_metrics.csv"))

# Raw-value profile at all levels (figure panel b), with the mean fraction of
# box cells holding water in each model at that level.
profile <- model_pairs %>%
  group_by(site, variable, depth_m) %>%
  reframe(agreement(x, y),
          coverage_hycom = mean(coverage_hycom),
          coverage_cop   = mean(coverage_cop))
write_csv(profile, out_path(SECTION, "hycom_copernicus_profile.csv"))

# -----------------------------------------------------------------------------
# 6. Figure
# -----------------------------------------------------------------------------
site_labels <- set_names(paste0(SITES$site_name, " (", SITES$site, ")"), SITES$site)
theme_s1 <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", hjust = 0),
        legend.position  = "top")

# (a) Surface temperature time series: model lines, MODIS days as points.
ts_models <- bind_rows(
  cop   %>% filter(depth_m == surface_m) %>% mutate(product = "Copernicus"),
  hycom %>% filter(depth_m == surface_m) %>% mutate(product = "HYCOM")
) %>%
  semi_join(hycom %>% distinct(site, date), by = c("site", "date"))
ts_modis <- modis %>%
  filter(coverage >= MODIS_MIN_COVERAGE) %>%
  semi_join(ts_models, by = c("site", "date"))

p_ts <- ggplot(ts_models, aes(date, temperature, colour = product)) +
  geom_line(linewidth = 0.5) +
  geom_point(data = ts_modis, aes(y = sst, colour = "MODIS"), size = 0.9) +
  facet_wrap(~ factor(site, SITE_ORDER, site_labels[SITE_ORDER]),
             ncol = 1, scales = "free_x") +
  scale_colour_manual(values = PRODUCT_COLORS, breaks = names(PRODUCT_COLORS),
                      name = NULL) +
  scale_x_date(date_labels = "%b %Y") +
  labs(x = NULL, y = "Surface temperature (\u00B0C)", tag = "a") +
  theme_s1

# (b) HYCOM - Copernicus bias and RMSE profiles; species depths dashed.
profile_long <- profile %>%
  pivot_longer(c(bias, rmse), names_to = "metric") %>%
  mutate(panel = factor(
    paste(str_to_sentence(variable), if_else(metric == "bias", "bias", "RMSE")),
    levels = c("Temperature bias", "Temperature RMSE", "Salinity bias", "Salinity RMSE")))

p_prof <- ggplot(profile_long, aes(value, depth_m, colour = site)) +
  geom_hline(yintercept = unique(depth_map$depth_m), colour = "grey75",
             linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_path(linewidth = 0.6) +
  geom_point(data = semi_join(profile_long, depth_map, by = "depth_m"), size = 1.2) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  scale_y_continuous(trans = scales::compose_trans("sqrt", "reverse"),
                     breaks = c(0, 50, 100, 250, 500, 1000, 1500)) +
  scale_colour_manual(values = SITE_COLORS, breaks = SITE_ORDER, name = NULL) +
  labs(x = "HYCOM \u2212 Copernicus: bias and RMSE (\u00B0C or psu)",
       y = "Depth (m)", tag = "b") +
  theme_s1

fig <- p_ts / p_prof + plot_layout(heights = c(1.3, 1))
ggsave(out_path(SECTION, "FigS3_hycom_copernicus.png"), fig,
       width = 180, height = 210, units = "mm", dpi = 600)
ggsave(out_path(SECTION, "FigS3_hycom_copernicus.pdf"), fig,
       width = 180, height = 210, units = "mm", device = cairo_pdf)

print(metrics, n = Inf)
