# =============================================================================
# Script  : layer_boundaries.R
# Section : 2.2 Deep oceanography / 2.5 Environmental data (depth layers)
# Purpose : Find depth boundaries that split the water column at each site
#           into surface, mid-water and deep layers, from two independent
#           lines of evidence, so that layer-averaged predictors can replace
#           species-specific depths:
#             (1) coherence: depths are grouped by how similarly their daily
#                 temperature and salinity anomalies vary over the 10-year
#                 record, with a depth-constrained clustering (each group is
#                 a continuous depth range; CONISS-type, Grimm 1987);
#             (2) water-mass structure from the 10-year mean profiles: deepest
#                 winter mixed layer, Winter Water temperature minimum, main
#                 pycnocline (maximum density gradient below the mixed layer),
#                 and the warm-water (Circumpolar Deep Water) temperature
#                 maximum.
# Inputs  : F_COP_PHY_DEPTH (copernicus_phy_by_depth.csv: daily box means of
#             temperature and salinity at every model level, 10-yr window)
#           F_COP_BGC_DEPTH (copernicus_bgc_by_depth.csv: oxygen, for the
#             mean profiles)
#           F_COP_SURFACE   (copernicus_surface.csv: mixed layer depth)
#           config.R (SITES, SPECIES depths)
# Outputs : outputs/02_deep_oceanography/layer_boundaries.csv
#             site, method, boundary, depth_m
#           outputs/02_deep_oceanography/layer_boundaries.png / .pdf
#             per site: mean T, S, O2 and sigma-theta profiles with the
#             water-mass markers; inter-depth anomaly correlation matrix with
#             the coherence boundaries; species depths marked
# Author  : Natalie Posdaljian
# Date    : 2026-11-19
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Depth range MAX_DEPTH_M = 1000 m (layer averages stop there;
#             below ~1060 m less than half of the KGI and CI boxes has water).
#           - Profiles are interpolated to a regular DZ_M grid before the
#             clustering, so the dense model levels near the surface do not
#             outweigh the sparse deep levels.
#           - Anomalies: daily value minus the site x level monthly mean
#             (seasonal cycle removed), standardized per depth; temperature
#             and salinity series are concatenated, so a depth's "signature"
#             is its joint T-S variability.
#           - The number of coherence layers is fixed at N_LAYERS = 3
#             (surface, mid-water, deep); the 2- and 4-layer cuts are also
#             written for comparison.
#           - sigma-theta at zero pressure (UNESCO EOS-80, oce), as in
#             deep_oceanography.R.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(oce)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
MAX_DEPTH_M <- 1000
DZ_M        <- 5          # regular depth grid for the clustering
N_LAYERS    <- 3
WINTER_MLD_QUANTILE <- 0.9  # "deepest winter mixed layer" = 90th percentile
                            # of daily MLD in the 10-yr record
out_section <- "02_deep_oceanography"

# -----------------------------------------------------------------------------
# 2. Data
# -----------------------------------------------------------------------------
phy <- read_csv(F_COP_PHY_DEPTH, show_col_types = FALSE) |>
  select(date, site, depth_m, temperature, salinity) |>
  filter(depth_m <= MAX_DEPTH_M + 100, !is.na(temperature), !is.na(salinity)) |>
  mutate(date = as.Date(date), month = as.integer(format(date, "%m")))
o2 <- read_csv(F_COP_BGC_DEPTH, show_col_types = FALSE) |>
  select(date, site, depth_m, o2) |>
  filter(depth_m <= MAX_DEPTH_M + 100, !is.na(o2))
mld <- read_csv(F_COP_SURFACE, show_col_types = FALSE) |>
  select(date, site, mixed_layer)

grid_depths <- seq(DZ_M, MAX_DEPTH_M, by = DZ_M)

# -----------------------------------------------------------------------------
# 3. Line of evidence 1: coherence clustering
# -----------------------------------------------------------------------------
# Standardized anomaly series on the regular grid: one column per depth,
# rows = days (temperature block on top of salinity block).
anomaly_matrix <- function(d) {
  long <- d |>
    group_by(depth_m, month) |>
    mutate(across(c(temperature, salinity), ~ .x - mean(.x))) |>
    ungroup() |>
    arrange(date, depth_m)
  interp_var <- function(v) {
    long |>
      group_by(date) |>
      filter(n() >= 2) |>
      reframe(depth = grid_depths,
              value = approx(depth_m, .data[[v]], xout = grid_depths, rule = 1)$y) |>
      pivot_wider(names_from = depth, values_from = value) |>
      select(-date) |>
      as.matrix()
  }
  m <- rbind(interp_var("temperature"), interp_var("salinity"))
  keep <- colSums(is.na(m)) == 0           # depths with water every day
  scale(m[, keep])
}

# Depth-constrained (contiguous) agglomerative clustering, Ward criterion:
# at each step merge the adjacent pair of groups whose merge increases the
# within-group sum of squares least (CONISS; Grimm 1987). Returns the
# boundaries (index of the last depth in each upper group) for k groups.
coniss_boundaries <- function(x, k_values) {
  n <- ncol(x)
  groups <- lapply(seq_len(n), function(i) list(idx = i, n = 1, mean = x[, i]))
  out <- list()
  while (length(groups) > 1) {
    if (length(groups) %in% k_values) {
      out[[as.character(length(groups))]] <-
        vapply(groups[-length(groups)], function(g) max(g$idx), numeric(1))
    }
    cost <- vapply(seq_len(length(groups) - 1), function(i) {
      a <- groups[[i]]; b <- groups[[i + 1]]
      (a$n * b$n / (a$n + b$n)) * sum((a$mean - b$mean)^2)
    }, numeric(1))
    i <- which.min(cost)
    a <- groups[[i]]; b <- groups[[i + 1]]
    groups[[i]] <- list(idx = c(a$idx, b$idx), n = a$n + b$n,
                        mean = (a$n * a$mean + b$n * b$mean) / (a$n + b$n))
    groups[[i + 1]] <- NULL
  }
  out
}

coherence <- map(set_names(SITE_ORDER), function(s) {
  m <- anomaly_matrix(filter(phy, site == s))
  depths <- as.numeric(colnames(m))
  cuts <- coniss_boundaries(m, k_values = 2:4)
  list(
    cor = cor(m),
    depths = depths,
    boundaries = imap(cuts, \(b, k) tibble(
      site = s, method = paste0("coherence (", k, " layers)"),
      boundary = paste0("b", seq_along(b)),
      # boundary placed halfway between the last depth above and the first below
      depth_m = (depths[b] + depths[b + 1]) / 2)) |>
      list_rbind()
  )
})

# -----------------------------------------------------------------------------
# 4. Line of evidence 2: water-mass structure from mean profiles
# -----------------------------------------------------------------------------
mean_profiles <- phy |>
  group_by(site, depth_m) |>
  summarise(temperature = mean(temperature), salinity = mean(salinity),
            .groups = "drop") |>
  mutate(sigma_theta = swSigmaTheta(salinity, temperature, pressure = 0,
                                    eos = "unesco")) |>
  left_join(o2 |> group_by(site, depth_m) |>
              summarise(o2 = mean(o2), .groups = "drop") |>
              rename(o2_depth = depth_m) |>
              group_by(site) |>
              nest(o2_prof = c(o2_depth, o2)),
            by = "site") |>
  mutate(o2 = map2_dbl(depth_m, o2_prof, \(z, p)
                       approx(p$o2_depth, p$o2, xout = z, rule = 1)$y)) |>
  select(-o2_prof)

# Summer (Dec-Mar) mean temperature profile, for the Winter Water minimum,
# which is a summer feature (the remnant of the previous winter mixed layer).
summer_temp <- phy |>
  filter(month %in% c(12, 1, 2, 3)) |>
  group_by(site, depth_m) |>
  summarise(temperature = mean(temperature), .groups = "drop")

water_mass <- map(SITE_ORDER, function(s) {
  prof <- filter(mean_profiles, site == s) |> arrange(depth_m)
  summ <- filter(summer_temp, site == s) |> arrange(depth_m)
  winter_mld <- quantile(filter(mld, site == s)$mixed_layer, WINTER_MLD_QUANTILE,
                         na.rm = TRUE)
  ww <- summ |> filter(depth_m > 10, depth_m < 250) |> slice_min(temperature, n = 1)
  # Pycnocline: maximum d(sigma-theta)/dz below the winter mixed layer
  dz <- diff(prof$depth_m)
  grad <- tibble(depth_m = (head(prof$depth_m, -1) + prof$depth_m[-1]) / 2,
                 dsig = diff(prof$sigma_theta) / dz) |>
    filter(depth_m > winter_mld, depth_m < MAX_DEPTH_M)
  pycno <- grad |> slice_max(dsig, n = 1)
  tmax <- prof |> filter(depth_m > 100, depth_m <= MAX_DEPTH_M) |>
    slice_max(temperature, n = 1)
  tibble(site = s, method = "water mass",
         boundary = c("winter mixed layer (90th pct)", "Winter Water T-min (summer)",
                      "main pycnocline", "warm-water T-max"),
         depth_m = c(as.numeric(winter_mld), ww$depth_m, pycno$depth_m, tmax$depth_m))
}) |>
  list_rbind()

boundaries <- bind_rows(map(coherence, "boundaries") |> list_rbind(), water_mass) |>
  mutate(depth_m = round(depth_m))
write_csv(boundaries, out_path(out_section, "layer_boundaries.csv"))

# -----------------------------------------------------------------------------
# 5. Figure
# -----------------------------------------------------------------------------
species_depths <- enframe(SPECIES, name = "species", value = "info") |>
  mutate(depth_m = map(info, "depths_m")) |>
  select(species, depth_m) |>
  unnest(depth_m) |>
  filter(species %in% MODELED_SPECIES, depth_m > 0, depth_m <= MAX_DEPTH_M + 700)

sqrt_rev <- scales::transform_compose(scales::transform_sqrt(),
                                      scales::transform_reverse())
depth_breaks <- c(0, 50, 100, 200, 400, 600, 800, 1000)
coh3 <- boundaries |> filter(method == "coherence (3 layers)")
wm   <- boundaries |> filter(method == "water mass")
wm_colors <- c("winter mixed layer (90th pct)" = "#1b9e77",
               "Winter Water T-min (summer)" = "#7570b3",
               "main pycnocline" = "#d95f02",
               "warm-water T-max" = "#e7298a")

profile_panel <- function(s) {
  prof <- mean_profiles |>
    filter(site == s) |>
    pivot_longer(c(temperature, salinity, o2, sigma_theta), names_to = "var") |>
    mutate(var = factor(var, levels = c("temperature", "salinity", "sigma_theta", "o2"),
                        labels = c("Temperature (°C)", "Salinity", "σθ (kg m⁻³)",
                                   "Oxygen (mmol m⁻³)")))
  ggplot(prof, aes(value, depth_m)) +
    geom_hline(data = filter(coh3, site == s), aes(yintercept = depth_m),
               linewidth = 0.8, colour = "black") +
    geom_hline(data = filter(wm, site == s), aes(yintercept = depth_m, colour = boundary),
               linetype = "dashed", linewidth = 0.6) +
    geom_path(na.rm = TRUE) +
    facet_wrap(~ var, nrow = 1, scales = "free_x") +
    scale_y_continuous(transform = sqrt_rev, breaks = depth_breaks,
                       limits = c(MAX_DEPTH_M, 0), expand = c(0, 0)) +
    scale_colour_manual(values = wm_colors, name = "Water-mass marker") +
    labs(x = NULL, y = "Depth (m)", title = s) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

cor_panel <- function(s) {
  co <- coherence[[s]]
  df <- as_tibble(co$cor) |>
    mutate(d1 = co$depths) |>
    pivot_longer(-d1, names_to = "d2", values_to = "r") |>
    mutate(d2 = as.numeric(d2))
  ggplot(df, aes(d1, d2, fill = r)) +
    geom_tile(width = DZ_M, height = DZ_M) +
    geom_hline(data = filter(coh3, site == s), aes(yintercept = depth_m), linewidth = 0.6) +
    geom_vline(data = filter(coh3, site == s), aes(xintercept = depth_m), linewidth = 0.6) +
    geom_rug(data = species_depths, aes(x = depth_m), inherit.aes = FALSE,
             colour = "grey20", length = unit(0.06, "npc")) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-1, 1), name = "r") +
    scale_x_continuous(transform = "sqrt", breaks = depth_breaks, expand = c(0, 0)) +
    scale_y_continuous(transform = sqrt_rev, breaks = depth_breaks, expand = c(0, 0)) +
    coord_equal() +
    labs(x = "Depth (m)", y = "Depth (m)",
         title = paste(s, "- correlation of T-S anomalies between depths")) +
    theme_bw(base_size = 9)
}

fig <- wrap_plots(
  map(SITE_ORDER, \(s) (profile_panel(s) | cor_panel(s)) +
        plot_layout(widths = c(2.2, 1))),
  ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(out_path(out_section, "layer_boundaries.png"), fig, width = 13, height = 12,
       dpi = 300, bg = "white")
ggsave(out_path(out_section, "layer_boundaries.pdf"), fig, width = 13, height = 12,
       device = cairo_pdf)

print(boundaries |> pivot_wider(names_from = site, values_from = depth_m), n = Inf)
