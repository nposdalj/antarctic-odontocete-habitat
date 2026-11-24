# =============================================================================
# Script  : response_curves.R
# Section : 3.2 Results (GAMs) - figures
# Purpose : Partial-effect (response) curves of the final models: CENTERED
#           smooths on the log (link) scale with 95% pointwise confidence
#           intervals and a rug of the binned observations.
#             (1) Sea-ice models: one panel per species, sites as colours,
#                 significant smooths only.
#             (2) Habitat models: one figure per species x site (all terms).
#             (3) Habitat models, publication layout: a grid with one row per
#                 species x site model and one column per predictor; panels
#                 only where the predictor is in the final model.
#           The plotted values are saved as tables.
# Inputs  : <OUTPUT_ROOT>/06_statistical_analyses/sea_ice/models/*.rds
#           <OUTPUT_ROOT>/06_statistical_analyses/<RESULTS_STAGE>/models/*.rds
# Outputs : <OUTPUT_ROOT>/08_results_gams/
#             supplementary_sea_ice_response_curves.png / .pdf,
#               sea_ice_response_curves_data.csv
#             <label>_habitat_response_curves/<species>_<site>.png
#             <label>_habitat_response_curves_data.csv
#             <label>_habitat_response_curves_grid.png / .pdf
#           <label> = "main_stage3" (default) or "supplementary_<stage>"
# Author  : Natalie Posdaljian
# Date    : 2026-11-16 (updated 2026-11-19: significant sea-ice curves only,
#           rugs on all panels, curves limited to the observed range,
#           publication grid)
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Curves are centered partial effects (gratia::smooth_estimates),
#             NOT back-transformed (the original pipeline applied plogis(), a
#             logit inverse that is wrong for the Tweedie log link). 0 =
#             average log presence minutes; +0.69 = twice the average.
#           - Every curve is drawn only over the range of the predictor in
#             that model's data (no extrapolation). For the cyclic FSLV
#             orientation smooth the grid is restricted to the observed range.
#           - Confidence bands that extend far beyond the curves (sparse data
#             at the ends of the range) are clipped at the panel limits, which
#             are set from the curves themselves (CI_CLIP_FACTOR); the clipping
#             is stated in the figure legend.
#           - Sea-ice smooths with p > 0.05 are not drawn (decision N.
#             Posdaljian: no plots of non-significant terms); they remain in
#             the tables.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))
suppressPackageStartupMessages(library(patchwork))

# --- Parameters ---------------------------------------------------------------
RESULTS_STAGE <- Sys.getenv("ODONTOCETE_RESULTS_STAGE", unset = "stage3")
HABITAT_LABEL <- if (RESULTS_STAGE == "stage3") "main_stage3" else
  paste0("supplementary_", RESULTS_STAGE)
CURVE_GRID_N   <- 200  # evaluation points per smooth
CI_CLIP_FACTOR <- 1.5  # y limits = curve range expanded by this factor

# Predictor order (columns) and short labels for the publication grid
PREDICTOR_ORDER <- c("temperature_surface_anom", "temperature_mid_anom",
                     "temperature_deep_anom", "mixed_layer_anom", "FSLE",
                     "fsle_orient", "EKE_mid", "EKE_deep", "chla_lag3_anom")
PREDICTOR_SHORT <- c(
  temperature_surface_anom = "Temperature anomaly,\nsurface layer (°C)",
  temperature_mid_anom     = "Temperature anomaly,\nmid-water layer (°C)",
  temperature_deep_anom    = "Temperature anomaly,\ndeep layer (°C)",
  mixed_layer_anom         = "Mixed layer depth\nanomaly (m)",
  FSLE                     = "FSLE (day⁻¹)",
  fsle_orient              = "FSLV orientation (°)",
  EKE_mid                  = "EKE, mid-water\nlayer (cm² s⁻²)",
  EKE_deep                 = "EKE, deep layer\n(cm² s⁻²)",
  chla_lag3_anom           = "log10 chl-a std. anomaly\n(SD), 3-month lag")

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
hab_fig_dir <- file.path(RESULTS_DIR, paste0(HABITAT_LABEL, "_habitat_response_curves"))
if (!dir.exists(hab_fig_dir)) dir.create(hab_fig_dir, recursive = TRUE)

# Evaluate every reported smooth of a list of models (one row per grid point),
# restricted to the observed range of the predictor in that model's data.
model_curves <- function(models) {
  imap(models, function(m, id) {
    map(m$terms$variable, function(v) {
      obs <- range(m$data[[v]], na.rm = TRUE)
      smooth_grid(m$model, v, n = CURVE_GRID_N) |>
        filter(x >= obs[1], x <= obs[2]) |>
        mutate(variable = v)
    }) |>
      bind_rows() |>
      mutate(species = m$species, site = m$site, .before = 1)
  }) |> bind_rows()
}

# Rug data: the binned values of each term's predictor in each model
model_rugs <- function(models) {
  imap(models, function(m, id) {
    map(m$terms$variable, \(v) tibble(variable = v, x = m$data[[v]],
                                      detected = m$data[[m$species]] > 0)) |>
      bind_rows() |>
      mutate(species = m$species, site = m$site, .before = 1)
  }) |> bind_rows()
}

# Kernel density of the observed predictor values of the time bins, drawn
# along the bottom HIST_HEIGHT of a panel (replaces the rug; decision
# N. Posdaljian 2026-11-20): all bins (grey) and bins with detections (site
# colour). The detection density is multiplied by the proportion of bins with
# detections, so it lies under the all-bins curve and the two share one scale.
HIST_HEIGHT <- 0.22
density_curves <- function(x, detected, xr, yr) {
  grid <- seq(xr[1], xr[2], length.out = 256)
  dens <- function(v) if (length(v) < 2) rep(0, length(grid)) else
    approx(density(v, from = xr[1], to = xr[2], n = 256)$x,
           density(v, from = xr[1], to = xr[2], n = 256)$y, grid)$y
  # Each density is clipped to the range of the values it describes
  inside <- function(v) if (length(v) == 0) FALSE else grid >= min(v) & grid <= max(v)
  all <- dens(x) * inside(x)
  det <- dens(x[detected]) * mean(detected) * inside(x[detected])
  scale <- HIST_HEIGHT * diff(yr) / max(all, 1e-12)
  tibble(x = grid, ymin = yr[1], y_all = yr[1] + all * scale, y_det = yr[1] + det * scale)
}

# y limits from the curves (not the confidence bands)
clip_limits <- function(est) {
  r <- range(est, na.rm = TRUE)
  mid <- mean(r); half <- max(diff(r) / 2, 0.25) * CI_CLIP_FACTOR
  c(mid - half, mid + half)
}

p_label <- function(p) if_else(p < 0.001, "< 0.001", paste("=", sprintf("%.3f", p)))

# -----------------------------------------------------------------------------
# 1. Sea-ice models (significant smooths only)
# -----------------------------------------------------------------------------
# Sea-ice models contain sea-ice concentration and, at KGI and CI, the signed
# distance to the ice edge (ICE_EDGE_VAR; decision N. Posdaljian 2026-11-20).
# One column per ice variable; concentration in %, distance in km (negative =
# recorder inside the ice).
ICE_VAR_LABEL <- c(ice_conc = "Sea-ice concentration (%)",
                   dist_ice_edge_km = "Distance to ice edge (km; < 0 = inside the ice)")
ice_models <- read_stage_models("sea_ice")
ice_terms <- bind_rows(map(ice_models, "terms")) |>
  select(species, site, variable, p_value)
to_plot_x <- \(x, v) if_else(v == "ice_conc", 100 * x, x)
ice_curves <- model_curves(ice_models) |>
  left_join(ice_terms, by = c("species", "site", "variable")) |>
  mutate(species = factor(species, levels = MODELED_SPECIES),
         site = factor(site, levels = SITE_ORDER))
write_csv(ice_curves, file.path(RESULTS_DIR, "sea_ice_response_curves_data.csv"))

ice_sig <- filter(ice_curves, p_value <= P_THRESH)
# Robustness to residual autocorrelation (05_ar1_sensitivity.R on the sea-ice
# models): solid = effect still significant with an AR(1) process where AIC
# supports one; dashed = not robust. Without the AR(1) table all are solid.
ar1_file <- file.path(STATS_DIR, "sea_ice", "tables", "ar1_sensitivity_terms.csv")
robust <- if (file.exists(ar1_file)) {
  read_csv(ar1_file, show_col_types = FALSE) |>
    transmute(species, site, variable = gsub("^s[(]|[)]$", "", term),
              robust = coalesce(holds, TRUE))
} else distinct(ice_sig, species, site, variable) |> mutate(robust = TRUE)
ice_sig <- ice_sig |>
  left_join(mutate(robust, species = factor(species, levels = MODELED_SPECIES),
                   site = factor(site, levels = SITE_ORDER)),
            by = c("species", "site", "variable")) |>
  mutate(robust = factor(if_else(coalesce(robust, TRUE), "Robust to autocorrelation (AR1)",
                                 "Not robust (AR1)"),
                         levels = c("Robust to autocorrelation (AR1)", "Not robust (AR1)")))
ice_rug <- model_rugs(ice_models) |>
  semi_join(distinct(ice_sig, species, site, variable), by = c("species", "site", "variable")) |>
  mutate(species = factor(species, levels = MODELED_SPECIES),
         site = factor(site, levels = SITE_ORDER))
# Where a species was never detected above some ice concentration (complete
# separation; the beaked whales), the fitted decline beyond that point has an
# arbitrary steepness. Such a concentration curve (at least MIN_BINS_ABOVE bins
# above the highest detection) therefore stops at the highest bin-mean ice
# concentration with a detection at that site (decision N. Posdaljian,
# 2026-11-19); the figure legend states the truncation.
MIN_BINS_ABOVE <- 5
max_ice_detected <- imap(ice_models, function(m, id) {
  tibble(species = factor(m$species, levels = MODELED_SPECIES),
         site = factor(m$site, levels = SITE_ORDER),
         max_ice_det = max(m$data$ice_conc[m$data[[m$species]] > 0]),
         n_bins_above = sum(m$data$ice_conc > max_ice_det))
}) |>
  list_rbind()
write_csv(max_ice_detected, file.path(RESULTS_DIR, "sea_ice_max_ice_with_detection.csv"))
ice_sig <- ice_sig |>
  inner_join(max_ice_detected, by = c("species", "site")) |>
  mutate(truncated = variable == "ice_conc" & n_bins_above >= MIN_BINS_ABOVE) |>
  filter(!truncated | x <= max_ice_det + 1e-9) |>
  # Re-center truncated curves over the range shown (a centered smooth has
  # no fixed level; the original centering included the removed part).
  group_by(species, site, variable) |>
  mutate(shift = if_else(truncated, mean(.estimate), 0),
         across(c(.estimate, .lower_ci, .upper_ci), ~ .x - shift)) |>
  ungroup() |>
  select(-shift)

ice_limits <- ice_sig |>
  group_by(species) |>
  summarise(lo = clip_limits(.estimate)[1], hi = clip_limits(.estimate)[2])

# One faceted plot (rows = species, columns = ice variable). Bands are clipped
# to each species' limits in the data, since facets cannot have separate
# coord_cartesian() limits; the clipping is stated in the figure legend.
ice_plot_data <- ice_sig |>
  left_join(ice_limits, by = "species") |>
  mutate(.lower_ci = pmax(.lower_ci, lo), .upper_ci = pmin(.upper_ci, hi),
         x_plot = to_plot_x(x, variable),
         variable = factor(ICE_VAR_LABEL[variable], ICE_VAR_LABEL))
ice_rug <- ice_rug |>
  mutate(x_plot = to_plot_x(x, variable), variable = factor(ICE_VAR_LABEL[variable], ICE_VAR_LABEL))
ice_species <- intersect(MODELED_SPECIES, as.character(unique(ice_sig$species)))
n_ice_cols <- n_distinct(ice_plot_data$variable)
ice_plot <- ggplot(ice_plot_data, aes(x_plot, .estimate, colour = site, fill = site)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = robust), linewidth = 0.8) +
  scale_linetype_manual(values = c("Robust to autocorrelation (AR1)" = "solid",
                                   "Not robust (AR1)" = "22"), name = NULL, drop = FALSE) +
  geom_rug(data = ice_rug, aes(x = x_plot, colour = site), inherit.aes = FALSE,
           sides = "b", alpha = 0.6, length = unit(0.04, "npc")) +
  facet_grid(species ~ variable, scales = "free", switch = "x",
             labeller = labeller(species = as_labeller(\(x) species_label(x)))) +
  scale_colour_manual(values = SITE_COLORS, drop = FALSE, name = "Site") +
  scale_fill_manual(values = SITE_COLORS, drop = FALSE, name = "Site") +
  labs(x = NULL, y = "Centered partial effect (log scale)") +
  theme_bw(base_size = 9) +
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 1)) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.placement = "outside", strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
        legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
ice_h <- 1.9 * length(ice_species) + 0.8
ice_w <- 2.6 * n_ice_cols + 1.6
ggsave(file.path(RESULTS_DIR, "supplementary_sea_ice_response_curves.png"), ice_plot,
       width = ice_w, height = ice_h, dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(RESULTS_DIR, "supplementary_sea_ice_response_curves.pdf"), ice_plot,
       width = ice_w, height = ice_h, device = cairo_pdf)

# -----------------------------------------------------------------------------
# 2. Habitat models: one figure per model
# -----------------------------------------------------------------------------
# Models reported as "no significant terms" have no curves to show.
hab_models <- read_stage_models(RESULTS_STAGE) |>
  keep(\(m) nrow(m$terms) > 0)
hab_curves <- model_curves(hab_models)
hab_rugs   <- model_rugs(hab_models)
write_csv(hab_curves, file.path(RESULTS_DIR,
                                paste0(HABITAT_LABEL, "_habitat_response_curves_data.csv")))

plot_model_curves <- function(m) {
  term_info <- m$terms |>
    mutate(panel = paste0(variable_label(variable), "\nedf = ", round(edf, 2),
                          ", p ", p_label(p_value),
                          ", dev. expl. = ", round(100 * dev_expl_term, 1), "%"))
  d <- hab_curves |>
    filter(species == m$species, site == m$site) |>
    left_join(term_info |> select(variable, panel), by = "variable")
  rug <- hab_rugs |>
    filter(species == m$species, site == m$site) |>
    left_join(term_info |> select(variable, panel), by = "variable")

  ggplot(d, aes(x, .estimate)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = SITE_COLORS[[m$site]], alpha = 0.2) +
    geom_line(colour = SITE_COLORS[[m$site]], linewidth = 0.8) +
    geom_rug(data = rug, aes(x = x), inherit.aes = FALSE, sides = "b",
             alpha = 0.5, length = unit(0.03, "npc")) +
    facet_wrap(~panel, scales = "free_x") +
    coord_cartesian(ylim = clip_limits(d$.estimate)) +
    labs(x = NULL, y = "Centered partial effect (log scale)",
         title = paste0(species_label(m$species), " at ", m$site,
                        if (isTRUE(m$meta$low_sample)) " (low sample)" else ""),
         subtitle = paste0("Deviance explained ", round(100 * m$summary$dev_expl, 1),
                           "%; ", m$summary$n_bins, " bins of ", m$summary$bin_days,
                           " d (", m$summary$n_detection_bins, " with detections)")) +
    theme_bw(base_size = 9) +
    theme(strip.background = element_blank(), strip.text = element_text(hjust = 0))
}

for (id in names(hab_models)) {
  n_terms <- nrow(hab_models[[id]]$terms)
  ggsave(file.path(hab_fig_dir, paste0(id, ".png")), plot_model_curves(hab_models[[id]]),
         width = 3.2 * min(n_terms, 3) + 0.8, height = 3 * ceiling(n_terms / 3) + 0.8,
         dpi = 300, bg = "white")
}

# -----------------------------------------------------------------------------
# 3. Publication grid: rows = species x site models, columns = predictors
# -----------------------------------------------------------------------------
# Columns are the predictors that appear in at least one final model, in
# PREDICTOR_ORDER; each column shares its x range (pooled observed range) so
# shapes compare across species and sites. Each row has its own y range
# (curves of one model share a scale). Empty cells are left blank.
model_order <- expand_grid(species = MODELED_SPECIES, site = SITE_ORDER) |>
  semi_join(distinct(hab_curves, species, site), by = c("species", "site"))
cols <- intersect(PREDICTOR_ORDER, unique(hab_curves$variable))
x_ranges <- hab_rugs |>
  group_by(variable) |>
  summarise(lo = min(x, na.rm = TRUE), hi = max(x, na.rm = TRUE))
row_limits <- hab_curves |>
  group_by(species, site) |>
  summarise(lo = clip_limits(.estimate)[1], hi = clip_limits(.estimate)[2], .groups = "drop")

term_p <- bind_rows(map(hab_models, "terms")) |> select(species, site, variable, p_value)

grid_cell <- function(sp, st, v, first_col, last_row) {
  d <- filter(hab_curves, species == sp, site == st, variable == v)
  xr <- unlist(x_ranges[x_ranges$variable == v, c("lo", "hi")])
  yr <- unlist(row_limits[row_limits$species == sp & row_limits$site == st, c("lo", "hi")])
  theme_cell <- theme_bw(base_size = 7) +
    theme(panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2),
          axis.title.y = if (first_col) element_text(size = 7, face = "bold",
                                                     colour = SITE_COLORS[[st]]) else element_blank(),
          axis.text.x = element_text(size = 5.5),
          axis.title.x = if (last_row) element_text(size = 6.5) else element_blank())
  ylab <- paste0(sp, " ", st)
  if (nrow(d) == 0) {
    return(ggplot() + coord_cartesian(xlim = xr, ylim = yr) +
             labs(x = PREDICTOR_SHORT[[v]], y = ylab) + theme_cell +
             theme(panel.border = element_blank(), axis.ticks = element_blank(),
                   axis.text = element_blank()))
  }
  p <- term_p$p_value[term_p$species == sp & term_p$site == st & term_p$variable == v]
  obs <- filter(hab_rugs, species == sp, site == st, variable == v)
  h <- density_curves(obs$x, obs$detected, xr, yr)
  ggplot(d, aes(x, .estimate)) +
    geom_ribbon(data = h, aes(x = x, ymin = ymin, ymax = y_all), inherit.aes = FALSE,
                fill = "grey82", colour = NA) +
    geom_ribbon(data = h, aes(x = x, ymin = ymin, ymax = y_det), inherit.aes = FALSE,
                fill = SITE_COLORS[[st]], alpha = 0.55, colour = NA) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = SITE_COLORS[[st]], alpha = 0.2) +
    geom_line(colour = SITE_COLORS[[st]], linewidth = 0.6) +
    annotate("text", x = xr[2], y = yr[2], hjust = 1.05, vjust = 1.3, size = 1.9,
             label = paste0("p ", p_label(p))) +
    coord_cartesian(xlim = xr, ylim = yr) +
    labs(x = PREDICTOR_SHORT[[v]], y = ylab) +
    theme_cell
}


# Writes one grid for the given predictor columns and model rows.
save_grid <- function(cols, rows, file_stem) {
  cells <- list()
  for (r in seq_len(nrow(rows))) {
    for (j in seq_along(cols)) {
      cells[[length(cells) + 1]] <- grid_cell(rows$species[r], rows$site[r], cols[j],
                                              first_col = j == 1,
                                              last_row = r == nrow(rows))
    }
  }
  fig <- wrap_plots(cells, ncol = length(cols), byrow = TRUE)
  w <- 1.45 * length(cols) + 0.5
  h <- 1.05 * nrow(rows) + 0.6
  ggsave(file.path(RESULTS_DIR, paste0(file_stem, ".png")), fig, width = w, height = h,
         dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(RESULTS_DIR, paste0(file_stem, ".pdf")), fig, width = w, height = h,
         device = cairo_pdf)
}

# Main text: models with at least MIN_TERMS_MAIN retained terms, with all of
# their terms (decision N. Posdaljian 2026-11-20: the main figure shows the
# models that can be described ecologically; single-term models go to the
# supplementary grid).
MIN_TERMS_MAIN <- 2
model_terms <- distinct(hab_curves, species, site, variable)
n_terms <- count(model_terms, species, site)
split_grid <- function(rows, file_stem) {
  if (nrow(rows) == 0) return(invisible(NULL))
  used <- semi_join(model_terms, rows, by = c("species", "site"))$variable
  save_grid(intersect(cols, used), rows, file_stem)
}
split_grid(semi_join(model_order, filter(n_terms, n >= MIN_TERMS_MAIN), by = c("species", "site")),
           paste0(HABITAT_LABEL, "_habitat_response_curves_grid"))
split_grid(semi_join(model_order, filter(n_terms, n < MIN_TERMS_MAIN), by = c("species", "site")),
           paste0("supplementary_", RESULTS_STAGE, "_habitat_response_curves_single_term"))

message("Response curves (sea ice + habitat ", RESULTS_STAGE, ") written to ", RESULTS_DIR)
