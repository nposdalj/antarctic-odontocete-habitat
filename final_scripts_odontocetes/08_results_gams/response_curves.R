# =============================================================================
# Script  : response_curves.R
# Section : 3.2 GAM results
# Purpose : Partial-effect (response) curves: CENTERED smooths on the log
#           (link) scale with 95% pointwise confidence intervals.
#           Sea-ice models (main text): one figure, species in rows, sites as
#           colours, non-significant curves dashed.
#           Habitat models (main text = Stage 3; Stage 1 = supplementary): one
#           multi-panel figure per species x site with a rug of the binned
#           observations.
#           The plotted values are saved as tables for the final layout.
# Inputs  : <OUTPUT_ROOT>/06_statistical_analyses/sea_ice/models/*.rds
#           <OUTPUT_ROOT>/06_statistical_analyses/<RESULTS_STAGE>/models/*.rds
# Outputs : <OUTPUT_ROOT>/08_results_gams/
#             main_sea_ice_response_curves.png, sea_ice_response_curves_data.csv
#             <label>_habitat_response_curves/<species>_<site>.png
#             <label>_habitat_response_curves_data.csv
#           <label> = "main_stage3" or "supplementary_stage1"
# Author  : Natalie Posdaljian
# Date    : 2026-11-16 (updated 2026-11-17)
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Curves are centered partial effects (gratia::smooth_estimates),
#           NOT back-transformed: the original pipeline added the intercept
#           and applied plogis(), a logit inverse that is wrong for the
#           Tweedie log link. A value of 0 = average log presence minutes;
#           +0.69 = twice the average. Habitat models reported as "no
#           significant terms" have no curves. FSLE orientation is a cyclic
#           smooth plotted over -90..90 deg; both ends are the same direction
#           (the curve wraps).
#           The publication layout of the habitat curves is decided after
#           Stage 2: [INSERT layout] (e.g., species rows x theme columns).
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

# --- Parameters ---------------------------------------------------------------
RESULTS_STAGE <- Sys.getenv("ODONTOCETE_RESULTS_STAGE", unset = "stage3")  # or "stage1"
HABITAT_LABEL <- if (RESULTS_STAGE == "stage3") "main_stage3" else
  paste0("supplementary_", RESULTS_STAGE)
CURVE_GRID_N  <- 200   # evaluation points per smooth

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
hab_fig_dir <- file.path(RESULTS_DIR, paste0(HABITAT_LABEL, "_habitat_response_curves"))
if (!dir.exists(hab_fig_dir)) dir.create(hab_fig_dir, recursive = TRUE)

# Evaluate every reported smooth of a list of models (one row per grid point).
model_curves <- function(models) {
  imap(models, function(m, id) {
    map(m$terms$variable, \(v) smooth_grid(m$model, v, n = CURVE_GRID_N) |>
          mutate(variable = v)) |>
      bind_rows() |>
      mutate(species = m$species, site = m$site, .before = 1)
  }) |> bind_rows()
}

p_label <- function(p) if_else(p < 0.001, "< 0.001", paste("=", sprintf("%.3f", p)))

# -----------------------------------------------------------------------------
# 1. Sea-ice models (main text)
# -----------------------------------------------------------------------------
ice_models <- read_stage_models("sea_ice")
ice_curves <- model_curves(ice_models) |>
  left_join(bind_rows(map(ice_models, "terms")) |> select(species, site, p_value),
            by = c("species", "site")) |>
  mutate(significant = if_else(p_value <= P_THRESH, "p ≤ 0.05", "p > 0.05"),
         species = factor(species, levels = MODELED_SPECIES),
         site = factor(site, levels = SITE_ORDER))
write_csv(ice_curves, file.path(RESULTS_DIR, "sea_ice_response_curves_data.csv"))

ice_plot <- ggplot(ice_curves, aes(x, .estimate, colour = site, fill = site)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = significant), linewidth = 0.8) +
  facet_wrap(~species, ncol = 1, scales = "free_y",
             labeller = as_labeller(\(x) species_label(x))) +
  scale_colour_manual(values = SITE_COLORS, drop = FALSE) +
  scale_fill_manual(values = SITE_COLORS, drop = FALSE) +
  scale_linetype_manual(values = c("p ≤ 0.05" = "solid", "p > 0.05" = "dashed")) +
  labs(x = variable_label(ICE_MODEL_VAR), y = "Centered partial effect (log scale)",
       colour = "Site", fill = "Site", linetype = NULL) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold"))
ggsave(file.path(RESULTS_DIR, "main_sea_ice_response_curves.png"), ice_plot,
       width = 5, height = 1.8 * n_distinct(ice_curves$species) + 0.8, dpi = 300, bg = "white")

# -----------------------------------------------------------------------------
# 2. Habitat models: standard per-model panels
# -----------------------------------------------------------------------------
# Models reported as "no significant terms" have no curves to show.
hab_models <- read_stage_models(RESULTS_STAGE) |>
  keep(\(m) nrow(m$terms) > 0)
hab_curves <- model_curves(hab_models)
write_csv(hab_curves, file.path(RESULTS_DIR,
                                paste0(HABITAT_LABEL, "_habitat_response_curves_data.csv")))

# One panel per term (common y axis within a model so effect sizes compare);
# panel titles give edf, p and per-term deviance explained.
plot_model_curves <- function(m) {
  term_info <- m$terms |>
    mutate(panel = paste0(variable_label(variable), "\nedf = ", round(edf, 2),
                          ", p ", p_label(p_value),
                          ", dev. expl. = ", round(100 * dev_expl_term, 1), "%"))
  d <- hab_curves |>
    filter(species == m$species, site == m$site) |>
    left_join(term_info |> select(variable, panel), by = "variable")
  rug <- m$terms$variable |>
    map(\(v) tibble(variable = v, x = m$data[[v]])) |>
    bind_rows() |>
    left_join(term_info |> select(variable, panel), by = "variable")

  ggplot(d, aes(x, .estimate)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#2a78d6", alpha = 0.2) +
    geom_line(colour = "#2a78d6", linewidth = 0.8) +
    geom_rug(data = rug, aes(x = x), inherit.aes = FALSE, sides = "b",
             alpha = 0.5, length = unit(0.03, "npc")) +
    facet_wrap(~panel, scales = "free_x") +
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
# 3. Publication layout of the habitat curves
# -----------------------------------------------------------------------------
# [INSERT layout] after Stage 2/3: combine the curves in
# <label>_habitat_response_curves_data.csv into the manuscript figure.

message("Response curves (sea ice + habitat ", RESULTS_STAGE, ") written to ", RESULTS_DIR)
