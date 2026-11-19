# =============================================================================
# Script  : 06_compare_models_insitu.R
# Section : Supplementary Fig. S4 (2.5.1) - model validation against in-situ
#           profiles
# Purpose : Compare Copernicus and HYCOM with in-situ temperature, salinity
#           and (Copernicus only) dissolved-oxygen profiles from the World
#           Ocean Database in the wider study region, 2014-2016: model minus
#           observation at every model level to 1000 m and for the three
#           depth layers used by the habitat models, on the same casts and
#           levels for both models.
# Inputs  : outputs/05b_hycom_comparison/insitu_profiles.csv
#             (04_download_insitu_profiles.R)
#           DOWNLOAD_DIR/models_at_profiles/<model>/<date>.rds
#             (05_download_models_at_profiles.R; read directly, so the script
#             can be run while HYCOM days are still downloading)
#           config.R (SITES, LAYERS); 01_study_area/map_helpers.R (GMRT map)
# Outputs : outputs/05b_hycom_comparison/
#             insitu_matchups.csv        cast x level x variable: obs and
#                                        model values (both models)
#             insitu_metrics_by_level.csv model, variable, level: n, bias,
#                                        RMSE, r
#             insitu_metrics_by_layer.csv model, variable, layer: n casts,
#                                        bias, RMSE, r of layer means
#             FigS4_insitu_validation.png / .pdf
#               (a) map of the profiles by instrument type with the site
#                   boxes; (b) bias and RMSE profiles per variable and model
# Author  : Natalie Posdaljian
# Date    : 2026-11-19
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Quality control: only observations with WOD observation flag 0
#             AND profile flag 0 (accepted) are used.
#           - Vertical matching: each observed profile is linearly
#             interpolated to the model levels inside its observed depth
#             range (no extrapolation; casts with gaps > MAX_GAP_M between
#             observations around a level are not used at that level).
#           - Oxygen: WOD reports umol kg-1; converted to mmol m-3 with the
#             in-situ density from the observed T, S and depth (oce, EOS-80).
#           - Fair comparison: metrics use only cast x level x variable
#             combinations available in BOTH models (oxygen: Copernicus only).
#           - Layer means are thickness-weighted over the levels of each
#             layer, computed for a cast only if it covers at least
#             MIN_LAYER_COVERAGE of the layer's thickness.
#           - Both models assimilate some in-situ profiles (e.g. Argo), so
#             this is not a fully independent validation.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "01_study_area", "map_helpers.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggnewscale)
  library(oce)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
SECTION     <- "05b_hycom_comparison"
CACHE_DIR   <- file.path(DOWNLOAD_DIR, "models_at_profiles")
MAX_DEPTH_M <- max(LAYERS$bottom_m)
MAX_GAP_M   <- 50       # max distance between observations bracketing a level
MIN_LAYER_COVERAGE <- 0.8
REGION_LON  <- c(-60, -52)
REGION_LAT  <- c(-63, -60)
PROBE_LABELS <- c(apb = "Animal-borne (seal tag)", pfl = "Argo float",
                  ctd = "Ship CTD", osd = "Bottle / low-res. CTD",
                  xbt = "XBT (temperature only)")
MODEL_COLORS <- c(Copernicus = "#0072B2", HYCOM = "#D55E00")
VAR_LABELS   <- c(temperature = "Temperature (°C)", salinity = "Salinity",
                  oxygen = "Oxygen (mmol m⁻³)")

# -----------------------------------------------------------------------------
# 2. Observations (quality-controlled)
# -----------------------------------------------------------------------------
obs_raw <- read_csv(out_path(SECTION, "insitu_profiles.csv"), show_col_types = FALSE,
                    col_types = cols(platform = col_character(), .default = col_guess()))

casts <- obs_raw |> distinct(cast_id, probe, time, lat, lon)

obs <- obs_raw |>
  filter(depth_m <= MAX_DEPTH_M + 100) |>
  mutate(
    temperature = if_else(temperature_flag == 0 & coalesce(temperature_profile_flag, 0) == 0,
                          temperature, NA_real_),
    salinity    = if_else(salinity_flag == 0 & coalesce(salinity_profile_flag, 0) == 0,
                          salinity, NA_real_),
    oxygen_umolkg = if_else(oxygen_flag == 0 & coalesce(oxygen_profile_flag, 0) == 0,
                            oxygen, NA_real_),
    # umol kg-1 -> mmol m-3: x in-situ density (kg m-3) / 1000
    rho = swRho(salinity, temperature, pressure = depth_m, eos = "unesco"),
    oxygen = oxygen_umolkg * rho / 1000
  ) |>
  select(cast_id, depth_m, temperature, salinity, oxygen) |>
  pivot_longer(c(temperature, salinity, oxygen), names_to = "variable",
               values_to = "obs") |>
  filter(!is.na(obs))

# -----------------------------------------------------------------------------
# 3. Model profiles at the casts
# -----------------------------------------------------------------------------
models <- list.files(CACHE_DIR, pattern = "[.]rds$", recursive = TRUE,
                     full.names = TRUE) |>
  map(readRDS) |>
  list_rbind() |>
  filter(model_depth_m <= MAX_DEPTH_M)
message("Casts with model data: ",
        paste(names(table(distinct(models, model, cast_id)$model)),
              table(distinct(models, model, cast_id)$model), collapse = ", "))

# Copernicus model levels are the reference vertical grid for both models:
# HYCOM is interpolated to them (as in 03_compare_hycom_copernicus.R).
cop_levels <- models |> filter(model == "Copernicus") |> distinct(variable, model_depth_m)

interp_at <- function(z, v, levels, max_gap = Inf) {
  ok <- !is.na(v); z <- z[ok]; v <- v[ok]
  if (length(z) < 2) return(rep(NA_real_, length(levels)))
  o <- order(z); z <- z[o]; v <- v[o]
  out <- approx(z, v, xout = levels, rule = 1, ties = mean)$y
  # gap check: distance between the bracketing observations
  i <- findInterval(levels, z)
  gap <- ifelse(i >= 1 & i < length(z), z[pmin(i + 1, length(z))] - z[pmax(i, 1)], Inf)
  gap[levels %in% z] <- 0
  out[gap > max_gap] <- NA
  out
}

level_grid <- cop_levels |> group_by(variable) |>
  summarise(levels = list(sort(model_depth_m)), .groups = "drop")

obs_at_levels <- obs |>
  group_by(cast_id, variable) |>
  group_modify(function(d, key) {
    lv <- level_grid$levels[[match(key$variable, level_grid$variable)]]
    tibble(model_depth_m = lv, obs = interp_at(d$depth_m, d$obs, lv, MAX_GAP_M))
  }) |>
  ungroup() |>
  filter(!is.na(obs))

model_at_levels <- models |>
  group_by(model, cast_id, variable) |>
  group_modify(function(d, key) {
    lv <- level_grid$levels[[match(key$variable, level_grid$variable)]]
    tibble(model_depth_m = lv, value = interp_at(d$model_depth_m, d$value, lv))
  }) |>
  ungroup() |>
  filter(!is.na(value)) |>
  pivot_wider(names_from = model, values_from = value)

matchups <- obs_at_levels |>
  inner_join(model_at_levels, by = c("cast_id", "variable", "model_depth_m")) |>
  left_join(casts, by = "cast_id")

# Same casts and levels for both models (oxygen: Copernicus only)
both <- matchups |>
  filter(!is.na(Copernicus), variable == "oxygen" | !is.na(HYCOM %||% NA_real_))
write_csv(both, out_path(SECTION, "insitu_matchups.csv"))

# -----------------------------------------------------------------------------
# 4. Metrics by level and by layer
# -----------------------------------------------------------------------------
long <- both |>
  pivot_longer(any_of(names(MODEL_COLORS)), names_to = "model", values_to = "pred") |>
  filter(!is.na(pred))

metrics <- function(d) {
  summarise(d, n = n(), bias = mean(pred - obs), rmse = sqrt(mean((pred - obs)^2)),
            r = if (n() >= 5) cor(pred, obs) else NA_real_, .groups = "drop")
}

by_level <- long |> group_by(model, variable, model_depth_m) |> metrics()
write_csv(by_level, out_path(SECTION, "insitu_metrics_by_level.csv"))

# Layer means per cast: thickness-weighted over the levels in each layer
layer_weights <- function(levels) {
  edges <- c(0, (head(levels, -1) + levels[-1]) / 2, Inf)
  expand_grid(i = seq_along(levels), l = seq_len(nrow(LAYERS))) |>
    mutate(model_depth_m = levels[i], layer = LAYERS$layer[l],
           weight = pmax(0, pmin(edges[i + 1], LAYERS$bottom_m[l]) -
                            pmax(edges[i], LAYERS$top_m[l]))) |>
    filter(weight > 0) |>
    select(model_depth_m, layer, weight)
}
weights <- level_grid |>
  mutate(w = map(levels, layer_weights)) |>
  select(variable, w) |>
  unnest(w)
layer_thickness <- set_names(LAYERS$bottom_m - LAYERS$top_m, LAYERS$layer)

by_layer_cast <- long |>
  inner_join(weights, by = c("variable", "model_depth_m"), relationship = "many-to-many") |>
  group_by(model, variable, cast_id, layer) |>
  summarise(coverage = sum(weight) / layer_thickness[first(layer)],
            obs = sum(obs * weight) / sum(weight),
            pred = sum(pred * weight) / sum(weight), .groups = "drop") |>
  filter(coverage >= MIN_LAYER_COVERAGE)

by_layer <- by_layer_cast |>
  group_by(model, variable, layer = factor(layer, levels = LAYERS$layer)) |>
  metrics() |>
  arrange(variable, layer, model)
write_csv(by_layer, out_path(SECTION, "insitu_metrics_by_layer.csv"))
print(by_layer, n = Inf)

# -----------------------------------------------------------------------------
# 5. Figure S4
# -----------------------------------------------------------------------------
# (a) Map of the casts used (at least one matched level), by instrument type
used_casts <- casts |>
  semi_join(both, by = "cast_id") |>
  mutate(probe = factor(PROBE_LABELS[probe], levels = PROBE_LABELS))

boxes <- SITES |>
  transmute(site, xmin = box_lon_min, xmax = box_lon_max,
            ymin = box_lat_min, ymax = box_lat_max)

bathy <- load_bathymetry(xlim = REGION_LON, ylim = REGION_LAT)
map_panel <- ggplot() +
  bathymetry_layers(bathy, ocean_alpha = 0.6) +
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "black", linewidth = 0.6) +
  site_layers() +
  # site_layers() sets the site colour scale; the profile types get their own
  new_scale_colour() +
  geom_point(data = used_casts, aes(lon, lat, colour = probe, shape = probe),
             size = 1.4, alpha = 0.8) +
  scale_colour_manual(values = c("#e66101", "#5e3c99", "#1b7837", "#b2182b", "grey20"),
                      name = "Profile type", drop = TRUE) +
  scale_shape_manual(values = c(16, 17, 15, 18, 4), name = "Profile type", drop = TRUE) +
  map_frame(xlim = REGION_LON, ylim = REGION_LAT, corner = "bl") +
  # Legends below the map: horizontal depth bar and profile types in two rows
  guides(fill = guide_colourbar(direction = "horizontal", title.position = "top",
                                barwidth = unit(5, "cm"), barheight = unit(0.35, "cm"),
                                order = 1),
         colour = guide_legend(nrow = 3, title.position = "top", order = 2),
         shape = guide_legend(nrow = 3, title.position = "top", order = 2)) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8))

# (b) Bias and RMSE profiles: one row per metric, each with its own x scales
level_plot_data <- by_level |>
  filter(n >= 10) |>
  mutate(variable = factor(VAR_LABELS[variable], levels = VAR_LABELS))

sqrt_rev <- scales::transform_compose(scales::transform_sqrt(),
                                      scales::transform_reverse())
metric_row <- function(metric, xlab, zero_line) {
  ggplot(level_plot_data, aes(.data[[metric]], model_depth_m, colour = model)) +
    { if (zero_line) geom_vline(xintercept = 0, colour = "grey50") } +
    geom_hline(yintercept = LAYERS$bottom_m[-nrow(LAYERS)], linetype = "dashed",
               colour = "grey40") +
    geom_path() +
    geom_point(size = 0.9) +
    facet_wrap(~variable, nrow = 1, scales = "free_x") +
    scale_y_continuous(transform = sqrt_rev, breaks = c(0, 50, 100, 200, 500, 1000)) +
    scale_colour_manual(values = MODEL_COLORS, name = NULL) +
    labs(x = xlab, y = "Depth (m)") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"))
}
profile_panel <- (metric_row("bias", "Bias (model − observed)", TRUE) /
                    metric_row("rmse", "RMSE", FALSE)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

fig <- wrap_elements(map_panel) + wrap_elements(profile_panel) +
  plot_layout(widths = c(1.15, 1.5)) +
  plot_annotation(tag_levels = "a")

ggsave(out_path(SECTION, "FigS4_insitu_validation.png"), fig, width = 14, height = 8,
       dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(out_path(SECTION, "FigS4_insitu_validation.pdf"), fig, width = 14, height = 8,
       device = cairo_pdf)
