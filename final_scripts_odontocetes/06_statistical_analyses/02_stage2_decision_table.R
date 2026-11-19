# =============================================================================
# Script  : 02_stage2_decision_table.R
# Section : 2.6 Statistical analyses (Stage 2: input to manual harmonization)
# Purpose : Assemble the Stage 1 evidence into one decision table per
#           species x site x candidate variable, plus per-theme and
#           per-variable summaries and at-a-glance figures, so that one
#           variable per theme can be chosen by hand for the harmonized
#           Stage 3 models.
# Inputs  : <OUTPUT_ROOT>/06_statistical_analyses/stage1/tables/
#             model_inclusion.csv, candidates.csv, theme_selection.csv,
#             single_term_smooths.csv, final_model_terms.csv,
#             final_model_summary.csv, seasonal_variance.csv
# Outputs : <OUTPUT_ROOT>/06_statistical_analyses/stage2/
#             stage2_decision_table.csv  (one row per species x site x variable)
#             stage2_theme_summary.csv   (per species x theme and all species)
#             stage2_variable_summary.csv (per species x variable across sites)
#             stage2_selection_heatmap.png/.pdf
#             candidate_curves/<species>_<site>_single_term_curves.png
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Column definitions:
#             screening        - removed (missing values / zero variance / VIF)
#                                or VIF retained
#             competes         - entered the within-theme AIC comparison
#                                (VIF retained and single-term fit OK)
#             version          - raw (always), raw or anomaly; raw vs anomaly
#                                decided per site x variable from
#                                seasonal_var_frac (> SEASONAL_VAR_THRESH)
#             single_*         - single-term model y ~ s(x, k = 4)
#             delta_aic_theme  - AIC minus the best AIC among competitors in
#                                the theme (0 = theme winner)
#             final_*          - term in the Stage 1 final model
#             shape/direction  - classify_smooth_shape() rules (00_gam_functions.R:
#                                "not significant" if p > 0.05, else shape of
#                                the fitted centered curve), evaluated
#                                on the final smooth if selected, otherwise on
#                                the single-term smooth (shape_source says which)
#             model_status     - inclusion status; "no significant terms" when
#                                the last term left in backward elimination
#                                had p > 0.05
#           Run after 01_stage1_automated_selection.R.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

s1_tables <- file.path(STATS_DIR, "stage1", "tables")
out_dir   <- stage_dir("stage2")
curve_dir <- stage_dir("stage2", "candidate_curves")

read_s1 <- function(file) read_csv(file.path(s1_tables, file), show_col_types = FALSE)
inclusion   <- read_s1("model_inclusion.csv")
candidates  <- read_s1("candidates.csv")
theme_sel   <- read_s1("theme_selection.csv")
smooths     <- read_s1("single_term_smooths.csv")
final_terms <- read_s1("final_model_terms.csv")
final_summ  <- read_s1("final_model_summary.csv")
seasonal    <- read_s1("seasonal_variance.csv")

# -----------------------------------------------------------------------------
# 1. Decision table
# -----------------------------------------------------------------------------
# Start from every candidate of every modeled species x site (so variables
# removed before screening stay visible), then add single-term evidence,
# theme ranking and the final-model results.
single <- theme_sel |>
  group_by(species, site, theme) |>
  mutate(delta_aic_theme = if (any(competes)) AIC - min(AIC[competes], na.rm = TRUE)
                           else NA_real_) |>
  ungroup() |>
  transmute(species, site, variable, competes,
            single_p = p_value, single_dev_expl = dev_expl, single_edf = edf,
            single_aic = AIC, delta_aic_theme = if_else(competes, delta_aic_theme, NA_real_),
            aic_rank_theme, aic_rank_all, theme_winner,
            single_shape = shape, single_direction = direction)

final <- final_terms |>
  transmute(species, site, variable, selected_final = TRUE,
            final_p = p_value, final_dev_expl_term = dev_expl_term, final_edf = edf,
            final_shape = shape, final_direction = direction)

# raw vs anomaly per site x variable (environmental data only).
version_tbl <- seasonal |>
  transmute(site, source = variable, seasonal_var_frac)

decision <- candidates |>
  left_join(inclusion |> select(species, site, model_status = status),
            by = c("species", "site")) |>
  left_join(final_summ |> select(species, site, no_significant_terms),
            by = c("species", "site")) |>
  mutate(model_status = if_else(coalesce(no_significant_terms, FALSE),
                                "no significant terms", model_status),
         source = sub("_anom$", "", variable)) |>
  left_join(version_tbl, by = c("site", "source")) |>
  mutate(version = case_when(is_anom ~ "anomaly",
                             !is.na(seasonal_var_frac) ~ "raw",
                             TRUE ~ "raw (always)")) |>
  left_join(single, by = c("species", "site", "variable")) |>
  left_join(final,  by = c("species", "site", "variable")) |>
  mutate(variable_label = variable_label(variable),
         vif_retained = screening == "VIF retained",
         competes = coalesce(competes, FALSE),
         theme_winner = coalesce(theme_winner, FALSE),
         selected_final = coalesce(selected_final, FALSE),
         shape = if_else(selected_final, final_shape, single_shape),
         direction = if_else(selected_final, final_direction, single_direction),
         shape_source = case_when(selected_final ~ "final model",
                                  !is.na(single_shape) ~ "single-term model",
                                  TRUE ~ NA_character_),
         theme = factor(theme, levels = THEME_ORDER),
         species = factor(species, levels = MODELED_SPECIES),
         site = factor(site, levels = SITE_ORDER)) |>
  arrange(species, site, theme, base, depth_m, lag_months, is_anom) |>
  select(species, site, model_status, theme, variable, variable_label, base, depth_m,
         version, seasonal_var_frac,
         lag_months, is_anom, screening, vif_retained, competes,
         single_p, single_dev_expl, single_edf, single_aic, delta_aic_theme,
         aic_rank_theme, aic_rank_all, theme_winner,
         selected_final, final_p, final_dev_expl_term, final_edf,
         direction, shape, shape_source, single_direction, single_shape)

write_csv(decision, file.path(out_dir, "stage2_decision_table.csv"))

# -----------------------------------------------------------------------------
# 2. Summaries
# -----------------------------------------------------------------------------
# Per theme: how often the theme contributed a term to a final model, which
# variables did, and mean deviance explained (single-term for theme winners,
# drop-one for selected terms). Rows with species = "ALL" pool all species.
summarise_theme <- function(df) {
  df |>
    summarise(
      n_models = n_distinct(paste(species, site)),
      n_models_with_winner = n_distinct(paste(species, site)[theme_winner]),
      n_models_selected = n_distinct(paste(species, site)[selected_final]),
      selected_variables = paste(sort(unique(variable[selected_final])), collapse = "; "),
      selected_directions = paste(direction[selected_final], collapse = "; "),
      mean_single_dev_expl_winner = mean(single_dev_expl[theme_winner], na.rm = TRUE),
      mean_final_dev_expl_term = mean(final_dev_expl_term[selected_final], na.rm = TRUE),
      .groups = "drop")
}
theme_summary <- bind_rows(
  decision |> group_by(species, theme) |> summarise_theme() |>
    mutate(species = as.character(species)),
  decision |> group_by(theme) |> summarise_theme() |> mutate(species = "ALL")) |>
  mutate(across(where(is.numeric), \(x) if_else(is.nan(x), NA_real_, x))) |>
  arrange(factor(species, levels = c(MODELED_SPECIES, "ALL")), theme)
write_csv(theme_summary, file.path(out_dir, "stage2_theme_summary.csv"))

# Per variable within species (across sites): consistency of evidence.
variable_summary <- decision |>
  group_by(species, theme, variable, variable_label) |>
  summarise(
    n_sites_candidate = n(),
    n_sites_vif_retained = sum(vif_retained),
    n_sites_theme_winner = sum(theme_winner),
    n_sites_selected = sum(selected_final),
    n_sites_single_p_sig = sum(single_p <= P_THRESH, na.rm = TRUE),
    mean_single_dev_expl = mean(single_dev_expl, na.rm = TRUE),
    directions_by_site = paste0(site, ":", coalesce(direction, "NA"), collapse = "; "),
    shapes_by_site = paste0(site, ":", coalesce(shape, "NA"), collapse = "; "),
    .groups = "drop") |>
  mutate(mean_single_dev_expl = if_else(is.nan(mean_single_dev_expl), NA_real_,
                                        mean_single_dev_expl)) |>
  arrange(species, theme, desc(n_sites_selected), desc(mean_single_dev_expl))
write_csv(variable_summary, file.path(out_dir, "stage2_variable_summary.csv"))

# -----------------------------------------------------------------------------
# 3. At-a-glance heatmap
# -----------------------------------------------------------------------------
# One column per site within each species panel; one row per candidate
# variable grouped by theme. Fill = selection status (direction for selected
# terms); text = shape code, so identity never depends on colour alone.
status_levels <- c("Selected: positive", "Selected: negative", "Selected: mixed / cyclic",
                   "Selected: flat", "Theme winner, dropped", "Competed",
                   "Removed (VIF)", "Removed (no variance / NA)", "Not modeled")
status_colors <- c("#e34948", "#2a78d6", "#6250d6", "#77776f", "#a9a8a2", "#d6d5d0",
                   "#f1f0ec", "#ffffff", "#ffffff")
shape_codes <- c(increasing = "inc", decreasing = "dec", peaked = "pk",
                 `not significant` = "ns",
                 `U-shaped` = "U", complex = "cx", flat = "flat")

# Variables are ordered by theme, then name, and every species x site x
# variable cell is shown (sites without a model are marked "Not modeled").
heat <- decision |>
  mutate(status = case_when(
    selected_final ~ paste0("Selected: ", case_when(direction == "none" ~ "flat",
                                                    direction %in% c("cyclic", "mixed") ~ "mixed / cyclic",
                                                    TRUE ~ direction)),
    theme_winner ~ "Theme winner, dropped",
    competes ~ "Competed",
    screening == "removed (VIF)" ~ "Removed (VIF)",
    TRUE ~ "Removed (no variance / NA)"),
    label = case_when(selected_final & startsWith(shape, "cyclic") ~ "cyc",
                      selected_final ~ unname(shape_codes[shape]),
                      theme_winner ~ "w",
                      grepl("variance|missing", screening) ~ "nv",
                      TRUE ~ "")) |>
  select(species, site, theme, variable, status, label)
heat <- heat |>
  distinct(species, theme, variable) |>
  crossing(site = factor(SITE_ORDER, levels = SITE_ORDER)) |>
  left_join(heat, by = c("species", "theme", "variable", "site")) |>
  mutate(status = factor(coalesce(status, "Not modeled"), levels = status_levels),
         label = coalesce(label, ""),
         row = paste(theme, variable),
         variable_label = variable_label(variable)) |>
  arrange(theme, variable) |>
  mutate(row = factor(row, levels = rev(unique(row))))

heatmap_plot <- ggplot(heat, aes(site, row, fill = status)) +
  geom_tile(colour = "grey80", linewidth = 0.3) +
  geom_text(aes(label = label), size = 2.2,
            colour = if_else(grepl("^Selected", heat$status), "white", "grey20")) +
  facet_wrap(~species, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = setNames(status_colors, status_levels), drop = FALSE,
                    name = NULL) +
  scale_y_discrete(labels = \(x) variable_label(sub("^\\S+ ", "", x))) +
  labs(x = NULL, y = NULL,
       title = "Stage 1 predictor selection by species and site",
       subtitle = paste("Text: shape of selected smooth (inc, dec, pk = peaked, U, cx = complex, cyc = cyclic,",
                        "flat); w = theme winner removed in backward elimination;",
                        "nv = no variance / missing values; blank column = not modeled")) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        axis.text.y = element_text(size = 5.5), strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "stage2_selection_heatmap.png"), heatmap_plot,
       width = 16, height = 11, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "stage2_selection_heatmap.pdf"), heatmap_plot,
       width = 16, height = 11, bg = "white")

# -----------------------------------------------------------------------------
# 4. Single-term curves of all candidates (one figure per species x site)
# -----------------------------------------------------------------------------
# Centered single-term partial effects (log scale) let the shapes of
# alternative variables within a theme be compared before choosing one.
curve_status <- decision |>
  transmute(species = as.character(species), site = as.character(site), variable, theme,
            status = case_when(selected_final ~ "selected in final",
                               theme_winner ~ "theme winner",
                               competes ~ "competed",
                               TRUE ~ "removed (VIF)"))
curves <- smooths |>
  left_join(curve_status, by = c("species", "site", "variable")) |>
  mutate(panel = paste0(theme, ": ", variable_label(variable)),
         status = factor(status, levels = c("selected in final", "theme winner",
                                            "competed", "removed (VIF)")))

for (id in unique(paste(curves$species, curves$site, sep = "_"))) {
  d <- curves |> filter(paste(species, site, sep = "_") == id) |>
    arrange(factor(theme, levels = THEME_ORDER), variable) |>
    mutate(panel = factor(panel, levels = unique(panel)))
  p <- ggplot(d, aes(x, .estimate)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = status), alpha = 0.25) +
    geom_line(aes(colour = status), linewidth = 0.7) +
    facet_wrap(~panel, scales = "free", labeller = label_wrap_gen(32)) +
    scale_colour_manual(values = c("#e34948", "#2a78d6", "grey35", "grey65"), drop = FALSE) +
    scale_fill_manual(values = c("#e34948", "#2a78d6", "grey35", "grey65"), drop = FALSE) +
    labs(x = NULL, y = "Centered partial effect (log scale)", colour = NULL, fill = NULL,
         title = paste0(species_label(d$species[1]), " at ", d$site[1],
                        ": single-term models of all candidates")) +
    theme_bw(base_size = 7) +
    theme(legend.position = "bottom", strip.text = element_text(size = 5.5))
  ggsave(file.path(curve_dir, paste0(id, "_single_term_curves.png")), p,
         width = 14, height = 11, dpi = 200, bg = "white")
}

message("Stage 2 complete. Outputs in ", out_dir)
