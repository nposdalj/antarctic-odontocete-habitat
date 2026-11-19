# =============================================================================
# Script  : 06_box_size_comparison.R
# Section : 05_environmental_data
# Purpose : Compare the daily environmental predictors computed over the MAIN
#           +/- 20 km box with those over the +/- 40 km comparison box, per
#           site and predictor: Pearson correlation, mean difference and RMSE,
#           over the recording periods and over the whole climatology window.
#           Supports a one-sentence statement in the supplement on how
#           sensitive the predictors are to the box size.
# Inputs  : final_scripts_odontocetes/config.R
#           Main box (04_extract_copernicus.R, 03_extract_aviso_fsle.R):
#             F_COP_PHY_DEPTH, F_COP_BGC_DEPTH, F_COP_SURFACE, F_FSLE_DAILY
#           Comparison box: the same files with suffix "_box40"
# Outputs : outputs/05_environmental_data/box_size_comparison.csv
#             site, variable, period (recording | climatology), n_days, r,
#             mean_diff (main - box40), rmse, sd_box40, rmse_rel (rmse /
#             sd_box40)
#           outputs/05_environmental_data/box_size_comparison_r.png / .pdf
#             heatmap of r by site x variable, one panel per period
# Author  : Natalie Posdaljian
# Date    : 2026-11-18
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Predictors are those of the model input table (05_compile_
#             dataset.R) at the species depths: temperature, salinity, o2,
#             EKE and EKE_mad at every nominal depth in ALL_DEPTHS_M (nearest
#             model level), chl-a and productivity at nominal depths <=
#             CHL_PROD_MAX_DEPTH_M, and the surface predictors SSH,
#             mixed_layer, ice_conc, ice_thickness, FSLE and fsle_orient.
#             Lagged chl-a/productivity are window means of chla_0 /
#             productivity_0 and are not compared separately.
#           - EKE and EKE_mad are derived as in 05_compile_dataset.R (EKE:
#             anomalies from the 10-year mean per site, level and box;
#             EKE_mad = 0.5 (u_mad^2 + v_mad^2) x 1e4). Nominal depths whose
#             nearest level has no water in a box use the deepest level with
#             water at that site, as in 05.
#           - fsle_orient is an axial angle: differences are wrapped to
#             [-90, 90) degrees and r is the circular correlation
#             (Jammalamadaka & SenGupta 2001) of the doubled angles.
#           - Raw daily series (not deseasoned): seasonal cycles shared by
#             both boxes inflate r over the climatology window, so the
#             recording-period values are the more conservative ones.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

library(tidyverse)

F_BOX_COMPARISON     <- out_path("05_environmental_data", "box_size_comparison.csv")
F_BOX_COMPARISON_FIG <- out_path("05_environmental_data", "box_size_comparison_r")
BOX_SUFFIX <- c(main = "", box40 = "_box40")
with_suffix <- function(file, suffix) sub("\\.csv$", paste0(suffix, ".csv"), file)

# -----------------------------------------------------------------------------
# 1. Build the predictor table for one box (same rules as 05_compile_dataset)
# -----------------------------------------------------------------------------
nearest_level <- function(levels, d) levels[which.min(abs(levels - d))]

build_predictors <- function(suffix) {
  phy  <- read_csv(with_suffix(F_COP_PHY_DEPTH, suffix), show_col_types = FALSE)
  bgc  <- read_csv(with_suffix(F_COP_BGC_DEPTH, suffix), show_col_types = FALSE)
  sfc  <- read_csv(with_suffix(F_COP_SURFACE, suffix),   show_col_types = FALSE)
  fsle <- read_csv(with_suffix(F_FSLE_DAILY, suffix),    show_col_types = FALSE)

  # Nearest level per nominal depth; where it has no water in this box at a
  # site, the deepest level with water at that site (as in 05).
  site_map <- function(d) {
    levels <- sort(unique(d$depth_m))
    expand_grid(site = SITE_ORDER, nominal = ALL_DEPTHS_M) %>%
      rowwise() %>%
      mutate(available = list(unique(d$depth_m[d$site == site])),
             depth_m = {
               lv <- nearest_level(levels, nominal)
               if (lv %in% available) lv else max(available)
             }) %>%
      ungroup() %>%
      select(site, nominal, depth_m)
  }

  phy_long <- phy %>%
    inner_join(site_map(phy), by = c("site", "depth_m"), relationship = "many-to-many") %>%
    group_by(site, nominal) %>%
    mutate(EKE     = 0.5 * ((u - mean(u))^2 + (v - mean(v))^2) * 1e4,
           EKE_mad = 0.5 * (u_mad^2 + v_mad^2) * 1e4) %>%
    ungroup() %>%
    select(date, site, nominal, temperature, salinity, EKE, EKE_mad) %>%
    pivot_longer(-c(date, site, nominal), names_to = "var")

  bgc_long <- bgc %>%
    inner_join(site_map(bgc), by = c("site", "depth_m"), relationship = "many-to-many") %>%
    select(date, site, nominal, chla, o2, productivity) %>%
    pivot_longer(-c(date, site, nominal), names_to = "var") %>%
    filter(var == "o2" | nominal <= CHL_PROD_MAX_DEPTH_M)

  surface_long <- sfc %>%
    transmute(date, site, SSH = ssh, mixed_layer, ice_conc, ice_thickness) %>%
    left_join(select(fsle, date, site, FSLE, fsle_orient), by = c("date", "site")) %>%
    pivot_longer(-c(date, site), names_to = "variable")

  bind_rows(phy_long, bgc_long) %>%
    transmute(date, site, variable = paste0(var, "_", nominal), value) %>%
    bind_rows(surface_long)
}

paired <- inner_join(build_predictors(BOX_SUFFIX[["main"]]),
                     build_predictors(BOX_SUFFIX[["box40"]]),
                     by = c("date", "site", "variable"), suffix = c("_main", "_box40")) %>%
  filter(!is.na(value_main), !is.na(value_box40))

# -----------------------------------------------------------------------------
# 2. Periods: recording days and the whole climatology window
# -----------------------------------------------------------------------------
rec_days <- SITES %>%
  transmute(site, date = map2(as.Date(rec_start), as.Date(rec_end), seq, by = "day")) %>%
  unnest(date) %>%
  mutate(period = "recording")
periods <- bind_rows(mutate(paired, period = "climatology"),
                     inner_join(paired, select(rec_days, site, date, period),
                                by = c("site", "date")))

# -----------------------------------------------------------------------------
# 3. Agreement metrics
# -----------------------------------------------------------------------------
# Axial (180-degree periodic) helpers for fsle_orient.
wrap_axial <- function(d) (d + 90) %% 180 - 90
circ_cor_axial <- function(a, b) {
  a2 <- 2 * a * pi / 180
  b2 <- 2 * b * pi / 180
  ma <- atan2(mean(sin(a2)), mean(cos(a2)))
  mb <- atan2(mean(sin(b2)), mean(cos(b2)))
  sum(sin(a2 - ma) * sin(b2 - mb)) / sqrt(sum(sin(a2 - ma)^2) * sum(sin(b2 - mb)^2))
}

metrics <- periods %>%
  group_by(site, variable, period) %>%
  summarise(
    n_days    = n(),
    r         = if (first(variable) == "fsle_orient") circ_cor_axial(value_main, value_box40)
                else suppressWarnings(cor(value_main, value_box40)),
    diff      = list(if (first(variable) == "fsle_orient") wrap_axial(value_main - value_box40)
                     else value_main - value_box40),
    sd_box40  = sd(value_box40),
    .groups = "drop") %>%
  mutate(mean_diff = map_dbl(diff, mean),
         rmse      = map_dbl(diff, ~ sqrt(mean(.x^2))),
         rmse_rel  = rmse / sd_box40) %>%
  select(site, variable, period, n_days, r, mean_diff, rmse, sd_box40, rmse_rel) %>%
  arrange(factor(site, levels = SITE_ORDER), variable, period)

write_csv(metrics, F_BOX_COMPARISON)
message("Wrote ", F_BOX_COMPARISON)

# Console summary (r over variables, per site and period). Variables that
# are constant in a period (e.g., no sea ice) give NA r and are counted.
print(metrics %>%
        group_by(period, site) %>%
        summarise(n_vars = n(), n_r_na = sum(is.na(r)),
                  r_min = min(r, na.rm = TRUE), r_median = median(r, na.rm = TRUE),
                  n_r_below_0.9 = sum(r < 0.9, na.rm = TRUE), .groups = "drop"))
message("Variables with r < 0.9 during recording periods:")
print(metrics %>% filter(period == "recording", r < 0.9) %>% select(site, variable, r, rmse_rel),
      n = Inf)

# -----------------------------------------------------------------------------
# 4. Figure: heatmap of r (site x variable), one panel per period
# -----------------------------------------------------------------------------
# Variables ordered by type then nominal depth, so the depth structure reads
# along the axis.
var_order <- metrics %>%
  distinct(variable) %>%
  mutate(base  = str_remove(variable, "_\\d+$"),
         depth = as.numeric(str_extract(variable, "(?<=_)\\d+$"))) %>%
  arrange(factor(base, levels = c("SSH", "mixed_layer", "ice_conc", "ice_thickness", "FSLE",
                                  "fsle_orient", "temperature", "salinity", "o2", "EKE",
                                  "EKE_mad", "chla", "productivity")), depth) %>%
  pull(variable)

fig <- metrics %>%
  mutate(variable = factor(variable, levels = rev(var_order)),
         site = factor(site, levels = SITE_ORDER),
         period = factor(period, levels = c("recording", "climatology"),
                         labels = c("Recording periods", "Climatology window"))) %>%
  ggplot(aes(site, variable, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.2) +
  scale_fill_viridis_c(name = "Pearson r\n(±20 vs ±40 km)", limits = c(min(0, metrics$r, na.rm = TRUE), 1),
                       na.value = "grey85") +
  facet_wrap(~ period) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank())

ggsave(paste0(F_BOX_COMPARISON_FIG, ".png"), fig, width = 6.5, height = 9, dpi = 300, bg = "white")
ggsave(paste0(F_BOX_COMPARISON_FIG, ".pdf"), fig, width = 6.5, height = 9)
message("Wrote ", F_BOX_COMPARISON_FIG, ".png/.pdf")
