# =============================================================================
# Script  : 08_eke_per_cell_comparison.R
# Section : 2.5 Environmental data - choice of EKE definition (supplement)
# Purpose : Compare the per-cell EKE used in the models (Ducet et al. 2000)
#           with the previous box-mean EKE in the mid-water and deep layers:
#           correlation and size ratio, and a figure of both versions
#           against daily whale presence over each recording period.
# Inputs  : outputs/05_environmental_data/eke_layers.csv (04b_eke_per_cell.R)
#           F_PRESENCE_DAILY
# Outputs : outputs/05_environmental_data/eke_definition_comparison.csv
#           outputs/05_environmental_data/eke_definition_comparison.png/.pdf
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages(library(tidyverse))

EKE_LAYERS <- c("mid", "deep")
both <- read_csv(out_path("05_environmental_data", "eke_layers.csv"), show_col_types = FALSE) |>
  filter(layer %in% EKE_LAYERS) |>
  mutate(date = as.Date(date)) |>
  rename(EKE_cell = EKE)

presence <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) |>
  filter(rec_minutes > 0) |> mutate(date = as.Date(date))
both_rec <- semi_join(both, presence, by = c("site", "date"))

comparison <- bind_rows(mutate(both, period = "10-year record"),
                        mutate(both_rec, period = "recording period")) |>
  group_by(period, site, layer) |>
  summarise(n_days = n(),
            median_box = median(EKE_box), median_cell = median(EKE_cell),
            ratio_cell_to_box = median(EKE_cell) / median(EKE_box),
            r_pearson = cor(EKE_box, EKE_cell),
            r_spearman = cor(EKE_box, EKE_cell, method = "spearman"), .groups = "drop") |>
  mutate(across(where(is.double), \(x) signif(x, 3)))
write_csv(comparison, out_path("05_environmental_data", "eke_definition_comparison.csv"))
print(comparison, n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# 3. Figure: both EKE versions (log scale) and daily presence per species,
#    over each site's recording period
# -----------------------------------------------------------------------------
eke_long <- both_rec |>
  pivot_longer(c(EKE_box, EKE_cell), names_to = "version", values_to = "EKE") |>
  mutate(EKE = log10(EKE)) |>
  mutate(version = recode(version, EKE_box = "box-mean (current)", EKE_cell = "per-cell (standard)"),
         row = factor(paste0("log10 EKE, ", ifelse(layer == "mid", "mid-water", "deep"), " layer"),
                      c("log10 EKE, mid-water layer", "log10 EKE, deep layer")))
pres_long <- presence |>
  pivot_longer(all_of(MODELED_SPECIES), names_to = "species", values_to = "minutes") |>
  mutate(row = factor("Presence (min per day)"), species = factor(species, MODELED_SPECIES))
site_lab <- \(d) mutate(d, site = factor(site, SITE_ORDER))

p <- ggplot() +
  geom_col(data = site_lab(pres_long), aes(date, minutes, fill = species), width = 1) +
  geom_line(data = site_lab(eke_long), aes(date, EKE, colour = version), linewidth = 0.35) +
  facet_grid(factor(row, c("Presence (min per day)", levels(eke_long$row))) ~ site,
             scales = "free", switch = "y") +
  scale_y_continuous() +
  scale_colour_manual(values = c("box-mean (current)" = "grey45", "per-cell (standard)" = "#E66100"),
                      name = "EKE (cm2 s-2)") +
  scale_fill_brewer(palette = "Dark2", name = "Species") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.placement = "outside", panel.grid.minor = element_blank())
for (ext in c("png", "pdf")) ggsave(out_path("05_environmental_data", paste0("eke_definition_comparison.", ext)),
                                    p, width = 9, height = 6.5, dpi = 300)
