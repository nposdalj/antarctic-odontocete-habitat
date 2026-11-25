# =============================================================================
# Script  : 11_permutation_figure.R
# Section : Supplement - time-shift permutation test of the habitat models
# Purpose : Figure and table of the time-shift permutation test (mode C of
#           07_permutation_hypotheses.R): for each species x site, the null
#           distribution of deviance explained (selection repeated on
#           predictors shifted by whole years +/- up to 60 days) with the
#           observed value, and the permutation p-values.
# Inputs  : outputs/06_statistical_analyses/permutation/C_hypothesis_set_runs.csv
#           outputs/06_statistical_analyses/permutation/C_hypothesis_set_summary.csv
# Outputs : outputs/06_statistical_analyses/permutation/permutation_null.png/.pdf
# Author  : Natalie Posdaljian
# Date    : 2026-11-25
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Shifted runs in which no term was significant have deviance
#           explained 0 (no model selected).
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages(library(tidyverse))

perm_dir <- out_path("06_statistical_analyses", "permutation")
runs <- read_csv(file.path(perm_dir, "C_hypothesis_set_runs.csv"), show_col_types = FALSE) |>
  filter(!is.na(test))
summ <- read_csv(file.path(perm_dir, "C_hypothesis_set_summary.csv"), show_col_types = FALSE)

lab <- summ |>
  mutate(species = factor(species, names(SPECIES)), site = factor(site, SITE_ORDER)) |>
  arrange(species, site) |>
  mutate(panel = factor(paste(species, site), unique(paste(species, site))),
         txt = sprintf("obs. %.0f%%\np = %.2f", observed_dev_expl, p_perm))
null <- runs |> filter(run > 0) |>
  mutate(panel = factor(paste(species, site), levels(lab$panel)), site = factor(site, SITE_ORDER))

p <- ggplot(null, aes(dev_expl)) +
  geom_histogram(aes(fill = site), binwidth = 5, boundary = 0, alpha = 0.6) +
  geom_vline(data = lab, aes(xintercept = observed_dev_expl), linewidth = 0.6) +
  geom_text(data = lab, aes(x = 2, y = Inf, label = txt), hjust = 0, vjust = 1.2, size = 2.4) +
  facet_wrap(~ panel, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = SITE_COLORS, guide = "none") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Deviance explained by the selected model (%)", y = "Shifted runs") +
  theme_bw(base_size = 8) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank())
for (ext in c("png", "pdf")) ggsave(file.path(perm_dir, paste0("permutation_null.", ext)), p,
                                    width = 7.5, height = 5.5, dpi = 300)
print(select(lab, species, site, observed_terms, observed_dev_expl, null_median_dev_expl,
             null_95th_dev_expl, p_perm, p_bh), width = Inf)
