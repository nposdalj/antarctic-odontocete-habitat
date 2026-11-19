# =============================================================================
# Script  : sea_ice_occurrence.R
# Section : 3.2 Results (sea ice) - main-text figure and table
# Purpose : Describe acoustic occurrence across sea-ice conditions without a
#           model: for each species and site, the percentage of time bins
#           with detections and the mean daily presence minutes in each
#           sea-ice concentration class, with the number of bins per class.
#           Complements the sea-ice GAMs (supplement), whose inference is
#           limited by strong residual autocorrelation.
# Inputs  : <OUTPUT_ROOT>/06_statistical_analyses/stage1/binned_data/
#             <species>_binned.csv (all bins, before the ice filter; the same
#             bins as the sea-ice models)
#           config.R (SPECIES, SITE_ORDER, SITE_COLORS)
# Outputs : <OUTPUT_ROOT>/08_results_gams/
#             main_sea_ice_occurrence.png / .pdf
#             main_sea_ice_occurrence_table.csv
#               species, site, ice_class, n_bins, n_detection_bins,
#               pct_bins_detected, mean_presence_min_per_day
# Author  : Natalie Posdaljian
# Date    : 2026-11-19
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Ice classes ICE_BREAKS (fraction): 0-5%, 5-15%, 15-50%, >50%
#             (15% is the conventional sea-ice edge threshold).
#           - All species x site combinations with recordings are shown,
#             including those not modeled (fewer than 10 bins with
#             detections), since this figure is descriptive.
#           - Bin length per species as in the models (4 days; killer whale
#             2 days).
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

ICE_BREAKS <- c(-Inf, 0.05, 0.15, 0.50, Inf)
ICE_LABELS <- c("0–5", "5–15", "15–50", ">50")

binned <- map(MODELED_SPECIES, function(sp) {
  read_csv(file.path(STATS_DIR, "stage1", "binned_data", paste0(sp, "_binned.csv")),
           show_col_types = FALSE) |>
    transmute(species = sp, site, ice_conc, presence = .data[[sp]])
}) |>
  list_rbind()

occ <- binned |>
  mutate(ice_class = cut(ice_conc, ICE_BREAKS, labels = ICE_LABELS),
         species = factor(species, levels = MODELED_SPECIES),
         site = factor(site, levels = SITE_ORDER)) |>
  group_by(species, site, ice_class, .drop = FALSE) |>
  summarise(n_bins = n(), n_detection_bins = sum(presence > 0),
            pct_bins_detected = if (n() > 0) 100 * mean(presence > 0) else NA_real_,
            mean_presence_min_per_day = if (n() > 0) mean(presence) else NA_real_,
            .groups = "drop")
write_csv(occ, out_path("08_results_gams", "main_sea_ice_occurrence_table.csv"))

# Figure: rows = species, columns = sites; bars = % of bins with detections,
# labels = number of bins in the class (no bar where there were no bins).
plot_data <- occ |>
  mutate(label = if_else(n_bins > 0, paste0("n=", n_bins), "no bins"),
         y_label = coalesce(pct_bins_detected, 0))

fig <- ggplot(plot_data, aes(ice_class, pct_bins_detected, fill = site)) +
  geom_col(width = 0.75, na.rm = TRUE) +
  geom_text(aes(y = y_label, label = label), vjust = -0.4, size = 2.1, colour = "grey25") +
  facet_grid(species ~ site, labeller = labeller(species = \(x) species_label(x))) +
  scale_fill_manual(values = SITE_COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, 115), breaks = c(0, 50, 100),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Sea-ice concentration (%)", y = "Bins with detections (%)") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"))

ggsave(out_path("08_results_gams", "main_sea_ice_occurrence.png"), fig,
       width = 7.5, height = 6.5, dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(out_path("08_results_gams", "main_sea_ice_occurrence.pdf"), fig,
       width = 7.5, height = 6.5, device = cairo_pdf)
