# =============================================================================
# Script  : ice_edge_detections.R
# Section : 3.2 Sea ice - detections relative to the ice edge (descriptive)
# Purpose : For every recording day, the signed distance from the recorder to
#           the nearest AMSR2 15% ice edge (negative = recorder inside the
#           ice) and whether each species was detected. Shows how far inside
#           the ice edge each species was detected, and how often it was
#           detected on days when the recorder was inside the edge.
#           Descriptive only. The 6.25 km AMSR2 grid does not resolve leads
#           or small openings, so "inside the ice edge" does not mean the
#           whales were under unbroken ice.
# Inputs  : outputs/03_sea_ice/amsr2_daily_site.csv (amsr2_site_series.R)
#           F_PRESENCE_DAILY
# Outputs : outputs/08_results_gams/ice_edge_detections.png/.pdf
#           outputs/08_results_gams/ice_edge_detections_table.csv
# Author  : Natalie Posdaljian
# Date    : 2026-11-25
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : A day counts as a detection day if the species was detected in
#           any recorded minute. Days without an AMSR2 value are dropped.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages(library(tidyverse))

SP <- c(Pm = "Sperm whale", Gm = "Long-finned pilot whale", Oo = "Killer whale",
        BW29 = "Southern bottlenose whale (BW29)", BW37 = "BW37", BW58 = "BW58")
SECTION <- "08_results_gams"

ice <- read_csv(out_path("03_sea_ice", "amsr2_daily_site.csv"), show_col_types = FALSE) |>
  mutate(date = as.Date(date)) |> select(site, date, ice_conc, dist_ice_edge_km)
pres <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) |>
  mutate(date = as.Date(date)) |>
  filter(rec_minutes > 0) |>
  pivot_longer(all_of(names(SP)), names_to = "species", values_to = "minutes") |>
  mutate(detected = minutes > 0) |>
  inner_join(ice, by = c("site", "date")) |>
  filter(!is.na(dist_ice_edge_km)) |>
  mutate(site = factor(site, SITE_ORDER), species = factor(species, names(SP)))

# -----------------------------------------------------------------------------
# Table: per species and site
# -----------------------------------------------------------------------------
tab <- pres |>
  group_by(species, site) |>
  summarise(recording_days = n(),
            days_inside = sum(dist_ice_edge_km < 0),
            detection_days = sum(detected),
            detection_days_inside = sum(detected & dist_ice_edge_km < 0),
            pct_inside_days_with_detection = ifelse(days_inside > 0, 100 * detection_days_inside / days_inside, NA),
            pct_outside_days_with_detection = 100 * sum(detected & dist_ice_edge_km >= 0) / sum(dist_ice_edge_km >= 0),
            max_km_inside = ifelse(detection_days_inside > 0, -min(dist_ice_edge_km[detected]), NA),
            date_max_inside = if (any(detected & dist_ice_edge_km < 0))
              as.character(date[detected][which.min(dist_ice_edge_km[detected])]) else NA_character_,
            ice_pct_at_max_inside = ifelse(detection_days_inside > 0,
              100 * ice_conc[detected][which.min(dist_ice_edge_km[detected])], NA),
            .groups = "drop") |>
  mutate(across(where(is.double), \(x) round(x, 1)))
write_csv(tab, out_path(SECTION, "ice_edge_detections_table.csv"))

# -----------------------------------------------------------------------------
# Figure: all recording days (grey) and detection days (site colour) against
# signed distance to the ice edge, by species (rows) and site (columns)
# -----------------------------------------------------------------------------
lim <- range(pres$dist_ice_edge_km)
p <- ggplot(pres, aes(dist_ice_edge_km, species)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#DCE9F5") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(data = filter(pres, !detected), colour = "grey70", size = 0.5,
             position = position_jitter(height = 0.18, width = 0, seed = 1)) +
  geom_point(data = filter(pres, detected), aes(colour = site), size = 1.1,
             position = position_jitter(height = 0.18, width = 0, seed = 2)) +
  facet_wrap(~ site, ncol = 3) +
  scale_colour_manual(values = SITE_COLORS, guide = "none") +
  scale_y_discrete(limits = rev(names(SP)), labels = SP) +
  coord_cartesian(xlim = lim) +
  labs(x = "Distance to the 15% ice edge (km; negative = recorder inside the ice edge)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank())
for (ext in c("png", "pdf")) ggsave(out_path(SECTION, paste0("ice_edge_detections.", ext)), p,
                                    width = 8, height = 3.8, dpi = 300)
print(tab, n = Inf, width = Inf)
