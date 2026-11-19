# =============================================================================
# Script  : deep_oceanography.R
# Section : 2.2 Deep oceanography
# Purpose : Describe the full-water-column (0.5 to ~1684 m model levels)
#           temperature, salinity and dissolved oxygen at the HARP sites, as
#           context for the depths used by each species:
#             (a) one continuous depth x time section of temperature,
#                 salinity and oxygen, Jan 2014 - Dec 2016, stitched from the site that
#                 was (or was about to be) recording in each period, with
#                 species-relevant depths marked, each site's recording period
#                 outlined, and a sea-ice concentration strip on top;
#             (b) temperature-salinity (T-S) diagram, one panel per site
#                 (recording period only), points coloured by depth, with
#                 potential-density (sigma-theta) contours;
#             (c) table of mean +/- SD temperature, salinity and oxygen at each
#                 species depth (nearest model level) per site over its
#                 recording period.
# Inputs  : F_COP_PHY_DEPTH  (outputs/05_environmental_data/copernicus_phy_by_depth.csv)
#             date, site, depth_m, temperature, salinity, coverage, ...
#           F_PRESENCE_DAILY (outputs/04_acoustic_detection/presence_daily.csv)
#             species detected at each site (for the oxygen-depth notes)
#           F_COP_BGC_DEPTH  (outputs/05_environmental_data/copernicus_bgc_by_depth.csv)
#             date, site, depth_m, o2, coverage, ...
#           F_COP_SURFACE    (outputs/05_environmental_data/copernicus_surface.csv)
#             date, site, ice_conc, ...
#           F_COP_COVERAGE   (outputs/05_environmental_data/copernicus_depth_coverage.csv)
#             site, product, depth_m, n_cells, coverage
#           final_scripts_odontocetes/config.R (SITES, SPECIES, SITE_ORDER,
#             SITE_COLORS, BOX_HALF_KM, ...)
# Outputs : outputs/02_deep_oceanography/ts_sections.png / .pdf
#           outputs/02_deep_oceanography/ts_diagram.png / .pdf
#           outputs/02_deep_oceanography/ts_at_species_depths.csv
#           outputs/02_deep_oceanography/box_coverage_by_level.csv
# Author  : Natalie Posdaljian
# Date    : 2026-11-16 (revised 2026-11-18: main +/- BOX_HALF_KM box,
#           continuous 2014-2016 section, oxygen panel, captions moved to
#           the manuscript figure legends)
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - Values are daily means over each site's main analysis box
#             (+/- BOX_HALF_KM km). Temperature and salinity are from
#             GLOBAL_MULTIYEAR_PHY_001_030 (1/12 deg); oxygen is from
#             GLOBAL_MULTIYEAR_BGC_001_029 (1/4 deg, only a few cells per
#             box) and is missing below the deepest water in the BGC box.
#             Deep levels exist in only part of each box (the shelf is
#             shallower), so deep means represent only the cells deep enough;
#             the depth below which less than half of the box has water is
#             not drawn on the sections (it was easily confused with the
#             species-depth lines) but is saved in low_coverage_depth.csv
#             for the figure legend.
#           - The input files span each site's 10-year window, so the
#             continuous section can show the model at a site outside its
#             recording period; the actual recording periods are outlined.
#           - All days are shown (no sea-ice filter).
#           - Copernicus `thetao` is potential temperature and `so` practical
#             salinity, so sigma-theta is computed at zero pressure.
#           - No water-mass classification is applied.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggnewscale)
  library(oce)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
# Deepest level shown: the deepest species depth is 1665 m (sperm whale),
# whose nearest model level is ~1684 m.
max_depth_m <- 1700

# Box-coverage threshold used to report where deep means become
# representative of less than half of the box.
coverage_threshold <- 0.5

# Segments of the continuous section: which site's box supplies each period.
# Boundaries fall between deployments (KGI recording ends 2016-01-29, CI
# starts 2016-02-03), so every recording period lies inside its own segment.
section_segments <- tribble(
  ~site, ~seg_start,            ~seg_end,
  "EI",  as.Date("2014-01-01"), as.Date("2014-12-31"),
  "KGI", as.Date("2015-01-01"), as.Date("2016-01-31"),
  "CI",  as.Date("2016-02-01"), as.Date("2016-12-31")
) |>
  mutate(site = factor(site, levels = SITE_ORDER))

# Interpolation grid for the sections: regular in sqrt(depth) so the upper
# ocean (where most species depths are) gets more vertical room.
n_depth_grid <- 300

# Contour intervals drawn on the sections
temperature_contour_step <- 0.5   # deg C
salinity_contour_step    <- 0.1
oxygen_contour_step      <- 20    # mmol m-3

# Special characters built from code points so the script stays plain ASCII
degree_c           <- paste0(intToUtf8(176), "C")
section_label_temp <- paste0("Potential temperature (", degree_c, ")")
section_label_sal  <- "Salinity"
section_label_o2   <- expression("Oxygen (mmol m"^-3*")")   # plotmath: renders in PDF too

out_section <- "02_deep_oceanography"

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------
# Recording days: calendar days with at least MIN_DAY_CLOCK_MIN minutes of
# clock time inside the recording period, the same rule as in
# 04_acoustic_detection/build_presence.R (drops CI 2016-02-03, which has only
# 1 s of recording after the 23:59:59 start).
MIN_DAY_CLOCK_MIN <- 1
minutes_after_midnight <- function(t) {
  as.numeric(difftime(t, as.POSIXct(as.Date(t), tz = "UTC"), units = "mins"))
}

recording <- SITES |>
  transmute(
    site = factor(site, levels = SITE_ORDER),
    rec_start_date = as.Date(rec_start) +
      (24 * 60 - minutes_after_midnight(rec_start) < MIN_DAY_CLOCK_MIN),
    rec_end_date   = as.Date(rec_end) -
      (minutes_after_midnight(rec_end) < MIN_DAY_CLOCK_MIN)
  )

phy_all <- read_csv(F_COP_PHY_DEPTH, show_col_types = FALSE) |>
  select(date, site, depth_m, temperature, salinity) |>
  mutate(date = as.Date(date), site = factor(site, levels = SITE_ORDER)) |>
  filter(depth_m <= max_depth_m)

bgc_all <- read_csv(F_COP_BGC_DEPTH, show_col_types = FALSE) |>
  select(date, site, depth_m, oxygen = o2) |>
  mutate(date = as.Date(date), site = factor(site, levels = SITE_ORDER)) |>
  filter(depth_m <= max_depth_m)

ice_all <- read_csv(F_COP_SURFACE, show_col_types = FALSE) |>
  select(date, site, ice_conc) |>
  mutate(date = as.Date(date), site = factor(site, levels = SITE_ORDER))

# Main-box coverage of each model level (fraction of box cells with water),
# for the PHY grid (temperature, salinity) and the BGC grid (oxygen)
coverage_all <- read_csv(F_COP_COVERAGE, show_col_types = FALSE) |>
  filter(depth_m <= max_depth_m) |>
  transmute(site = factor(site, levels = SITE_ORDER), product, depth_m,
            coverage)
coverage_by_level <- coverage_all |>
  filter(product == "phy") |>
  select(-product)

# Each site's recording period, used by the T-S diagram and the table
clip_to_recording <- function(df) {
  df |>
    inner_join(recording, by = "site") |>
    filter(date >= rec_start_date, date <= rec_end_date) |>
    select(-rec_start_date, -rec_end_date)
}

# Each site's segment of the continuous 2014-2016 section
clip_to_segment <- function(df) {
  df |>
    inner_join(section_segments, by = "site") |>
    filter(date >= seg_start, date <= seg_end) |>
    select(-seg_start, -seg_end)
}

phy <- phy_all |>
  clip_to_recording() |>
  left_join(coverage_by_level, by = c("site", "depth_m"))

bgc <- clip_to_recording(bgc_all)

ice_section <- clip_to_segment(ice_all)

# Each variable is used only where it is defined (a level can be NA if no
# box cell is that deep).
phy_section_long <- bind_rows(
  phy_all |>
    clip_to_segment() |>
    pivot_longer(c(temperature, salinity), names_to = "variable",
                 values_to = "value"),
  bgc_all |>
    clip_to_segment() |>
    pivot_longer(oxygen, names_to = "variable", values_to = "value")
) |>
  filter(!is.na(value))

# Every day of the section must come from exactly one site
stopifnot(!anyDuplicated(distinct(phy_section_long, date, site)$date))

# -----------------------------------------------------------------------------
# 3. Species depths
# -----------------------------------------------------------------------------
# Nominal species depths from config (surface excluded from the section
# markers: it is the top of the plot). Species sharing a depth are combined
# into one label.
species_depths <- enframe(SPECIES, name = "species", value = "info") |>
  mutate(depth_m = map(info, "depths_m")) |>
  select(species, depth_m) |>
  unnest(depth_m) |>
  filter(species %in% MODELED_SPECIES)

depth_markers <- species_depths |>
  filter(depth_m > 0) |>
  group_by(depth_m) |>
  summarise(label = paste0(paste(species, collapse = "/"), " ", depth_m, " m"),
            .groups = "drop") |>
  arrange(depth_m)

# Axis labels for the markers: depths that would overprint each other on the
# square-root depth axis (e.g., 11 and 16 m) share one combined label placed
# at their mean depth; every depth still gets its own line.
min_label_gap_sqrt_m <- 1.5
marker_labels <- depth_markers |>
  mutate(label_group = cumsum(c(TRUE,
                                diff(sqrt(depth_m)) > min_label_gap_sqrt_m))) |>
  group_by(label_group) |>
  summarise(depth_m = mean(depth_m),
            label = paste(label, collapse = ", "),
            .groups = "drop")

# -----------------------------------------------------------------------------
# 4. Box coverage per site
# -----------------------------------------------------------------------------
# Shallowest model level at which less than `coverage_threshold` of the main
# box has water, per site and product (PHY for temperature and salinity, BGC
# for oxygen). Saved for the figure legend; the full site x level coverage
# table is saved too.
low_coverage_depth <- coverage_all |>
  filter(coverage < coverage_threshold) |>
  group_by(site, product) |>
  summarise(depth_m = min(depth_m), .groups = "drop")

write_csv(low_coverage_depth, out_path(out_section, "low_coverage_depth.csv"))

write_csv(
  coverage_all |>
    pivot_wider(names_from = site, values_from = coverage) |>
    rename_with(~ paste0("coverage_", .x), -c(product, depth_m)),
  out_path(out_section, "box_coverage_by_level.csv")
)

# -----------------------------------------------------------------------------
# 5. Interpolate profiles onto a regular depth grid
# -----------------------------------------------------------------------------
# Model levels are irregular (~1 m apart at the surface, ~230 m at depth).
# Each daily profile is linearly interpolated onto a grid regular in
# sqrt(depth), within the range of levels where the variable is defined (no
# extrapolation), so the section can be drawn as a raster.
sqrt_grid <- seq(sqrt(min(phy_section_long$depth_m)),
                 sqrt(max(phy_section_long$depth_m)),
                 length.out = n_depth_grid)
depth_grid <- sqrt_grid^2

interpolate_profile <- function(level_depth, level_value) {
  tibble(depth_m = depth_grid,
         value   = approx(level_depth, level_value, xout = depth_grid,
                          rule = 1)$y)
}

sections <- phy_section_long |>
  group_by(date, variable) |>
  filter(n() >= 2) |>
  summarise(profile = list(interpolate_profile(depth_m, value)),
            .groups = "drop") |>
  unnest(profile) |>
  filter(!is.na(value))

# -----------------------------------------------------------------------------
# 6. Figure (a): continuous depth x time sections with sea-ice strip
# -----------------------------------------------------------------------------
# Four stacked panels sharing the time axis: sea-ice concentration, then
# temperature, salinity and oxygen sections. Each section panel has its own
# colour bar, drawn beside the panel at the panel's full height (the bar's
# title runs vertically to its left). Depth axis: reversed square-root scale
# (surface at top). The secondary axis on the right labels the species depths
# drawn as dashed lines. Thin vertical lines mark where the source site
# changes; site-coloured outlines mark each site's recording period (no
# legend: the outlines match the site labels above the ice strip). Grey
# background: no water in the box at that depth (or below the deepest level
# of that product).
sqrt_reverse <- scales::transform_compose(scales::transform_sqrt(),
                                          scales::transform_reverse())
depth_breaks <- c(0, 50, 100, 200, 400, 800, 1200, 1600)

section_theme <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey85"),
        legend.position = "right",
        legend.key.height = unit(1, "null"),    # bar spans the panel height
        legend.key.width = unit(0.35, "cm"),
        legend.title.position = "left",
        legend.title = element_text(angle = 90, hjust = 0.5),
        legend.margin = margin(0, 0, 0, 0),
        axis.text.y.right = element_text(size = 7),
        axis.ticks.y.right = element_blank())

section_range <- c(min(section_segments$seg_start),
                   max(section_segments$seg_end) + 1)
date_scale <- scale_x_date(date_breaks = "3 months",
                           date_labels = "%b\n%Y", expand = c(0, 0))

# Site switches (start of every segment after the first)
site_switches <- section_segments$seg_start[-1]

# Recording-period outlines (same layer on every panel). Rectangles span the
# recording days inclusively.
recording_outline <- function(y_range = c(-Inf, Inf)) {
  list(
    geom_rect(data = recording,
              aes(xmin = rec_start_date, xmax = rec_end_date + 1,
                  ymin = y_range[1], ymax = y_range[2], color = site),
              fill = NA, linewidth = 0.7, inherit.aes = FALSE),
    scale_color_manual(values = SITE_COLORS, guide = "none")
  )
}

# One section panel.
section_panel <- function(var, label, colours, contour_step,
                          show_dates = FALSE) {
  section_data <- filter(sections, variable == var)

  panel <- ggplot(section_data, aes(x = date, y = depth_m)) +
    geom_raster(aes(fill = value)) +
    geom_contour(aes(z = value), binwidth = contour_step, color = "white",
                 linewidth = 0.15, alpha = 0.7) +
    scale_fill_gradientn(colours = colours, name = label) +
    # Species depths
    geom_hline(data = depth_markers, aes(yintercept = depth_m),
               linetype = "dashed", linewidth = 0.3, color = "black") +
    # Site switches and recording periods
    geom_vline(xintercept = site_switches, linewidth = 0.4,
               color = "grey15") +
    recording_outline(y_range = range(depth_grid)) +
    scale_y_continuous(transform = sqrt_reverse, breaks = depth_breaks,
                       expand = c(0, 0),
                       sec.axis = dup_axis(breaks = marker_labels$depth_m,
                                           labels = marker_labels$label,
                                           name = NULL)) +
    date_scale +
    coord_cartesian(xlim = section_range,
                    ylim = rev(range(depth_grid))) +
    labs(x = NULL, y = "Depth (m)") +
    section_theme

  if (!show_dates) panel <- panel + theme(axis.text.x = element_blank())
  panel
}

temperature_plot <- section_panel("temperature", section_label_temp,
                                  oceColorsTemperature(256),
                                  temperature_contour_step)
salinity_plot    <- section_panel("salinity", section_label_sal,
                                  oceColorsSalinity(256),
                                  salinity_contour_step)
oxygen_plot      <- section_panel("oxygen", section_label_o2,
                                  viridisLite::viridis(256),
                                  oxygen_contour_step, show_dates = TRUE)

# Oxygen-depth notes: where the 1/4-degree BGC box has no water below some
# level, species depths below it take oxygen from the deepest available level
# (as in the model predictors, 05_compile_dataset.R). Each affected segment is
# labelled in its grey (no-data) area. Gaps smaller than `o2_gap_min_m` (e.g.
# EI, deepest level 1653 m vs Pm 1665 m) are only the nearest-level rounding
# used everywhere else and are not labelled. Only species detected at a site
# are listed for it (e.g. no sperm whale detections at KGI, so no Pm model).
o2_gap_min_m <- 50
o2_deepest <- bgc_all |>
  filter(!is.na(oxygen)) |>
  group_by(site) |>
  summarise(deepest_m = max(depth_m), .groups = "drop")

detected_at_site <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) |>
  pivot_longer(any_of(names(SPECIES)), names_to = "species",
               values_to = "minutes") |>
  group_by(site, species) |>
  summarise(detected = any(minutes > 0), .groups = "drop") |>
  filter(detected) |>
  mutate(site = factor(site, levels = SITE_ORDER))

o2_notes <- species_depths |>
  filter(depth_m > 0) |>
  cross_join(o2_deepest) |>
  semi_join(detected_at_site, by = c("site", "species")) |>
  filter(depth_m > deepest_m + o2_gap_min_m) |>
  arrange(depth_m) |>
  group_by(site, deepest_m) |>
  summarise(species = paste0(species, " ", depth_m, " m", collapse = ", "),
            .groups = "drop") |>
  inner_join(section_segments, by = "site") |>
  mutate(
    mid = seg_start + (seg_end - seg_start) / 2,
    # vertical centre of the grey area, in square-root depth
    y = ((sqrt(deepest_m) + sqrt(max(depth_grid))) / 2)^2,
    label = paste0("No model oxygen below ", round(deepest_m), " m\n",
                   species, ": oxygen from ", round(deepest_m), " m")
  )

oxygen_plot <- oxygen_plot +
  geom_text(data = o2_notes, aes(x = mid, y = y, label = label),
            inherit.aes = FALSE, size = 3, color = "grey20",
            lineheight = 0.95)

# Sea-ice strip: daily box-mean ice concentration of the same source site,
# with the source-site labels above the panel.
site_labels <- section_segments |>
  mutate(mid = seg_start + (seg_end - seg_start) / 2)

ice_plot <- ggplot(ice_section, aes(x = date, y = ice_conc * 100)) +
  geom_area(aes(group = site), fill = "#9ecae1", color = "#3182bd",
            linewidth = 0.3) +
  geom_vline(xintercept = site_switches, linewidth = 0.4, color = "grey15") +
  recording_outline() +
  geom_text(data = site_labels, aes(x = mid, y = 118, label = site,
                                    color = site),
            inherit.aes = FALSE, fontface = "bold", size = 3.6) +
  scale_y_continuous(breaks = c(0, 50, 100),
                     expand = c(0, 0),
                     sec.axis = dup_axis(labels = NULL, name = NULL)) +
  date_scale +
  coord_cartesian(xlim = section_range, ylim = c(0, 100), clip = "off") +
  labs(x = NULL, y = "Ice (%)") +
  section_theme +
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        plot.margin = margin(26, 5.5, 2, 5.5))

# Each panel keeps its own legend (guides are not collected), so every colour
# bar sits next to its own section.
sections_figure <- wrap_plots(ice_plot, temperature_plot, salinity_plot,
                              oxygen_plot, ncol = 1,
                              heights = c(0.14, 1, 1, 1))

ggsave(out_path(out_section, "ts_sections.png"), sections_figure,
       width = 12, height = 9.5, dpi = 600, device = ragg::agg_png,
       bg = "white")
ggsave(out_path(out_section, "ts_sections.pdf"), sections_figure,
       width = 12, height = 9.5, device = cairo_pdf)

# -----------------------------------------------------------------------------
# 7. Figure (b): T-S diagram per site with sigma-theta contours
# -----------------------------------------------------------------------------
# Only levels where both temperature and salinity are defined are plotted.
ts_points <- phy |>
  filter(!is.na(temperature), !is.na(salinity)) |>
  arrange(depth_m)   # deep points drawn last so they are not hidden

# Potential-density anomaly on a T-S grid covering all sites (same axes in
# every panel). Inputs are potential temperature and practical salinity, so
# sigma-theta is evaluated at zero pressure (UNESCO EOS-80).
ts_pad <- 0.05
sigma_grid <- expand_grid(
  salinity    = seq(min(ts_points$salinity) - ts_pad,
                    max(ts_points$salinity) + ts_pad, length.out = 200),
  temperature = seq(min(ts_points$temperature) - ts_pad * 4,
                    max(ts_points$temperature) + ts_pad * 4, length.out = 200)
) |>
  mutate(sigma_theta = swSigmaTheta(salinity, temperature, pressure = 0,
                                    eos = "unesco"))

# Contour labels: every second isopycnal (0.2 kg m-3), placed where the
# isopycnal crosses a line just below the top of the panel. The crossing
# salinity is found by root-finding on the same equation of state.
sigma_contour_step <- 0.1
label_temperature  <- max(sigma_grid$temperature) - 0.25
salinity_range     <- range(sigma_grid$salinity)
sigma_at <- function(s) swSigmaTheta(s, label_temperature, pressure = 0,
                                     eos = "unesco")

sigma_labels <- tibble(
  sigma_theta = seq(ceiling(min(sigma_grid$sigma_theta) * 5) / 5,
                    max(sigma_grid$sigma_theta), by = 2 * sigma_contour_step)
) |>
  filter(sigma_theta > sigma_at(salinity_range[1]),
         sigma_theta < sigma_at(salinity_range[2])) |>
  mutate(
    salinity = map_dbl(sigma_theta, \(level) uniroot(
      \(s) sigma_at(s) - level, interval = salinity_range)$root),
    temperature = label_temperature,
    label = formatC(sigma_theta, format = "f", digits = 1)
  )

# Surface freezing temperature, for reference at the cold end of the diagram
freezing_line <- tibble(
  salinity = range(sigma_grid$salinity),
  temperature = swTFreeze(salinity, pressure = 0, eos = "unesco")
)

ts_plot <- ggplot(ts_points, aes(x = salinity, y = temperature)) +
  geom_contour(data = sigma_grid, aes(z = sigma_theta),
               binwidth = sigma_contour_step, color = "grey75",
               linewidth = 0.3) +
  geom_label(data = sigma_labels, aes(label = label), size = 2.3,
             color = "grey35", fill = "white", linewidth = 0,
             label.padding = unit(0.06, "lines")) +
  geom_line(data = freezing_line, linetype = "dashed", color = "grey40",
            linewidth = 0.3) +
  geom_point(aes(color = depth_m), size = 0.5, alpha = 0.5, shape = 16) +
  scale_color_viridis_c(option = "mako", direction = -1,
                        transform = "sqrt", breaks = depth_breaks[-1],
                        name = "Depth (m)") +
  facet_wrap(~ site, nrow = 1) +
  scale_x_continuous(breaks = seq(33, 35, by = 0.5)) +
  coord_cartesian(xlim = range(sigma_grid$salinity),
                  ylim = range(sigma_grid$temperature), expand = FALSE) +
  labs(x = "Salinity", y = section_label_temp) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.8, "lines"),
        legend.key.height = unit(1.1, "cm"))

ggsave(out_path(out_section, "ts_diagram.png"), ts_plot,
       width = 11, height = 4.2, dpi = 600, device = ragg::agg_png,
       bg = "white")
ggsave(out_path(out_section, "ts_diagram.pdf"), ts_plot,
       width = 11, height = 4.2, device = cairo_pdf)

# -----------------------------------------------------------------------------
# 8. Table: T, S and oxygen at species depths
# -----------------------------------------------------------------------------
# For each site and nominal species depth (including the surface), the
# nearest PHY level where temperature is defined (T, S) and the nearest BGC
# level where oxygen is defined are used, and daily values over the recording
# period are summarised as mean +/- SD. Where a species depth lies below the
# deepest oxygen level in the BGC box, oxygen comes from that deepest level
# (o2_model_depth_m shows which).
nearest_level <- function(target, levels) {
  levels[which.min(abs(levels - target))]
}

site_levels <- phy |>
  filter(!is.na(temperature), !is.na(salinity)) |>
  distinct(site, depth_m)

species_depth_table <- species_depths |>
  group_by(depth_m) |>
  summarise(species = paste(species, collapse = "/"), .groups = "drop") |>
  rename(nominal_depth_m = depth_m) |>
  cross_join(tibble(site = factor(SITE_ORDER, levels = SITE_ORDER))) |>
  rowwise() |>
  mutate(model_depth_m = nearest_level(
    nominal_depth_m, site_levels$depth_m[site_levels$site == site]
  )) |>
  ungroup()

ts_summary <- species_depth_table |>
  left_join(phy, by = c("site", "model_depth_m" = "depth_m")) |>
  group_by(site, species, nominal_depth_m, model_depth_m) |>
  summarise(n_days           = sum(!is.na(temperature)),
            box_coverage     = mean(coverage),
            temperature_mean = mean(temperature, na.rm = TRUE),
            temperature_sd   = sd(temperature, na.rm = TRUE),
            salinity_mean    = mean(salinity, na.rm = TRUE),
            salinity_sd      = sd(salinity, na.rm = TRUE),
            .groups = "drop") |>
  arrange(site, nominal_depth_m)

o2_levels <- bgc |>
  filter(!is.na(oxygen)) |>
  distinct(site, depth_m)

o2_summary <- species_depth_table |>
  select(site, nominal_depth_m) |>
  rowwise() |>
  mutate(o2_model_depth_m = nearest_level(
    nominal_depth_m, o2_levels$depth_m[o2_levels$site == site]
  )) |>
  ungroup() |>
  left_join(bgc, by = c("site", "o2_model_depth_m" = "depth_m"),
            relationship = "many-to-many") |>   # depths can share a level
  group_by(site, nominal_depth_m, o2_model_depth_m) |>
  summarise(o2_mean = mean(oxygen, na.rm = TRUE),
            o2_sd   = sd(oxygen, na.rm = TRUE),
            .groups = "drop")

ts_summary <- ts_summary |>
  left_join(o2_summary, by = c("site", "nominal_depth_m")) |>
  mutate(across(c(model_depth_m, o2_model_depth_m, box_coverage),
                ~ round(.x, 2)),
         across(starts_with(c("temperature", "salinity")), ~ round(.x, 3)),
         across(c(o2_mean, o2_sd), ~ round(.x, 1)))

write_csv(ts_summary, out_path(out_section, "ts_at_species_depths.csv"))
