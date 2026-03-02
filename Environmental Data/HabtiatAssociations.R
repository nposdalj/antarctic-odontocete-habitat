# =============================
# Response-first diagnostics for allData.csv (WIDE -> LONG)
# 1) Histograms: env distributions by presence vs absence
# 4) Lag heatmap: association strength across lags
# 6) Phase space: mean vs SD colored by presence
# =============================

library(tidyverse)
library(lubridate)
library(slider)

# -----------------------------
# USER SETTINGS (edit these)
# -----------------------------
csv_path <- "/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/allData.csv"
out_dir <- "F:/Antarctica"

date_col <- "date"
site_col <- "Site"

species_cols <- c("BW29","BW37","BW58","Gm","Pm","Oo")

# If species columns are minutes/day, keep threshold at 0
# If they're already 0/1, threshold still works fine.
presence_threshold <- 0

# Pick which environmental columns you want to use for histograms / phase space
# (These already exist in your file)
hist_vars <- c(
  "chla_stl", "npp_stl", "sst_stl", "sss_stl", "mlayer_stl",
  "chla_anomaly", "npp_anomaly", "sst_anomaly", "sss_anomaly", "mlayer_anomaly",
  "chla_sd_0", "productivity_sd_0", "mixed_layer_sd", "ssh_sd"
)

# LAG GROUPS (already computed in your file, 1–6 months)
lag_groups <- list(
  chla = paste0("chla_", 1:6, "mon"),
  npp  = paste0("productivity_", 1:6, "mon"),
  temp = paste0("temp_", 1:6, "mon"),
  sal  = paste0("salinity_", 1:6, "mon"),
  eke  = paste0("EKE_", 1:6, "mon")
)

# For phase space (mean vs variability):
# choose a "mean-like" x variable and an "SD-like" y variable
phase_x <- "chla_stl"     # or "chla_3mon", or "npp_stl", etc.
phase_y <- "chla_sd_0"    # or "productivity_sd_0", "mixed_layer_sd", etc.

# ============================
# LOAD
# ============================
df_raw <- readr::read_csv(csv_path, show_col_types = FALSE)

# Drop the annoying index column if it exists
df_raw <- df_raw %>% select(-any_of(c("...1")))

# Keep only variables that actually exist (so script never breaks)
hist_vars <- intersect(hist_vars, names(df_raw))
lag_groups <- lapply(lag_groups, intersect, names(df_raw))

stopifnot(all(c(date_col, site_col) %in% names(df_raw)))
stopifnot(all(species_cols %in% names(df_raw)))

df <- df_raw %>%
  mutate(
    .date = as.Date(.data[[date_col]]),
    .site = as.character(.data[[site_col]])
  )

# ============================
# WIDE -> LONG for species response
# ============================
dat_long <- df %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = ".species",
    values_to = ".response"
  ) %>%
  mutate(
    .response = as.numeric(.response),
    .present = .response > presence_threshold
  )

# ============================
# 1) HISTOGRAMS (presence vs absence)
# ============================
plot_hist_by_presence <- function(data, xcol) {
  ggplot(data, aes(x = .data[[xcol]], fill = .present)) +
    geom_histogram(position = "identity", bins = 40, alpha = 0.45, na.rm = TRUE) +
    facet_grid(.species ~ .site, scales = "free_y") +
    labs(
      title = paste0("Histogram: ", xcol, " (present vs absent)"),
      x = xcol, y = "Count", fill = "Present"
    ) +
    theme_bw() +
    theme(legend.position = "top")
}

# Make a list of histogram plots
hist_plots <- lapply(hist_vars, \(v) plot_hist_by_presence(dat_long, v))
names(hist_plots) <- hist_vars

# Print one to check
if (length(hist_plots) > 0) print(hist_plots[[1]])

# ============================
# 4) LAG HEATMAP (for one lag group at a time)
# Metric = Δ mean (z): present - absent, standardized within site/species/lag col
# ============================
lag_heatmap <- function(data, lag_cols, title_prefix = "") {
  lag_cols <- intersect(lag_cols, names(data))
  stopifnot(length(lag_cols) > 0)
  
  long <- data %>%
    select(.site, .species, .present, all_of(lag_cols)) %>%
    pivot_longer(cols = all_of(lag_cols), names_to = "lag_name", values_to = "x") %>%
    mutate(
      lag_months = as.numeric(stringr::str_extract(lag_name, "(?<=_)\\d+(?=mon$)"))
    ) %>%
    group_by(.site, .species, lag_name) %>%
    mutate(xz = as.numeric(scale(x))) %>%
    summarise(
      n = sum(!is.na(xz)),
      diff_means_z = mean(xz[.present], na.rm = TRUE) - mean(xz[!.present], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(lag_months = as.numeric(stringr::str_extract(lag_name, "\\d+")))
  
  ggplot(long, aes(x = lag_months, y = .species, fill = diff_means_z)) +
    geom_tile(na.rm = TRUE) +
    facet_wrap(~ .site, nrow = 1) +
    scale_x_continuous(breaks = sort(unique(long$lag_months))) +
    labs(
      title = paste0(title_prefix, " lag heatmap (Δ mean z: present - absent)"),
      x = "Lag (months)", y = "Species", fill = "Δ mean (z)"
    ) +
    theme_bw()
}

# Example: Chl-a lags heatmap
if (length(lag_groups$chla) > 0) {
  p_lag_chla <- lag_heatmap(dat_long, lag_groups$chla, "Chl-a")
  print(p_lag_chla)
}

# Example: NPP/productivity lags heatmap
if (length(lag_groups$npp) > 0) {
  p_lag_npp <- lag_heatmap(dat_long, lag_groups$npp, "NPP (productivity)")
  print(p_lag_npp)
}

# ============================
# 6) PHASE SPACE (mean vs variability)
# ============================
phase_space_plot <- function(data, xcol, ycol) {
  stopifnot(all(c(xcol, ycol) %in% names(data)))
  
  ggplot(data, aes(x = .data[[xcol]], y = .data[[ycol]], color = .present)) +
    geom_point(alpha = 0.35, size = 1, na.rm = TRUE) +
    facet_grid(.species ~ .site, scales = "free") +
    labs(
      title = paste0("Phase space: ", xcol, " vs ", ycol),
      x = xcol, y = ycol, color = "Present"
    ) +
    theme_bw() +
    theme(legend.position = "top")
}

if (all(c(phase_x, phase_y) %in% names(dat_long))) {
  p_phase <- phase_space_plot(dat_long, phase_x, phase_y)
  print(p_phase)
}

# -----------------------------
# Optional: Save plots
# -----------------------------
if (!is.null(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (nm in names(hist_plots)) {
    ggsave(filename = file.path(out_dir, paste0("hist_", nm, ".png")),
           plot = hist_plots[[nm]], width = 11, height = 7, dpi = 300)
  }
  
  if (exists("p_lag")) {
    ggsave(file.path(out_dir, paste0("lag_heatmap_", feat_for_heatmap, ".png")),
           p_lag, width = 12, height = 4, dpi = 300)
  }
  
  if (exists("p_phase")) {
    ggsave(file.path(out_dir, paste0("phase_", mean_col, "_vs_", sd_col, ".png")),
           p_phase, width = 11, height = 7, dpi = 300)
  }
}
