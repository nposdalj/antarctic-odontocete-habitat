library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(sf)
library(units)

species <- 'Bp' # Bp, Bm, Mn

# ------------------------------

if (species == 'Bp') {
  data <- read.csv("D:/Mystecedes Time Series/allData_fin.csv")
}

if (species == 'Bm') {
  data <- read.csv("D:/Mystecedes Time Series/allData_blueZ (5).csv")
}

if (species == 'Mn') {
  data <- read.csv("D:/Mystecedes Time Series/allData_humpback (2).csv")
  data <- data %>%
    mutate(date = as.Date(parse_date_time(date, orders = c("mdy", "ymd", "dmy"))))
}

# -----------------------------

data_list <- readRDS("D:/CCAMLR Data/756/R/756_USA_2026-01-02.Rds")
C1        <- data_list[["C1"]]
C1_CATCH  <- data_list[["C1_CATCH"]]

PANEL_LEVELS <- c(species, "Krill Catch", "productivity_0", "chla_0",
                  "EKE_0", "o2_0", "salinity_0", "temperature_0")


ACOUSTIC_PANELS <- c(species, "productivity_0", "chla_0",
                     "EKE_0", "o2_0", "salinity_0", "temperature_0")
KRILL_PANELS    <- "Krill Catch"


daily_env_raw <- data %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2014-01-01")) %>%
  group_by(Site, date) %>%
  summarise(
    temperature_0    = mean(temperature_0,  na.rm = TRUE),
    salinity_0       = mean(salinity_0,     na.rm = TRUE),
    o2_0             = mean(o2_0,           na.rm = TRUE),
    EKE_0            = mean(EKE_0,          na.rm = TRUE),
    chla_0           = mean(chla_0,         na.rm = TRUE),
    productivity_0   = mean(productivity_0, na.rm = TRUE),
    species_sum      = sum(.data[[species]], na.rm = TRUE),
    recording_effort = n(),
    .groups = "drop"
  )

site_levels <- unique(daily_env_raw$Site)
date_range  <- seq(min(daily_env_raw$date), max(daily_env_raw$date), by = "day")

daily_env_full <- expand.grid(
  date = date_range,
  Site = site_levels,
  stringsAsFactors = FALSE
) %>%
  mutate(date = as.Date(date)) %>%
  left_join(daily_env_raw, by = c("date", "Site")) %>%
  mutate(recording_effort = replace_na(recording_effort, 0))


daily_krill <- C1_CATCH %>%
  filter(taxon_code == "KRI") %>%
  left_join(C1 %>% select(c1_id, date_catchperiod_start), by = "c1_id") %>%
  mutate(date = as.Date(date_catchperiod_start)) %>%
  filter(date >= as.Date("2014-01-01")) %>%
  group_by(date) %>%
  summarise(value = sum(greenweight_caught_kg, na.rm = TRUE), .groups = "drop") %>%
  mutate(variable = factor("Krill Catch", levels = PANEL_LEVELS),
         plot_type = "bar")


daily_env_long <- daily_env_raw %>%
  select(Site, date, temperature_0, salinity_0, o2_0, EKE_0, chla_0, productivity_0) %>%
  pivot_longer(
    cols      = c(temperature_0, salinity_0, o2_0, EKE_0, chla_0, productivity_0),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  mutate(plot_type = "line",
         variable  = factor(variable, levels = PANEL_LEVELS))

daily_whale <- daily_env_full %>%
  select(Site, date, species_sum) %>%
  rename(value = species_sum) %>%
  mutate(variable  = factor(species, levels = PANEL_LEVELS),
         plot_type = "bar")

daily_combined <- bind_rows(daily_env_long, daily_whale, daily_krill)


effort_acoustic <- daily_env_full %>%
  select(Site, date, recording_effort)

effort_krill <- C1_CATCH %>%
  filter(taxon_code == "KRI") %>%
  left_join(C1 %>% select(c1_id, date_catchperiod_start), by = "c1_id") %>%
  mutate(date = as.Date(date_catchperiod_start)) %>%
  filter(date >= as.Date("2014-01-01")) %>%
  group_by(date) %>%
  summarise(recording_effort = n(), .groups = "drop") %>%
  right_join(tibble(date = date_range), by = "date") %>%
  mutate(recording_effort = replace_na(recording_effort, 0),
         Site = "All")

rescale_effort <- function(effort_df, panel_ranges_df) {
  e_min <- min(effort_df$recording_effort)
  e_max <- max(effort_df$recording_effort)
  cross_join(panel_ranges_df, effort_df) %>%
    mutate(
      effort_rescaled = if (e_max > e_min) {
        (recording_effort - e_min) / (e_max - e_min) * (y_max - y_min) + y_min
      } else {
        y_min
      },
      variable = factor(variable, levels = PANEL_LEVELS)
    )
}

panel_ranges_env_species <- daily_combined %>%
  filter(variable %in% factor(ACOUSTIC_PANELS, levels = PANEL_LEVELS)) %>%
  group_by(variable) %>%
  summarise(y_min = min(value, na.rm = TRUE),
            y_max = max(value, na.rm = TRUE),
            .groups = "drop")

panel_ranges_krill <- daily_combined %>%
  filter(variable == "Krill Catch") %>%
  group_by(variable) %>%
  summarise(y_min = min(value, na.rm = TRUE),
            y_max = max(value, na.rm = TRUE),
            .groups = "drop")

effort_scaled_env_species <- rescale_effort(effort_acoustic, panel_ranges_env_species)
effort_scaled_krill       <- rescale_effort(effort_krill,    panel_ranges_krill)

effort_scaled <- bind_rows(effort_scaled_env_species, effort_scaled_krill)

# ------------------------------
# ZERO-EFFORT SHADING — panel-specific
# Acoustic panels: shade when ALL sites have zero acoustic effort
# Krill panel:     shade when there are zero hauls that day
# ------------------------------
n_sites <- length(site_levels)


acoustic_zero_runs <- effort_acoustic %>%
  group_by(date) %>%
  summarise(n_present = n(),
            n_zero    = sum(recording_effort == 0),
            .groups   = "drop") %>%
  filter(n_present == n_sites, n_zero == n_sites) %>%
  arrange(date) %>%
  mutate(run_id = cumsum(c(1, diff(as.numeric(date)) > 1))) %>%
  group_by(run_id) %>%
  summarise(xmin = min(date) - 0.5,
            xmax = max(date) + 0.5,
            .groups = "drop") %>%
  select(-run_id)


krill_zero_runs <- effort_krill %>%
  filter(recording_effort == 0) %>%
  arrange(date) %>%
  mutate(run_id = cumsum(c(1, diff(as.numeric(date)) > 1))) %>%
  group_by(run_id) %>%
  summarise(xmin = min(date) - 0.5,
            xmax = max(date) + 0.5,
            .groups = "drop") %>%
  select(-run_id)

message("Acoustic zero-effort periods: ", nrow(acoustic_zero_runs))
message("Krill zero-effort periods: ",    nrow(krill_zero_runs))


zero_effort_rects <- bind_rows(
  cross_join(
    acoustic_zero_runs,
    tibble(variable = factor(ACOUSTIC_PANELS, levels = PANEL_LEVELS))
  ),
  cross_join(
    krill_zero_runs,
    tibble(variable = factor(KRILL_PANELS, levels = PANEL_LEVELS))
  )
)


color_line <- c(
  temperature_0  = "#E69F00",
  salinity_0     = "#56B4E9",
  o2_0           = "#009E73",
  EKE_0          = "#CC79A7",
  chla_0         = "#2ca02c",
  productivity_0 = "#D55E00"
)

fill_bar <- c(setNames("#7851A9", species), "Krill Catch" = "lightpink")


ggplot(daily_combined, aes(x = date, y = value)) +
  
  geom_rect(
    data        = zero_effort_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey90", alpha = 0.6
  ) +
  
  geom_col(
    data     = subset(daily_combined, plot_type == "bar"),
    aes(fill = variable),
    position = "dodge", width = 1, alpha = 0.7
  ) +
  
  geom_line(
    data      = subset(daily_combined, plot_type == "line"),
    aes(color = variable, group = interaction(variable, Site)),
    linewidth = 0.7,
    na.rm     = TRUE
  ) +
  
  geom_line(
    data = effort_scaled,
    aes(x = date, y = effort_rescaled,
        group = interaction(variable, Site),
        linetype = Site),
    color = "grey40", linewidth = 0.4, alpha = 0.8
  ) +
  
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  scale_color_manual(values = color_line) +
  scale_fill_manual(values = fill_bar) +
  scale_linetype_manual(
    values = c("KGI" = "solid", "EI" = "dashed", "CI" = "dotted", "All" = "solid")
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y") +
  labs(
    title    = paste("Daily Environmental Conditions and", species, "Detections by Site"),
    subtitle = "Grey shading = zero effort  |  Grey lines = recording/fishing effort per site (rescaled)",
    x = "Date", y = "Value"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold"),
    plot.title      = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "grey50"),
    legend.position = "none"
  )

