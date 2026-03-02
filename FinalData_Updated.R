# ============================================================
# Antarctic Odontocete Habitat – Environmental + Presence Plots
# FIXES:
#   1) Presence plotted as MINUTES per day (auto-detect scaling)
#   2) Plot SURFACE ONLY for all variables (depth == 0.5)
#   3) SD series are single LINE (no filled look)
#   4) De-seasoned titles indicate STL method
# ============================================================

library(tidyverse)
library(patchwork)
library(lubridate)

# -------------------- USER SETTINGS --------------------

# Output folder for plots
plot_dir <- "/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/plots/timeseries_3panel"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# What to plot
sites_to_plot   <- "ALL"   # c("EI","KGI","CI") or "ALL"
species_to_plot <- "ALL"   # c("Oo","Pm","Gm","BW29","BW37","BW58") or "ALL"

# Save settings
save_plots <- TRUE
dpi        <- 300
width_in   <- 10
height_in  <- 12

# -------------------- VARIABLE GROUPS (3 FIGURES) --------------------
vars_set1 <- c("FSLE","fsle_orient","fsle_sd","EKE","EKE_mad","SSH","ssh_sd","o2","o2_sd")
vars_set2 <- c("temperature","temp_sd","sst_anomaly","sst_stl","salinity","salinity_sd","sss_anomaly","sss_stl",
               "ice_conc","ice_thickness","ice_diff")
vars_set3 <- c("chla","chla_sd","chla_anomaly","chla_stl","productivity","productivity_sd","npp_anomaly","npp_stl",
               "mixed_layer")

var_sets <- list(set1 = vars_set1, set2 = vars_set2, set3 = vars_set3)
var_set_titles <- c(
  set1 = "FSLE / Dynamics / SSH / Oxygen",
  set2 = "Temperature / Salinity / Sea Ice",
  set3 = "Chlorophyll / NPP / Mixed Layer"
)

# -------------------- HELPERS --------------------

make_var_colors <- function(vars) {
  setNames(scales::hue_pal()(length(vars)), vars)
}

name <- function(abbrev) {
  if (abbrev == "CI") return("Clarence Island")
  if (abbrev == "KGI") return("King George Island")
  if (abbrev == "EI") return("Elephant Island")
  if (abbrev == "BW29") return("Southern Bottlenose Whale")
  if (abbrev == "BW37") return("Gray's and Strap-toothed Whale BW37")
  if (abbrev == "BW58") return("Gray's and Strap-toothed Whale BW58")
  if (abbrev == "Gm") return("Long-finned Pilot Whale")
  if (abbrev == "Oo") return("Killer Whale")
  if (abbrev == "Pm") return("Sperm Whale")
  return(abbrev)
}

# Labels (with STL callout for de-seasoned series)
var_label <- function(var) {
  lbl <- var
  if (var == "SSH") lbl <- "Sea Surface Height (m)"
  if (var == "ssh_sd") lbl <- "SD Sea Surface Height (m)"
  if (var == "FSLE") lbl <- "FSLE"
  if (var == "fsle_orient") lbl <- "Orientation of FSLE Vector"
  if (var == "fsle_sd") lbl <- "SD of FSLE"
  if (var == "EKE") lbl <- "Eddy Kinetic Energy"
  if (var == "EKE_mad") lbl <- "EKE Median Absolute Deviation"
  if (var == "o2") lbl <- "Oxygen (mmol/m3)"
  if (var == "o2_sd") lbl <- "SD Oxygen (mmol/m3)"

  if (var == "temperature") lbl <- "Temperature (°C)"
  if (var == "temp_sd") lbl <- "SD Temperature (°C)"
  if (var == "sst_anomaly") lbl <- "Temperature Anomaly (°C)"
  if (var == "sst_stl") lbl <- "De-seasoned Temperature (STL, °C)"

  if (var == "salinity") lbl <- "Salinity (psu)"
  if (var == "salinity_sd") lbl <- "SD Salinity (psu)"
  if (var == "sss_anomaly") lbl <- "Salinity Anomaly (psu)"
  if (var == "sss_stl") lbl <- "De-seasoned Salinity (STL, psu)"

  if (var == "ice_conc") lbl <- "Sea Ice Concentration"
  if (var == "ice_thickness") lbl <- "Sea Ice Thickness (m)"
  if (var == "ice_diff") lbl <- "Daily Change in Sea Ice Concentration"

  if (var == "chla") lbl <- "Chlorophyll (mg/m3)"
  if (var == "chla_sd") lbl <- "SD Chlorophyll (mg/m3)"
  if (var == "chla_anomaly") lbl <- "Chlorophyll Anomaly"
  if (var == "chla_stl") lbl <- "De-seasoned Chlorophyll (STL)"

  if (var == "productivity") lbl <- "Net Primary Production"
  if (var == "productivity_sd") lbl <- "SD Net Primary Production"
  if (var == "npp_anomaly") lbl <- "NPP Anomaly"
  if (var == "npp_stl") lbl <- "De-seasoned NPP (STL)"

  if (var == "mixed_layer") lbl <- "Mixed Layer Depth (m)"
  lbl
}

# ------------------------------------------------------------
# Presence conversion: convert whatever you have into minutes/day
# Auto rules:
#   - if max <= 1.05 : treat as proportion of day (0–1) -> *1440
#   - else if max <= 288 : treat as count of 5-min bins -> *5
#   - else if max <= 1440 : assume already minutes
#   - else: leave as-is (but you’ll notice)
# ------------------------------------------------------------
presence_to_minutes <- function(x) {
  x <- as.numeric(x)
  mx <- suppressWarnings(max(x, na.rm = TRUE))

  if (is.infinite(mx) || is.na(mx)) return(x)

  if (mx <= 1.05) {
    return(x * 1440)
  } else if (mx <= 288) {
    return(x * 5)
  } else if (mx <= 1440) {
    return(x)
  } else {
    return(x)
  }
}

# -------------------- IMPORTANT FILTER: SURFACE ONLY --------------------
# This prevents the “filled SD” look by ensuring only ONE value per date
surface_only <- function(df) {
  df %>%
    filter(depth == 0.5) %>%           # SURFACE ONLY
    group_by(date, Site) %>%
    summarise(across(everything(), ~ dplyr::first(.x)), .groups = "drop")
}

# ============================================================
# Antarctic Odontocete Habitat – Environmental + Presence Plots
# FIXES:
#   1) Presence plotted as MINUTES per day (auto-detect scaling)
#   2) Plot SURFACE ONLY for all variables (depth == 0.5)
#   3) SD series are single LINE (no filled look)
#   4) De-seasoned titles indicate STL method
# ============================================================

library(tidyverse)
library(patchwork)
library(lubridate)

# -------------------- USER SETTINGS --------------------

# Output folder for plots
plot_dir <- "/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/plots/timeseries_3panel"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# What to plot
sites_to_plot   <- "ALL"   # c("EI","KGI","CI") or "ALL"
species_to_plot <- "ALL"   # c("Oo","Pm","Gm","BW29","BW37","BW58") or "ALL"

# Save settings
save_plots <- TRUE
dpi        <- 300
width_in   <- 10
height_in  <- 12

# -------------------- VARIABLE GROUPS (3 FIGURES) --------------------
vars_set1 <- c("FSLE","fsle_orient","fsle_sd","EKE","EKE_mad","SSH","ssh_sd","o2","o2_sd")
vars_set2 <- c("temperature","temp_sd","sst_anomaly","sst_stl","salinity","salinity_sd","sss_anomaly","sss_stl",
               "ice_conc","ice_thickness","ice_diff")
vars_set3 <- c("chla","chla_sd","chla_anomaly","chla_stl","productivity","productivity_sd","npp_anomaly","npp_stl",
               "mixed_layer")

var_sets <- list(set1 = vars_set1, set2 = vars_set2, set3 = vars_set3)
var_set_titles <- c(
  set1 = "FSLE / Dynamics / SSH / Oxygen",
  set2 = "Temperature / Salinity / Sea Ice",
  set3 = "Chlorophyll / NPP / Mixed Layer"
)

# -------------------- HELPERS --------------------

make_var_colors <- function(vars) {
  setNames(scales::hue_pal()(length(vars)), vars)
}

name <- function(abbrev) {
  if (abbrev == "CI") return("Clarence Island")
  if (abbrev == "KGI") return("King George Island")
  if (abbrev == "EI") return("Elephant Island")
  if (abbrev == "BW29") return("Southern Bottlenose Whale")
  if (abbrev == "BW37") return("Gray's and Strap-toothed Whale BW37")
  if (abbrev == "BW58") return("Gray's and Strap-toothed Whale BW58")
  if (abbrev == "Gm") return("Long-finned Pilot Whale")
  if (abbrev == "Oo") return("Killer Whale")
  if (abbrev == "Pm") return("Sperm Whale")
  return(abbrev)
}

# Labels (with STL callout for de-seasoned series)
var_label <- function(var) {
  lbl <- var
  if (var == "SSH") lbl <- "Sea Surface Height (m)"
  if (var == "ssh_sd") lbl <- "SD Sea Surface Height (m)"
  if (var == "FSLE") lbl <- "FSLE"
  if (var == "fsle_orient") lbl <- "Orientation of FSLE Vector"
  if (var == "fsle_sd") lbl <- "SD of FSLE"
  if (var == "EKE") lbl <- "Eddy Kinetic Energy"
  if (var == "EKE_mad") lbl <- "EKE Median Absolute Deviation"
  if (var == "o2") lbl <- "Oxygen (mmol/m3)"
  if (var == "o2_sd") lbl <- "SD Oxygen (mmol/m3)"
  
  if (var == "temperature") lbl <- "Temperature (°C)"
  if (var == "temp_sd") lbl <- "SD Temperature (°C)"
  if (var == "sst_anomaly") lbl <- "Temperature Anomaly (°C)"
  if (var == "sst_stl") lbl <- "De-seasoned Temperature (STL, °C)"
  
  if (var == "salinity") lbl <- "Salinity (psu)"
  if (var == "salinity_sd") lbl <- "SD Salinity (psu)"
  if (var == "sss_anomaly") lbl <- "Salinity Anomaly (psu)"
  if (var == "sss_stl") lbl <- "De-seasoned Salinity (STL, psu)"
  
  if (var == "ice_conc") lbl <- "Sea Ice Concentration"
  if (var == "ice_thickness") lbl <- "Sea Ice Thickness (m)"
  if (var == "ice_diff") lbl <- "Daily Change in Sea Ice Concentration"
  
  if (var == "chla") lbl <- "Chlorophyll (mg/m3)"
  if (var == "chla_sd") lbl <- "SD Chlorophyll (mg/m3)"
  if (var == "chla_anomaly") lbl <- "Chlorophyll Anomaly"
  if (var == "chla_stl") lbl <- "De-seasoned Chlorophyll (STL)"
  
  if (var == "productivity") lbl <- "Net Primary Production"
  if (var == "productivity_sd") lbl <- "SD Net Primary Production"
  if (var == "npp_anomaly") lbl <- "NPP Anomaly"
  if (var == "npp_stl") lbl <- "De-seasoned NPP (STL)"
  
  if (var == "mixed_layer") lbl <- "Mixed Layer Depth (m)"
  lbl
}

# ------------------------------------------------------------
# Presence conversion: convert whatever you have into minutes/day
# Auto rules:
#   - if max <= 1.05 : treat as proportion of day (0–1) -> *1440
#   - else if max <= 288 : treat as count of 5-min bins -> *5
#   - else if max <= 1440 : assume already minutes
#   - else: leave as-is (but you’ll notice)
# ------------------------------------------------------------
presence_to_minutes <- function(x) {
  x <- as.numeric(x)
  mx <- suppressWarnings(max(x, na.rm = TRUE))
  
  if (is.infinite(mx) || is.na(mx)) return(x)
  
  if (mx <= 1.05) {
    return(x * 1440)
  } else if (mx <= 288) {
    return(x * 5)
  } else if (mx <= 1440) {
    return(x)
  } else {
    return(x)
  }
}

# -------------------- IMPORTANT FILTER: SURFACE ONLY --------------------
# This prevents the “filled SD” look by ensuring only ONE value per date
surface_only <- function(df) {
  df %>%
    filter(depth == 0.5) %>%           # SURFACE ONLY
    group_by(date, Site) %>%
    summarise(across(everything(), ~ dplyr::first(.x)), .groups = "drop")
}

# -------------------- PLOTTING FUNCTIONS --------------------

makePlot_surface <- function(data, var, var_colors, show_x = FALSE) {
  if (!var %in% names(data)) stop("Variable not found in data: ", var)
  
  label <- var_label(var)
  col   <- var_colors[[var]]
  
  x_scale <- if (show_x) scale_x_date(date_labels = "%b %Y") else scale_x_date(labels = NULL)
  
  ggplot(data, aes(x = date, y = .data[[var]])) +
    geom_line(color = col, linewidth = 1) +
    labs(y = NULL, x = NULL, title = label) +
    x_scale +
    theme(
      plot.margin = unit(c(0, 0.5, 0.3, 0.5), units = "line"),
      plot.title  = element_text(size = 10, margin = margin(t = 0, b = 0), face = "bold"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(color = "gray")
    )
}

presencePlot_minutes <- function(data, species_code, show_x = TRUE) {
  y_minutes <- presence_to_minutes(data[[species_code]])
  
  p <- ggplot(data, aes(x = date, y = y_minutes)) +
    geom_col(width = 1, color = "darkmagenta", fill = "mediumorchid") +
    labs(y = NULL, x = NULL, title = paste0(species_code, " Presence (minutes/day)")) +
    theme(
      plot.margin = unit(c(0, 0.5, 0, 0.5), units = "line"),
      plot.title  = element_text(size = 10, margin = margin(t = 0, b = 0), face = "bold"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(color = "gray")
    )
  
  if (show_x) {
    p <- p + scale_x_date(date_labels = "%b %Y")
  } else {
    p <- p + scale_x_date(labels = NULL)
  }
  
  p
}

aggregatePlot_surface <- function(data_surface, vars, species_code, site, fig_subtitle, var_colors) {
  
  env_plots <- map(vars, ~{
    v <- .x
    show_x <- (v == tail(vars, 1))
    makePlot_surface(data_surface, v, var_colors, show_x = show_x)
  })
  
  pres <- presencePlot_minutes(data_surface, species_code, show_x = TRUE)
  
  main_title <- paste0(name(species_code), " at ", name(site))
  
  wrap_plots(c(env_plots, list(pres)), ncol = 1, guides = "collect") &
    plot_annotation(title = main_title, subtitle = fig_subtitle) &
    theme(
      legend.position = "none"   # no depth legend since we’re surface-only
    )
}

save_plot_file <- function(plot, site, species, set_nm, out_dir,
                           width_in = 10, height_in = 12, dpi = 300) {
  fname <- paste0(site, "_", species, "_", set_nm, "_surface.png")
  fpath <- file.path(out_dir, fname)
  
  ggsave(
    filename = fpath,
    plot = plot,
    width = width_in,
    height = height_in,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  fpath
}

timeseriesPlots_3panel_surface_all <- function(allData, sites, species, var_sets, var_set_titles,
                                               out_dir, save_plots = TRUE,
                                               width_in = 10, height_in = 12, dpi = 300) {
  
  # auto sites/species
  if (identical(sites, "ALL")) sites <- sort(unique(allData$Site))
  if (identical(species, "ALL")) {
    candidate_species <- c("BW29", "BW37", "BW58", "Gm", "Pm", "Oo")
    species <- candidate_species[candidate_species %in% names(allData)]
  }
  
  # keep only needed cols
  all_vars <- unique(unlist(var_sets))
  keep_cols <- unique(c("date", "Site", "depth", species, all_vars))
  
  df <- allData %>% select(any_of(keep_cols)) %>% filter(Site %in% sites)
  df$date <- as.Date(df$date)
  
  saved_files <- c()
  
  for (si in sites) {
    # surface-only slice for site (key fix for “filled SD”)
    site_surface <- df %>%
      filter(Site == si) %>%
      surface_only()
    
    if (nrow(site_surface) == 0) next
    
    for (sp in species) {
      if (!sp %in% names(site_surface)) next
      
      for (set_nm in names(var_sets)) {
        vars <- var_sets[[set_nm]]
        
        # skip if missing vars
        missing_vars <- vars[!vars %in% names(site_surface)]
        if (length(missing_vars) > 0) {
          message("Skipping ", si, " / ", sp, " / ", set_nm,
                  " (missing: ", paste(missing_vars, collapse = ", "), ")")
          next
        }
        
        var_colors <- make_var_colors(vars)
        
        p <- aggregatePlot_surface(
          data_surface = site_surface,
          vars = vars,
          species_code = sp,
          site = si,
          fig_subtitle = var_set_titles[[set_nm]],
          var_colors = var_colors
        )
        
        print(p)
        
        if (save_plots) {
          f <- save_plot_file(p, si, sp, set_nm, out_dir,
                              width_in = width_in, height_in = height_in, dpi = dpi)
          saved_files <- c(saved_files, f)
        }
      }
    }
  }
  
  invisible(saved_files)
}

# -------------------- RUN IT --------------------
saved <- timeseriesPlots_3panel_surface_all(
  allData = allData,
  sites = sites_to_plot,
  species = species_to_plot,
  var_sets = var_sets,
  var_set_titles = var_set_titles,
  out_dir = plot_dir,
  save_plots = save_plots,
  width_in = width_in,
  height_in = height_in,
  dpi = dpi
)

if (save_plots && length(saved) > 0) {
  writeLines(saved, con = file.path(plot_dir, "saved_plot_files_surface.txt"))
  message("Saved ", length(saved), " plot(s) to: ", plot_dir)
} else {
  message("No plots saved (either save_plots=FALSE or plots were skipped due to missing variables).")
}
