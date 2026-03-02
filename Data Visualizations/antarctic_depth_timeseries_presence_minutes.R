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

# ============================================================
# Antarctic Odontocete Habitat – DEPTH-SPECIFIC Env + Presence Plots
#
# What this does (per your request):
#  - For each Site × Species:
#     - choose DEPTHS based on species (your table)
#     - IGNORE surface depth (0.5)
#     - for EACH remaining depth, create ONE stacked figure that includes:
#         temperature, salinity, EKE, chla, productivity, o2,
#         temp_sd, salinity_sd, EKE_mad, chla_sd, productivity_sd, o2_sd
#       (each panel is a single LINE at that depth)
#     - presence at the bottom as MINUTES/DAY
#  - Save each figure separately with depth in the filename.
# ============================================================

# -------------------- USER SETTINGS --------------------
plot_dir <- "/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/plots/timeseries_depth_profiles"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

sites_to_plot   <- "ALL"  # c("EI","KGI","CI") or "ALL"
species_to_plot <- "ALL"  # c("Oo","Pm","Gm","BW29","BW37","BW58") or "ALL"

save_plots <- TRUE
dpi        <- 300
width_in   <- 10
height_in  <- 14

# -------------------- DEPTH-VARYING VARIABLES (ONE FIGURE) --------------------
depth_profile_vars <- c(
  # means
  "temperature", "salinity", "EKE", "chla", "productivity", "o2",
  # sd / mad
  "temp_sd", "salinity_sd", "EKE_mad", "chla_sd", "productivity_sd", "o2_sd"
)

depth_profile_title <- "Depth-varying variables (selected depth)"

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

# Your species→depth mapping (includes 0.5; we’ll drop it later)
species_depths <- function(sp) {
  if (sp == "BW29") return(c(0.5, 768.0))
  if (sp == "BW37" | sp == "BW58") return(c(0.5, 67.0, 920.0))
  if (sp == "Oo") return(c(0.5, 11.0, 455.0))
  if (sp == "Pm") return(c(0.5, 375.0, 1665.0))
  if (sp == "Gm") return(c(0.5, 16.0, 635.0))
  if (sp == "none") return(c(0.5))
  stop("Species code not valid: ", sp)
}

# Labels (add/adjust as needed)
var_label <- function(var) {
  if (var == "temperature") return("Temperature (°C)")
  if (var == "temp_sd") return("SD Temperature (°C)")
  if (var == "salinity") return("Salinity (psu)")
  if (var == "salinity_sd") return("SD Salinity (psu)")
  if (var == "EKE") return("Eddy Kinetic Energy")
  if (var == "EKE_mad") return("EKE Median Absolute Deviation")
  if (var == "chla") return("Chlorophyll (mg/m3)")
  if (var == "chla_sd") return("SD Chlorophyll (mg/m3)")
  if (var == "productivity") return("Net Primary Production")
  if (var == "productivity_sd") return("SD Net Primary Production")
  if (var == "o2") return("Oxygen (mmol/m3)")
  if (var == "o2_sd") return("SD Oxygen (mmol/m3)")
  var
}

# Convert presence to minutes/day (same auto-logic as before)
presence_to_minutes <- function(x) {
  x <- as.numeric(x)
  mx <- suppressWarnings(max(x, na.rm = TRUE))
  if (is.infinite(mx) || is.na(mx)) return(x)
  
  if (mx <= 1.05) {
    return(x * 1440)       # proportion -> minutes/day
  } else if (mx <= 288) {
    return(x * 5)          # 5-min bins -> minutes/day
  } else if (mx <= 1440) {
    return(x)              # already minutes/day
  } else {
    return(x)              # fallback
  }
}

# Keep ONE value per date after filtering to a single depth
# (avoids the “filled” look from repeated rows)
collapse_to_daily <- function(df) {
  df %>%
    arrange(date) %>%
    group_by(date, Site) %>%
    summarise(across(everything(), ~ dplyr::first(.x)), .groups = "drop")
}

# Nice filename depth string
depth_tag <- function(d) {
  # 768 -> "768m", 67 -> "67m", 1665 -> "1665m"
  paste0(gsub("\\.0+$", "", as.character(d)), "m")
}

# -------------------- PLOTTING --------------------
make_env_panel <- function(data, var, var_colors, show_x = FALSE) {
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

make_presence_panel <- function(data, species_code) {
  y_minutes <- presence_to_minutes(data[[species_code]])
  
  ggplot(data, aes(x = date, y = y_minutes)) +
    geom_col(width = 1, color = "darkmagenta", fill = "mediumorchid") +
    scale_x_date(date_labels = "%b %Y") +
    labs(y = NULL, x = NULL, title = paste0(species_code, " Presence (minutes/day)")) +
    theme(
      plot.margin = unit(c(0, 0.5, 0, 0.5), units = "line"),
      plot.title  = element_text(size = 10, margin = margin(t = 0, b = 0), face = "bold"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(color = "gray")
    )
}

aggregate_depth_profile_plot <- function(data_depth_daily, vars, species_code, site, depth_value, var_colors) {
  
  env_plots <- map(vars, ~{
    v <- .x
    show_x <- (v == tail(vars, 1))
    make_env_panel(data_depth_daily, v, var_colors, show_x = show_x)
  })
  
  pres <- make_presence_panel(data_depth_daily, species_code)
  
  main_title <- paste0(name(species_code), " at ", name(site), " — depth ", depth_tag(depth_value))
  subtitle   <- depth_profile_title
  
  wrap_plots(c(env_plots, list(pres)), ncol = 1) &
    plot_annotation(title = main_title, subtitle = subtitle) &
    theme(legend.position = "none")
}

save_depth_plot <- function(plot, site, species, depth_value, out_dir,
                            width_in = 10, height_in = 14, dpi = 300) {
  fname <- paste0(site, "_", species, "_depth_", depth_tag(depth_value), ".png")
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

# -------------------- MAIN DRIVER: loop sites × species × depth --------------------
timeseries_depth_profiles_all <- function(allData, sites, species, vars, out_dir,
                                          save_plots = TRUE, width_in = 10, height_in = 14, dpi = 300) {
  
  if (identical(sites, "ALL")) sites <- sort(unique(allData$Site))
  if (identical(species, "ALL")) {
    candidate_species <- c("BW29", "BW37", "BW58", "Gm", "Pm", "Oo")
    species <- candidate_species[candidate_species %in% names(allData)]
  }
  
  keep_cols <- unique(c("date", "Site", "depth", species, vars))
  df <- allData %>% select(any_of(keep_cols)) %>% filter(Site %in% sites)
  df$date <- as.Date(df$date)
  
  saved_files <- c()
  
  for (si in sites) {
    site_df <- df %>% filter(Site == si)
    if (nrow(site_df) == 0) next
    
    for (sp in species) {
      if (!sp %in% names(site_df)) next
      
      # Get species depths, drop surface (0.5) per your request
      depths <- setdiff(species_depths(sp), 0.5)
      
      # If a species only had surface, skip
      if (length(depths) == 0) next
      
      for (d in depths) {
        depth_df <- site_df %>% filter(depth == d)
        
        if (nrow(depth_df) == 0) {
          message("No data for ", si, " / ", sp, " at depth ", d)
          next
        }
        
        # Make it strictly daily (prevents “filled SD” look)
        depth_daily <- collapse_to_daily(depth_df)
        
        # Check missing vars
        missing_vars <- vars[!vars %in% names(depth_daily)]
        if (length(missing_vars) > 0) {
          message("Skipping ", si, " / ", sp, " / depth ", d,
                  " (missing: ", paste(missing_vars, collapse = ", "), ")")
          next
        }
        
        var_colors <- make_var_colors(vars)
        
        p <- aggregate_depth_profile_plot(
          data_depth_daily = depth_daily,
          vars = vars,
          species_code = sp,
          site = si,
          depth_value = d,
          var_colors = var_colors
        )
        
        print(p)
        
        if (save_plots) {
          f <- save_depth_plot(p, si, sp, d, out_dir,
                               width_in = width_in, height_in = height_in, dpi = dpi)
          saved_files <- c(saved_files, f)
        }
      }
    }
  }
  
  invisible(saved_files)
}

# -------------------- RUN --------------------
saved <- timeseries_depth_profiles_all(
  allData = allData,                  # assumes you already built this upstream (your merge steps)
  sites = sites_to_plot,
  species = species_to_plot,
  vars = depth_profile_vars,
  out_dir = plot_dir,
  save_plots = save_plots,
  width_in = width_in,
  height_in = height_in,
  dpi = dpi
)

if (save_plots && length(saved) > 0) {
  writeLines(saved, con = file.path(plot_dir, "saved_depth_profile_files.txt"))
  message("Saved ", length(saved), " depth-profile plot(s) to: ", plot_dir)
} else {
  message("No plots saved (either save_plots=FALSE or plots were skipped).")
}
