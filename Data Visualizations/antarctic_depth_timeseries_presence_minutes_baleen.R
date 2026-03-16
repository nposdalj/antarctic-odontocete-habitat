library(tidyverse)
library(patchwork)
library(lubridate)

# -------------------- USER SETTINGS --------------------
plot_dir <- "/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/plots/timeseries_3panel" #CHANGE THIS TO WHERE YOU WANT TO SAVE YOUR PLOTS
input_file <- "L:/.shortcut-targets-by-id/16QptdYF6Cj57YnsWH6TVUR-l-QkP8voO/Baleen Whalesssss/mysticetes/allData_fin.csv" #CHANGE THIS TO YOUR DATA FOLDER
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

sites_to_plot   <- "ALL"   # or c("EI","KGI","CI")
species_to_plot <- "Bp"    # or "ALL" #CHANGE SPECIES
save_plots <- TRUE
dpi <- 300
width_in <- 10
height_in <- 12

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

# Variables that require depth-specific columns like temperature_54
depth_dependent_vars <- c(
  "temperature", "salinity", "EKE", "chla", "productivity", "o2",
  "temp_sd", "salinity_sd", "EKE_mad", "chla_sd", "productivity_sd", "o2_sd"
)

# -------------------- LOAD DATA --------------------
allData <- read.csv(input_file)
allData <- allData %>% subset(select=-X)
allData$date <- as.Date(allData$date, "%Y-%m-%d")

# -------------------- LOOKUPS --------------------
pretty_name <- c(
  CI   = "Clarence Island",
  KGI  = "King George Island",
  EI   = "Elephant Island",
  BW29 = "Southern Bottlenose Whale",
  BW37 = "Gray's and Strap-toothed Whale BW37",
  BW58 = "Gray's and Strap-toothed Whale BW58",
  Gm   = "Long-finned Pilot Whale",
  Oo   = "Killer Whale",
  Pm   = "Sperm Whale",
  Mn   = "Humpback",
  Bm   = "Blue Whale",
  Bp   ="Fin Whale"
)

#CHANGE DEPTHS HERE IF NEEDED
species_depth_lookup <- list(
  BW29 = 768,
  BW37 = 67,
  BW58 = 67,
  Oo   = 11,
  Pm   = 375,
  Gm   = 16,
  Mn   = 54,
  Bm   = 54,
  Bp   = 54
)

get_pretty_name <- function(x) {
  if (x %in% names(pretty_name)) pretty_name[[x]] else x
}

get_species_depth <- function(sp) {
  if (!sp %in% names(species_depth_lookup)) {
    stop("Species code not found in species_depth_lookup: ", sp)
  }
  species_depth_lookup[[sp]]
}

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

make_var_colors <- function(vars) {
  setNames(scales::hue_pal()(length(vars)), vars)
}

normalize_selection <- function(x, available_values) {
  if (length(x) == 1 && identical(x, "ALL")) {
    sort(unique(available_values))
  } else {
    intersect(x, available_values)
  }
}

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

depth_tag <- function(d) {
  paste0(gsub("\\.0+$", "", as.character(d)), "m")
}

# -------------------- DATA EXTRACTION --------------------
# For a requested variable:
# - if depth-dependent, pull var_depth
# - otherwise pull var directly
# Then rename everything back to the generic var name
build_plot_df <- function(df, vars, species_col, depth_value) {
  depth_suffix <- as.character(depth_value)
  
  source_cols <- map_chr(vars, function(v) {
    if (v %in% depth_dependent_vars) {
      paste0(v, "_", depth_suffix)
    } else {
      v
    }
  })
  
  keep <- tibble(var = vars, source = source_cols) %>%
    filter(source %in% names(df))
  
  if (nrow(keep) == 0) return(NULL)
  
  out <- df %>%
    select(any_of(c("date", "Site", species_col, keep$source))) %>%
    rename(!!!setNames(keep$source, keep$var)) %>%
    arrange(date) %>%
    group_by(date, Site) %>%
    summarise(across(everything(), ~ dplyr::first(.x)), .groups = "drop")
  
  out
}

# -------------------- PLOTTING --------------------
make_env_panel <- function(data, var, var_colors, show_x = FALSE) {
  if (!var %in% names(data)) stop("Variable not found in data: ", var)
  
  label <- var_label(var)
  col   <- var_colors[[var]]
  
  x_scale <- if (show_x) scale_x_date(date_labels = "%b %Y") else scale_x_date(labels = NULL)
  
  ggplot(data, aes(x = date, y = .data[[var]])) +
    geom_line(color = col, linewidth = 1, na.rm = TRUE) +
    labs(y = NULL, x = NULL, title = label) +
    x_scale +
    theme(
      plot.margin = unit(c(0, 0.5, 0.3, 0.5), units = "line"),
      plot.title  = element_text(size = 10, margin = margin(t = 0, b = 0), face = "bold"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(color = "gray"),
      panel.grid.minor = element_blank()
    )
}

make_presence_panel <- function(data, species_code) {
  y_minutes <- presence_to_minutes(data[[species_code]])
  
  ggplot(data, aes(x = date, y = y_minutes)) +
    geom_col(width = 1, color = "darkmagenta", fill = "mediumorchid", na.rm = TRUE) +
    scale_x_date(date_labels = "%b %Y") +
    labs(y = NULL, x = NULL, title = paste0(species_code, " Presence (minutes/day)")) +
    theme(
      plot.margin = unit(c(0, 0.5, 0, 0.5), units = "line"),
      plot.title  = element_text(size = 10, margin = margin(t = 0, b = 0), face = "bold"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(color = "gray"),
      panel.grid.minor = element_blank()
    )
}

make_set_plot <- function(data_plot, vars, species_code, site, set_title, depth_value) {
  vars_present <- vars[vars %in% names(data_plot)]
  vars_present <- vars_present[!map_lgl(data_plot[vars_present], ~ all(is.na(.x)))]
  
  if (length(vars_present) == 0) return(NULL)
  
  var_colors <- make_var_colors(vars_present)
  
  env_plots <- map(vars_present, ~{
    v <- .x
    show_x <- (v == tail(vars_present, 1))
    make_env_panel(data_plot, v, var_colors, show_x = show_x)
  })
  
  pres <- make_presence_panel(data_plot, species_code)
  
  wrap_plots(c(env_plots, list(pres)), ncol = 1) &
    plot_annotation(
      title = paste0(get_pretty_name(species_code), " at ", get_pretty_name(site), " — depth ", depth_tag(depth_value)),
      subtitle = set_title
    ) &
    theme(legend.position = "none")
}

save_set_plot <- function(plot, site, species, set_name, depth_value, out_dir,
                          width_in = 10, height_in = 12, dpi = 300) {
  fname <- paste0(site, "_", species, "_", set_name, "_depth_", depth_tag(depth_value), ".png")
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

# -------------------- MAIN DRIVER --------------------
timeseries_three_sets_all <- function(allData, sites, species, var_sets, var_set_titles, out_dir,
                                      save_plots = TRUE, width_in = 10, height_in = 12, dpi = 300) {
  
  available_sites <- unique(allData$Site)
  available_species <- intersect(names(species_depth_lookup), names(allData))
  
  sites <- normalize_selection(sites, available_sites)
  species <- normalize_selection(species, available_species)
  
  saved_files <- character(0)
  
  for (si in sites) {
    site_df <- allData %>% filter(Site == si)
    if (nrow(site_df) == 0) next
    
    for (sp in species) {
      if (!sp %in% names(site_df)) next
      
      depth_value <- get_species_depth(sp)
      
      for (set_name in names(var_sets)) {
        vars <- var_sets[[set_name]]
        
        plot_df <- build_plot_df(
          df = site_df,
          vars = vars,
          species_col = sp,
          depth_value = depth_value
        )
        
        if (is.null(plot_df)) {
          message("No usable columns for ", si, " / ", sp, " / ", set_name)
          next
        }
        
        p <- make_set_plot(
          data_plot = plot_df,
          vars = vars,
          species_code = sp,
          site = si,
          set_title = var_set_titles[[set_name]],
          depth_value = depth_value
        )
        
        if (is.null(p)) {
          message("No plottable data for ", si, " / ", sp, " / ", set_name)
          next
        }
        
        print(p)
        
        if (save_plots) {
          saved_files <- c(
            saved_files,
            save_set_plot(
              plot = p,
              site = si,
              species = sp,
              set_name = set_name,
              depth_value = depth_value,
              out_dir = out_dir,
              width_in = width_in,
              height_in = height_in,
              dpi = dpi
            )
          )
        }
      }
    }
  }
  
  invisible(saved_files)
}

# -------------------- RUN --------------------
saved <- timeseries_three_sets_all(
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
  writeLines(saved, con = file.path(plot_dir, "saved_three_set_files.txt"))
  message("Saved ", length(saved), " plot(s) to: ", plot_dir)
} else {
  message("No plots saved.")
}
