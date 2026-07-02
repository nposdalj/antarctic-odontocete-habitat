library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)
library(purrr)
library(dplyr)
library(lubridate)

# =========================================================
# CONFIG
# =========================================================

plot_save_dir <- "L:/.shortcut-targets-by-id/16QptdYF6Cj57YnsWH6TVUR-l-QkP8voO/Baleen Whalesssss/mysticetes/GAMs/humpback"
allData <- read.csv("L:/.shortcut-targets-by-id/16QptdYF6Cj57YnsWH6TVUR-l-QkP8voO/Baleen Whalesssss/mysticetes/allData_humpback.csv")

# Choose species as a SINGLE string
species <- "Mn"   # options: "BW29", "BW37", "BW58", "Gm", "Pm", "Oo", "Bp", "Mn", "Bm"

# Sites to run
sites <- c("EI", "CI", "KGI")

# Define species-specific depths here ONLY
# 0-54 m: 111054
# 130-222 m: 130222
species_depths <- list(
  BW29 = c(0, 375, 1665),
  BW37 = c(0, 375, 1665),
  BW58 = c(0, 375, 1665),
  Gm   = c(0, 375, 1665),   # keep/edit as needed
  Pm   = c(0, 375, 1665),
  Oo   = c(0, 375, 1665),
  Bp   = c(111054, 130222),
  Mn   = c(111054, 130222),
  Bm   = c(111054, 130222)
)

# Pull depths automatically from selected species
depths <- species_depths[[species]]

if (is.null(depths)) {
  stop(paste("No depth configuration found for species:", species))
}

odontocete_species <- c("BW29", "BW37", "BW58", "Gm", "Pm", "Oo")
baleen_species    <- c("Bp", "Mn", "Bm")

include_ice_predictors <- species %in% baleen_species
apply_ice_filter <- species %in% odontocete_species

# -------------------------
# ACF block-model settings
# Define the variables to use for the ACF binning GAM at each site.
# Change ONLY here when you want different inputs.
# -------------------------
acf_terms <- list(
  EI = c("julian_day", "productivity_111054", "o2_130222"),
  CI = c("julian_day", "FSLE", "chla_111054", "o2_111054"),
  KGI = c("julian_day", "chla_111054", "FSLE", "temperature_111054", "o2_111054")
)

# -------------------------
# Predictor screening settings
# -------------------------
drop_depths <- numeric(0)
keep_lags <- c(3, 6)
chla_prod_depth_max <- 130222

themes_keywords <- list(
  mesoscale = c("FSLE", "SSH", "ssh", "EKE"),
  front_orientation = c("fsle_orient"),
  productivity = c("productivity", "chla"),
  oxygen = c("o2"),
  salinity = c("salinity"),
  stratification = c("mixed_layer"),
  temperature = c("temp", "temperature")
)

if (include_ice_predictors) {
  themes_keywords$ice <- c("ice_conc", "ice_thickness")
}

# Ice threshold quantile
ice_threshold_q <- 0.95

# =========================================================
# HELPERS
# =========================================================

name <- function(abbrev) {
  if (abbrev == "CI") {
    fullname <- "Clarence Island"
  } else if (abbrev == "KGI") {
    fullname <- "King George Island"
  } else if (abbrev == "EI") {
    fullname <- "Elephant Island"
  } else if (abbrev == "BW29") {
    fullname <- "Southern Bottlenose Whale"
  } else if (abbrev == "BW37") {
    fullname <- "Gray's and Strap-toothed Whale BW37"
  } else if (abbrev == "BW58") {
    fullname <- "Gray's and Strap-toothed Whale BW58"
  } else if (abbrev == "Gm") {
    fullname <- "Long-finned Pilot Whale"
  } else if (abbrev == "Oo") {
    fullname <- "Killer Whale"
  } else if (abbrev == "Pm") {
    fullname <- "Sperm Whale"
  } else if (abbrev == "Mn") {
    fullname <- "Humpback Whale"
  } else if (abbrev == "Bm") {
    fullname <- "Blue Whale"
  } else if (abbrev == "Bp") {
    fullname <- "Fin Whale"
  } else {
    fullname <- abbrev
  }
  return(fullname)
}

depth_regex <- function(depths) {
  paste0("_(", paste(depths, collapse = "|"), ")$")
}

build_vars_to_anom <- function(depths) {
  depth_bases <- c("temperature", "salinity", "o2", "chla", "productivity", "EKE")
  
  c(
    "SSH",
    "mixed_layer",
    "FSLE",
    "fsle_orient",
    as.vector(outer(depth_bases, depths, paste, sep = "_"))
  )
}

remove_zero_var <- function(df, vars) {
  vars[sapply(df[vars], function(x) sd(x, na.rm = TRUE) > 0)]
}

# =========================================================
# STEP 1: LOAD DATA
# =========================================================
# remove index-like columns if present
drop_cols <- intersect(c("X", "X.1"), names(allData))
if (length(drop_cols) > 0) {
  allData <- allData %>% dplyr::select(-dplyr::all_of(drop_cols))
}

allData$date <- as.Date(allData$date, "%Y-%m-%d")

all_species <- c("BW29", "BW37", "BW58", "Gm", "Pm", "Oo", "Ba", "Bm", "Bp")
drop_species <- setdiff(all_species, species)

keep_depth_regex <- depth_regex(depths)

sp_specific <- allData %>%
  dplyr::select(-dplyr::any_of(drop_species)) %>%
  dplyr::select(
    date,
    Site,
    julian_day,
    dplyr::all_of(species),
    dplyr::matches(keep_depth_regex),
    AAO:last_col()
  )

# =========================================================
# STEP 2: ACF BINNING
# =========================================================
precheck_save_dir <- file.path(plot_save_dir, "pre_model_checks")
dir.create(precheck_save_dir, recursive = TRUE, showWarnings = FALSE)

fit_acf_block_model <- function(data, site, species, predictors,
                                k = 4,
                                sp_val = 0.1,
                                lag.max = 1500) {
  
  site_data <- data %>%
    filter(Site == site)
  
  smooth_terms <- paste0("s(", predictors, ",k=", k, ",sp=", sp_val, ")")
  form <- as.formula(
    paste(species, "~", paste(smooth_terms, collapse = " + "))
  )
  
  mod <- gam(
    form,
    family = tw(link = "log", a = 1.1, b = 1.9),
    data = site_data,
    method = "REML"
  )
  
  ACF <- acf(residuals(mod), lag.max = lag.max, plot = FALSE)
  CI <- ggfortify:::confint.acf(ACF)
  ACFidx <- which(ACF[["acf"]] < CI, arr.ind = TRUE)
  acf_val <- ACFidx[1]
  
  list(
    model = mod,
    acf = ACF,
    acf_val = acf_val
  )
}

acf_results <- list()
acf_table <- tibble(site = character(), acf_val = numeric())

for (site in sites) {
  acf_results[[site]] <- fit_acf_block_model(
    data = allData,
    site = site,
    species = species,
    predictors = acf_terms[[site]]
  )
  
  acf_table <- bind_rows(
    acf_table,
    tibble(
      site = site,
      acf_val = acf_results[[site]]$acf_val
    )
  )
}

acfVal <- function(site) {
  acf_table %>%
    filter(site == !!site) %>%
    pull(acf_val)
}

binByACF <- function(data, site, bin_days, species, depths,
                     id_cols = c("date","Site"),
                     group_cols = c("bin_start","Site"),
                     exclude_cols = c("X","X.1")) {
  
  stopifnot(is.data.frame(data))
  stopifnot(length(species) == 1)
  
  sp_filtered <- data %>%
    filter(.data$Site == site) %>%
    mutate(bin_start = floor_date(.data$date, unit = paste(bin_days, "days")))
  
  mean_col <- function(col) rlang::expr(mean(.data[[!!col]], na.rm = TRUE))
  
  depth_regex <- paste0("_(", paste(depths, collapse = "|"), ")$")
  depth_cols <- names(sp_filtered)[grepl(depth_regex, names(sp_filtered))]
  
  species_cols <- intersect(species, names(sp_filtered))
  
  all_species <- c("BW29","BW37","BW58","Gm","Pm","Oo","Bp","Ba","Bm")
  other_species <- setdiff(all_species, species)
  
  non_depth_cols <- setdiff(
    names(sp_filtered),
    c(id_cols, other_species, exclude_cols, depth_cols)
  )
  
  non_depth_cols <- non_depth_cols[sapply(sp_filtered[non_depth_cols], is.numeric)]
  
  if ("julian_day" %in% names(sp_filtered) && !"julian_day" %in% non_depth_cols) {
    non_depth_cols <- c("julian_day", non_depth_cols)
  }
  
  species_expr <- set_names(
    lapply(species_cols, \(x) rlang::expr(mean(.data[[x]], na.rm = TRUE))),
    species_cols
  )
  
  non_depth_expr <- set_names(
    lapply(non_depth_cols, mean_col),
    non_depth_cols
  )
  
  depth_expr <- set_names(
    lapply(depth_cols, mean_col),
    depth_cols
  )
  
  all_summaries <- c(species_expr, non_depth_expr, depth_expr)
  
  sp_binned <- sp_filtered %>%
    group_by(across(all_of(intersect(group_cols, names(sp_filtered))))) %>%
    summarise(!!!all_summaries, .groups = "drop") %>%
    mutate(
      ice_regime = case_when(
        is.na(.data$ice_diff) ~ NA_character_,
        .data$ice_diff == 0 ~ "none",
        .data$ice_diff <= -0.01 ~ "decreasing",
        .data$ice_diff >=  0.01 ~ "increasing",
        TRUE ~ "stable"
      )
    )
  
  sp_binned
}

binned_data <- list()

for (site in sites) {
  binned_data[[site]] <- binByACF(
    data     = sp_specific,
    site     = site,
    bin_days = acfVal(site),
    species  = species,
    depths   = depths
  )
}

for (site in sites) {
  assign(
    x = paste0(site, "_binned"),
    value = binned_data[[site]],
    envir = .GlobalEnv
  )
}

save_acf_outputs <- function(acf_results, acf_table, species_name, save_dir) {
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Main site summary table
  write.csv(
    acf_table,
    file = file.path(save_dir, paste0("ACF_", species_name, "_site_summary.csv")),
    row.names = FALSE
  )
  
  # Per-site ACF values
  acf_values_all <- list()
  
  for (site in names(acf_results)) {
    current_acf <- acf_results[[site]]$acf
    
    acf_df_site <- tibble(
      site = site,
      lag = as.numeric(current_acf$lag),
      acf = as.numeric(current_acf$acf)
    )
    
    acf_values_all[[site]] <- acf_df_site
    
    write.csv(
      acf_df_site,
      file = file.path(save_dir, paste0("ACF_", species_name, "_", site, "_values.csv")),
      row.names = FALSE
    )
  }
  
  acf_values_all <- bind_rows(acf_values_all)
  
  write.csv(
    acf_values_all,
    file = file.path(save_dir, paste0("ACF_", species_name, "_all_sites_values.csv")),
    row.names = FALSE
  )
}

# =========================================================
# STEP 3: PLOT PRESENCE TIMESERIES
# =========================================================

binnedTimeseries <- function(data,site, bin) { # Function to create a timeseries plot
  # Making timeseries 
  ggplot(data = data, mapping = aes(x = bin_start, y = get(species))) + geom_col(width = 1, color = "slateblue") +
    scale_x_date(date_labels = "%b %Y")+
    labs(subtitle = name(site), y = NULL, x = NULL) + 
    theme(plot.subtitle = element_text(size = 9, face = "bold"), 
          plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), units = "line"))
  
}
binned_plot <- grid.arrange(binnedTimeseries(EI_binned,'EI',EI_acf), binnedTimeseries(CI_binned,'CI',CI_acf), nrow=3, 
                            top = paste('ACF Binned Species Presence for ', name(species), sep=''))

# =========================================================
# STEP 4: ICE FILTERING
# =========================================================

get_ice_thresholds <- function(df, response,
                               site_col = "Site",
                               q = 0.95) {
  df %>%
    filter(!is.na(ice_conc), !is.na(.data[[response]])) %>%
    mutate(present = .data[[response]] > 0) %>%
    group_by(.data[[site_col]]) %>%
    summarise(
      thr = quantile(ice_conc[present], probs = q, na.rm = TRUE),
      n_present = sum(present, na.rm = TRUE),
      .groups = "drop"
    )
}

save_ice_filter_outputs <- function(ice_thresholds,
                                    binned_data_before,
                                    binned_data_after,
                                    species_name,
                                    save_dir) {
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  ice_tbl_all <- list()
  
  for (site in names(ice_thresholds)) {
    thr_df <- ice_thresholds[[site]]
    if (is.null(thr_df) || nrow(thr_df) == 0) next
    
    n_before <- if (site %in% names(binned_data_before)) nrow(binned_data_before[[site]]) else NA_integer_
    n_after  <- if (site %in% names(binned_data_after))  nrow(binned_data_after[[site]])  else NA_integer_
    
    thr_df <- thr_df %>%
      mutate(
        species = species_name,
        n_before = n_before,
        n_after = n_after,
        n_removed = n_before - n_after
      ) %>%
      dplyr::relocate(species, .before = 1)
    
    ice_tbl_all[[site]] <- thr_df
  }
  
  if (length(ice_tbl_all) > 0) {
    ice_tbl_all <- bind_rows(ice_tbl_all)
    
    write.csv(
      ice_tbl_all,
      file = file.path(save_dir, paste0("IceFilter_", species_name, "_summary.csv")),
      row.names = FALSE
    )
    
    # Optional txt summary
    txt_lines <- c(
      paste0("Ice filtering summary for species: ", species_name),
      "",
      capture.output(print(ice_tbl_all))
    )
    
    writeLines(
      txt_lines,
      con = file.path(save_dir, paste0("IceFilter_", species_name, "_summary.txt"))
    )
  }
}

binned_data_pre_ice <- binned_data
ice_thresholds <- list()

if (apply_ice_filter) {
  
  for (site in sites) {
    if (!site %in% names(binned_data)) next
    
    current_df <- binned_data[[site]]
    
    ice_thresholds[[site]] <- get_ice_thresholds(
      df = current_df,
      response = species,
      q = ice_threshold_q
    )
    
    thr_val <- ice_thresholds[[site]]$thr[1]
    
    binned_data[[site]] <- current_df %>%
      filter(ice_conc < thr_val)
  }
  
  # re-create site-specific objects like EI_binned, CI_binned, etc.
  for (site in sites) {
    if (site %in% names(binned_data)) {
      assign(
        x = paste0(site, "_binned"),
        value = binned_data[[site]],
        envir = .GlobalEnv
      )
    }
  }
}

if (apply_ice_filter) {
  save_ice_filter_outputs(
    ice_thresholds = ice_thresholds,
    binned_data_before = binned_data_pre_ice,
    binned_data_after = binned_data,
    species_name = species,
    save_dir = precheck_save_dir
  )
}
# =========================================================
# STEP 5: DESEASON
# =========================================================

vars_to_anom <- build_vars_to_anom(depths)

# Mid-month day-of-year anchors (non-leap-year reference grid) used to place each
# calendar month's mean at a single representative day for spline interpolation.
mid_month_yday <- lubridate::yday(as.Date(paste0("2001-", 1:12, "-15")))

# Fits a periodic cubic spline through a site's 12 monthly means, anchored at
# mid-month, so the seasonal baseline is continuous across every month boundary
# (including Dec->Jan) instead of a flat step per calendar month.
fit_periodic_climatology <- function(month_vec, value_vec) {
  monthly_means <- tapply(value_vec, month_vec, mean, na.rm = TRUE)
  present_months <- as.integer(names(monthly_means))
  monthly_means <- as.numeric(monthly_means)

  keep <- is.finite(monthly_means)
  present_months <- present_months[keep]
  monthly_means <- monthly_means[keep]

  if (length(present_months) < 2) {
    return(list(fun = NULL, anchor_min = NA_real_))
  }

  ord <- order(present_months)
  present_months <- present_months[ord]
  monthly_means <- monthly_means[ord]

  x <- mid_month_yday[present_months]
  y <- monthly_means
  # Close the loop: repeat the earliest month's value exactly one year later, so
  # splinefun's periodic method (which requires y[1] == y[n]) wraps smoothly.
  x_period <- c(x, x[1] + 365)
  y_period <- c(y, y[1])

  list(
    fun = stats::splinefun(x_period, y_period, method = "periodic"),
    anchor_min = x[1]
  )
}

make_monthly_anoms_binned <- function(
    df,
    vars,
    date_col = "bin_start",
    site_col = "Site",
    standardize = FALSE
) {
  stopifnot(date_col %in% names(df), site_col %in% names(df))

  if (!inherits(df[[date_col]], c("Date", "POSIXt"))) {
    df[[date_col]] <- as.Date(df[[date_col]])
    if (any(is.na(df[[date_col]]))) {
      stop(sprintf("Column '%s' could not be coerced to Date.", date_col))
    }
  }

  vars <- intersect(vars, names(df))
  if (length(vars) == 0) stop("None of the requested vars are present in the data.")

  month_num <- lubridate::month(df[[date_col]])
  yday_num  <- lubridate::yday(df[[date_col]])

  for (v in vars) {
    anom_name <- paste0(v, "_anom")
    df[[anom_name]] <- NA_real_

    for (s in unique(df[[site_col]])) {
      site_rows <- which(df[[site_col]] == s)

      clim <- fit_periodic_climatology(month_num[site_rows], df[[v]][site_rows])
      if (is.null(clim$fun)) next

      yday_wrapped <- ifelse(
        yday_num[site_rows] < clim$anchor_min,
        yday_num[site_rows] + 365,
        yday_num[site_rows]
      )

      df[[anom_name]][site_rows] <- df[[v]][site_rows] - clim$fun(yday_wrapped)
    }
  }

  if (isTRUE(standardize)) {
    for (v in vars) {
      anom_name <- paste0(v, "_anom")
      z_name <- paste0(v, "_z")
      s <- stats::sd(df[[anom_name]], na.rm = TRUE)
      df[[z_name]] <- if (is.na(s) || s == 0) 0 else df[[anom_name]] / s
    }
  }

  df
}

binned_deseasoned <- list()

for (site in names(binned_data)) {
  binned_deseasoned[[site]] <- make_monthly_anoms_binned(
    binned_data[[site]],
    vars_to_anom,
    date_col = "bin_start",
    standardize = FALSE
  )
}

save_acf_outputs(
  acf_results = acf_results,
  acf_table = acf_table,
  species_name = species,
  save_dir = precheck_save_dir
)

# =========================================================
# STEP 6: STACKED ENVIRONMENTAL + PRESENCE TIMESERIES PLOTS
# (daily resolution and ACF-binned resolution, per site)
# =========================================================

stack_save_dir <- file.path(plot_save_dir, "stacked_timeseries")
dir.create(stack_save_dir, recursive = TRUE, showWarnings = FALSE)

# Daily monthly-mean-centered anomalies, computed the same way as binned_deseasoned
# (STEP 5), so the daily and binned stacked plots use consistently-computed "_anom" columns.
daily_deseasoned <- list()

for (site in sites) {
  site_daily <- sp_specific %>% filter(Site == site)
  
  daily_deseasoned[[site]] <- make_monthly_anoms_binned(
    site_daily,
    vars_to_anom,
    date_col = "date",
    site_col = "Site",
    standardize = FALSE
  )
}

# Variables to plot, in panel order (top to bottom)
stack_vars <- tibble::tribble(
  ~var,                      ~label,
  "EKE_mad_111054",          "EKE MAD (0-54m, cm²/s²)",
  "EKE_mad_130222",          "EKE MAD (130-222m, cm²/s²)",
  "FSLE",                    "FSLE Magnitude",
  "fsle_orient",             "FSLE Orientation",
  "mixed_layer_anom",        "De-seasoned Mixed Layer Depth (m)",
  "temperature_111054_anom", "De-seasoned Temperature (0-54m, °C)",
  "temperature_130222_anom", "De-seasoned Temperature (130-222m, °C)",
  "salinity_111054_anom",    "De-seasoned Salinity (0-54m, psu)",
  "salinity_130222_anom",    "De-seasoned Salinity (130-222m, psu)",
  "o2_111054_anom",          "De-seasoned Oxygen (0-54m, mmol/m3)",
  "o2_130222_anom",          "De-seasoned Oxygen (130-222m, mmol/m3)",
  "chla_111054",             "Chlorophyll (0-54m, mg/m3)",
  "ice_conc",                "Sea Ice Concentration"
)

stack_var_colors <- setNames(scales::hue_pal()(nrow(stack_vars)), stack_vars$var)

make_stack_env_panel <- function(data, var, label, col, date_col, show_x = FALSE) {
  x_scale <- if (show_x) scale_x_date(date_labels = "%b %Y") else scale_x_date(labels = NULL)
  
  ggplot(data, aes(x = .data[[date_col]], y = .data[[var]])) +
    geom_line(color = col, linewidth = 0.7, na.rm = TRUE) +
    labs(y = NULL, x = NULL, title = label) +
    x_scale +
    theme(
      plot.margin = unit(c(0, 0.5, 0.3, 0.5), units = "line"),
      plot.title  = element_text(size = 9, margin = margin(t = 0, b = 0), face = "bold"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(color = "gray"),
      panel.grid.minor = element_blank()
    )
}

make_stack_presence_panel <- function(data, species_col, date_col) {
  ggplot(data, aes(x = .data[[date_col]], y = .data[[species_col]])) +
    geom_col(width = 1, color = "darkmagenta", fill = "mediumorchid", na.rm = TRUE) +
    scale_x_date(date_labels = "%b %Y") +
    labs(y = NULL, x = NULL, title = paste0(species_col, " Presence")) +
    theme(
      plot.margin = unit(c(0, 0.5, 0, 0.5), units = "line"),
      plot.title  = element_text(size = 9, margin = margin(t = 0, b = 0), face = "bold"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(color = "gray"),
      panel.grid.minor = element_blank()
    )
}

build_stacked_plot <- function(data, site, species_col, date_col, title_suffix, show_acf_bin = TRUE) {
  vars_present <- stack_vars %>% filter(var %in% names(data))
  
  missing_vars <- setdiff(stack_vars$var, vars_present$var)
  if (length(missing_vars) > 0) {
    message(
      "Skipping missing columns for ", site, " (", title_suffix, "): ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  env_plots <- purrr::pmap(vars_present, function(var, label) {
    show_x <- (var == tail(vars_present$var, 1))
    make_stack_env_panel(data, var, label, stack_var_colors[[var]], date_col, show_x = show_x)
  })
  
  pres_plot <- make_stack_presence_panel(data, species_col, date_col)
  
  acf_suffix <- if (isTRUE(show_acf_bin)) paste0(", ACF Bin = ", acfVal(site), " days") else ""
  
  wrap_plots(c(env_plots, list(pres_plot)), ncol = 1) &
    plot_annotation(
      title = paste0(name(species_col), " at ", name(site), " (", title_suffix, acf_suffix, ")")
    ) &
    theme(legend.position = "none")
}

for (site in sites) {
  
  # --- Daily version ---
  if (site %in% names(daily_deseasoned)) {
    daily_plot <- build_stacked_plot(
      data = daily_deseasoned[[site]],
      site = site,
      species_col = species,
      date_col = "date",
      title_suffix = "Daily",
      show_acf_bin = FALSE
    )
    print(daily_plot)
    ggsave(
      filename = file.path(stack_save_dir, paste0("StackedTimeseries_", species, "_", site, "_daily.png")),
      plot = daily_plot,
      width = 10, height = 16, dpi = 300, bg = "white"
    )
  }
  
  # --- ACF-binned version ---
  if (site %in% names(binned_deseasoned)) {
    binned_plot <- build_stacked_plot(
      data = binned_deseasoned[[site]],
      site = site,
      species_col = species,
      date_col = "bin_start",
      title_suffix = "ACF-Binned"
    )
    print(binned_plot)
    ggsave(
      filename = file.path(stack_save_dir, paste0("StackedTimeseries_", species, "_", site, "_binned.png")),
      plot = binned_plot,
      width = 10, height = 16, dpi = 300, bg = "white"
    )
  }
}

# =========================================================
# STEP 7: VIF SCREENING
# =========================================================

vif_stepwise_drop <- function(data, response, predictors,
                              vif_thresh = 10,
                              family = gaussian(),
                              max_iter = 200,
                              verbose = TRUE) {
  preds <- intersect(predictors, names(data))
  if (length(preds) < 2) stop("Need at least 2 predictors present in `data`.")
  
  get_vif <- function(fit) {
    v <- tryCatch(car::vif(fit), error = function(e) NA)
    if (all(is.na(v))) return(v)
    
    if (is.matrix(v)) {
      if ("GVIF^(1/(2*Df))" %in% colnames(v)) {
        v_num <- v[, "GVIF^(1/(2*Df))"]
      } else if ("GVIF" %in% colnames(v)) {
        v_num <- v[, "GVIF"]
      } else {
        v_num <- v[, 1]
      }
      return(setNames(as.numeric(v_num), rownames(v)))
    }
    return(v)
  }
  
  fallback_drop <- function(df, preds) {
    X <- df[, preds, drop = FALSE]
    X <- X[, sapply(X, is.numeric), drop = FALSE]
    if (ncol(X) < 2) return(preds[1])
    
    cm <- suppressWarnings(stats::cor(X, use = "pairwise.complete.obs"))
    cm[lower.tri(cm, diag = TRUE)] <- NA
    max_pair <- which(abs(cm) == max(abs(cm), na.rm = TRUE), arr.ind = TRUE)
    if (nrow(max_pair) == 0) return(colnames(X)[1])
    
    i <- max_pair[1, 1]
    j <- max_pair[1, 2]
    a <- colnames(X)[i]
    b <- colnames(X)[j]
    
    mean_abs_cor <- function(var) {
      mean(abs(stats::cor(X[[var]], X, use = "pairwise.complete.obs")), na.rm = TRUE)
    }
    
    if (mean_abs_cor(a) >= mean_abs_cor(b)) a else b
  }
  
  log <- list()
  iter <- 0
  
  while (iter < max_iter && length(preds) >= 2) {
    iter <- iter + 1
    
    fml <- stats::as.formula(paste(response, "~", paste(preds, collapse = " + ")))
    fit <- tryCatch(stats::glm(fml, data = data, family = family),
                    error = function(e) NULL)
    
    if (is.null(fit)) {
      drop_var <- fallback_drop(data, preds)
      log[[iter]] <- list(step = iter, action = "drop (model failed)", dropped = drop_var,
                          max_vif = NA_real_, vif = NA)
      preds <- setdiff(preds, drop_var)
      if (verbose) message("[", iter, "] model failed; dropped: ", drop_var)
      next
    }
    
    v <- get_vif(fit)
    
    if (all(is.na(v)) || any(is.na(v))) {
      drop_var <- fallback_drop(data, preds)
      log[[iter]] <- list(step = iter, action = "drop (vif NA)", dropped = drop_var,
                          max_vif = NA_real_, vif = v)
      preds <- setdiff(preds, drop_var)
      if (verbose) message("[", iter, "] VIF NA; dropped: ", drop_var)
      next
    }
    
    max_v <- max(v, na.rm = TRUE)
    max_var <- names(v)[which.max(v)]
    
    log[[iter]] <- list(step = iter, action = "evaluate", dropped = NA_character_,
                        max_vif = max_v, worst = max_var, vif = v)
    
    if (verbose) message("[", iter, "] max VIF = ", round(max_v, 3), " (", max_var, ")")
    
    if (max_v <= vif_thresh) break
    
    preds <- setdiff(preds, max_var)
    log[[iter]]$action <- "drop (max vif)"
    log[[iter]]$dropped <- max_var
  }
  
  final_fml <- stats::as.formula(paste(response, "~", paste(preds, collapse = " + ")))
  final_fit <- tryCatch(stats::glm(final_fml, data = data, family = family),
                        error = function(e) NULL)
  
  list(
    kept_predictors = preds,
    dropped_predictors = setdiff(predictors, preds),
    final_formula = final_fml,
    final_fit = final_fit,
    steps = log
  )
}

build_preds <- function(df,
                        species,
                        depths,
                        drop_depths = numeric(0),
                        keep_lags = c(3, 6),
                        chla_prod_depth_max = 130222,
                        core_allow = c("FSLE", "SSH", "mixed_layer", "fsle_orient"),
                        keep_anom_core = c("FSLE", "SSH", "mixed_layer"),
                        keep_anom_depth_bases = c("temperature", "salinity", "o2"),
                        keep_EKE_mad = TRUE,
                        include_ice_predictors = FALSE) {
  
  nms <- names(df)
  
  excluded_vars <- c(
    "bin_start", "Site", species,
    "julian_day", "AAO",
    "ice_diff", "ice_regime",
    "bathymetry"
  )
  
  # Only exclude ice_conc and ice_thickness when NOT including them
  if (!isTRUE(include_ice_predictors)) {
    excluded_vars <- c(excluded_vars, "ice_conc", "ice_thickness")
  }
  
  pred <- setdiff(nms, excluded_vars)
  
  pred <- pred[!grepl("(_anomaly$|_stl$)", pred)]
  pred <- pred[!grepl("_sd", pred, ignore.case = TRUE)]
  
  if (isTRUE(keep_EKE_mad)) {
    drop_mad <- pred[grepl("_mad", pred, ignore.case = TRUE) & !grepl("^EKE_mad_\\d+$", pred)]
    pred <- setdiff(pred, drop_mad)
  } else {
    pred <- pred[!grepl("_mad", pred, ignore.case = TRUE)]
  }
  
  lag_vars_all <- pred[grepl("_\\d+mon$", pred)]
  lag_vars_all <- lag_vars_all[grepl("^(productivity|chla)_", lag_vars_all)]
  lag_nums <- as.integer(sub(".*_(\\d+)mon$", "\\1", lag_vars_all))
  lag_keep_vars <- lag_vars_all[lag_nums %in% keep_lags]
  
  depth_vars_all <- pred[grepl("_\\d+$", pred)]
  depth_vars_all <- setdiff(depth_vars_all, pred[grepl("_anom$", pred)])
  
  depth_nums <- as.integer(sub(".*_(\\d+)$", "\\1", depth_vars_all))
  
  keep_depth <- !(depth_nums %in% drop_depths)
  
  is_chla_prod <- grepl("^(chla|productivity)_", depth_vars_all)
  keep_chla_prod <- !(is_chla_prod & depth_nums > chla_prod_depth_max)
  
  depth_keep_vars <- depth_vars_all[keep_depth & keep_chla_prod]
  
  anom_vars_all <- pred[grepl("_anom$", pred)]
  core_vars <- setdiff(pred, c(depth_vars_all, lag_vars_all, anom_vars_all))
  core_keep <- intersect(core_vars, core_allow)
  
  # optionally keep raw ice vars
  ice_keep <- character(0)
  if (isTRUE(include_ice_predictors)) {
    ice_keep <- intersect(c("ice_conc", "ice_thickness"), pred)
  }
  
  anom_keep <- character(0)
  anom_keep <- c(anom_keep, paste0(keep_anom_core, "_anom"))
  
  depth_anom_candidates <- as.vector(outer(keep_anom_depth_bases, depths, paste, sep = "_"))
  anom_keep <- c(anom_keep, paste0(depth_anom_candidates, "_anom"))
  anom_keep <- intersect(anom_keep, pred)
  
  final <- sort(unique(c(core_keep, depth_keep_vars, lag_keep_vars, anom_keep, ice_keep)))
  final
}

assign_themes <- function(kept_predictors, themes_keywords) {
  theme_list <- list()
  
  for (theme in names(themes_keywords)) {
    keywords <- themes_keywords[[theme]]
    
    matched <- kept_predictors[
      sapply(kept_predictors, function(var) {
        any(sapply(keywords, function(k)
          grepl(paste0("^", k), var, ignore.case = TRUE)))
      })
    ]
    
    if (theme == "mesoscale") {
      matched <- matched[!grepl("^fsle_orient", matched, ignore.case = TRUE)]
    }
    
    if (length(matched) > 0) {
      theme_list[[theme]] <- matched
    }
  }
  
  return(theme_list)
}

select_best_by_theme <- function(data,
                                 response,
                                 themes,
                                 family = tw(link = "log", a = 1.1, b = 1.9),
                                 prefer_anoms = TRUE,
                                 anom_suffix = "_anom$",
                                 date_col = "bin_start",
                                 site_col = "Site",
                                 plot = TRUE,
                                 plot_dir = NULL,
                                 plot_prefix = "theme_compare",
                                 aggregate_sites = TRUE,
                                 na.rm = TRUE) 
            {stopifnot(is.data.frame(data),
            is.character(response), length(response) == 1,
            is.list(themes),
            date_col %in% names(data),
            response %in% names(data))
  
  winners <- character(0)
  all_results <- data.frame()
  
  is_anom <- function(x) grepl(anom_suffix, x, ignore.case = TRUE)
  base_name <- function(x) sub(anom_suffix, "", x, ignore.case = TRUE)
  
  .coerce_date <- function(x) {
    if (inherits(x, c("Date", "POSIXt"))) return(as.Date(x))
    as.Date(x)
  }
  
  .plot_theme_vars <- function(df, theme, vars, response, date_col, site_col,
                               aggregate_sites = TRUE, na.rm = TRUE) {
    
    df[[date_col]] <- .coerce_date(df[[date_col]])
    
    if (isTRUE(aggregate_sites) && site_col %in% names(df)) {
      df_sum <- df |>
        dplyr::group_by(.data[[date_col]]) |>
        dplyr::summarise(
          !!response := if (na.rm) sum(.data[[response]], na.rm = TRUE) else sum(.data[[response]]),
          dplyr::across(dplyr::all_of(vars), ~ if (na.rm) mean(.x, na.rm = TRUE) else mean(.x)),
          .groups = "drop"
        )} else {
      df_sum <- df |>
        dplyr::select(dplyr::all_of(c(date_col, response, vars)))
    }
    
    p_top <- ggplot2::ggplot(df_sum, ggplot2::aes(x = .data[[date_col]], y = .data[[response]])) +
      ggplot2::geom_col() +
      ggplot2::labs(
        title = paste0("Theme: ", theme),
        y = response,
        x = NULL
      ) +
      ggplot2::theme_bw()
    
    df_long <- df_sum |>
      dplyr::select(dplyr::all_of(c(date_col, vars))) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(vars),
        names_to = "variable",
        values_to = "value"
      )
    
    df_long <- df_long |>
      dplyr::group_by(.data$variable) |>
      dplyr::mutate(
        value_z = {
          m <- mean(.data$value, na.rm = TRUE)
          s <- stats::sd(.data$value, na.rm = TRUE)
          if (is.na(s) || s == 0) rep(0, dplyr::n()) else (.data$value - m) / s
        }
      ) |>
      dplyr::ungroup()
    
    p_bottom <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = .data[[date_col]], y = .data$value_z, color = .data$variable)
    ) +
      ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
      ggplot2::labs(
        y = "Env vars (z-scored)",
        x = date_col,
        color = NULL
      ) +
      ggplot2::theme_bw()
    
    p_top / p_bottom + patchwork::plot_layout(heights = c(1, 1.3))
  }
  
  for (theme in names(themes)) {
    vars <- themes[[theme]]
    vars <- vars[vars %in% names(data)]
    if (length(vars) == 0) next
    
    if (isTRUE(prefer_anoms)) {
      bases_with_anom <- unique(base_name(vars[is_anom(vars)]))
      drop_raw <- vars[!is_anom(vars) & (vars %in% bases_with_anom)]
      vars <- setdiff(vars, drop_raw)
    }
    
    theme_results <- data.frame()
    
    for (v in vars) {
      fml <- stats::as.formula(paste0(response, " ~ s(", v, ", k=4)"))
      
      mod <- mgcv::gam(
        formula = fml,
        data = data,
        family = family,
        method = "REML"
      )
      
      summ <- summary(mod)
      
      if (!is.null(summ$s.table) && nrow(summ$s.table) >= 1) {
        edf <- summ$s.table[1, "edf"]
        F   <- summ$s.table[1, "F"]
        p   <- summ$s.table[1, "p-value"]
      } else {
        edf <- NA_real_
        F   <- NA_real_
        p   <- NA_real_
      }
      
      theme_results <- rbind(
        theme_results,
        data.frame(
          theme = theme,
          variable = v,
          is_anom = is_anom(v),
          base_variable = base_name(v),
          AIC = stats::AIC(mod),
          dev_expl = summ$dev.expl,
          edf = edf,
          F = F,
          p_value = p,
          REML = mod$gcv.ubre,
          stringsAsFactors = FALSE
        )
      )
    }
    
    theme_results <- theme_results[order(theme_results$AIC), ]
    
    best <- theme_results$variable[1]
    
    cat("\nTheme:", theme,
        "\n  Winner:", best,
        "\n  AIC:", round(theme_results$AIC[1], 2),
        "\n  DevExpl:", round(theme_results$dev_expl[1], 3), "\n")
    
    winners <- c(winners, best)
    all_results <- rbind(all_results, theme_results)
    
    if (isTRUE(plot)) {
      p <- .plot_theme_vars(
        df = data,
        theme = theme,
        vars = vars,
        response = response,
        date_col = date_col,
        site_col = site_col,
        aggregate_sites = aggregate_sites,
        na.rm = na.rm
      )
      
      print(p)
      
      if (!is.null(plot_dir)) {
        if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
        fn <- file.path(plot_dir, paste0(plot_prefix, "_", theme, ".png"))
        ggplot2::ggsave(fn, p, width = 11, height = 7, dpi = 300)
      }
    }
  }
  
  list(
    winners = winners,
    results = all_results
  )
}

predictor_results <- list()
theme_results <- list()
winner_results <- list()

for (site in names(binned_deseasoned)) {
  current_df <- binned_deseasoned[[site]]
  
  current_pred <- build_preds(
    df = current_df,
    species = species,
    depths = depths,
    drop_depths = drop_depths,
    include_ice_predictors = include_ice_predictors
  )
  
  current_pred <- remove_zero_var(current_df, current_pred)
  
  current_vif <- vif_stepwise_drop(
    data = current_df,
    response = species,
    predictors = current_pred,
    vif_thresh = 5,
    family = gaussian(),
    verbose = TRUE
  )
  
  predictor_results[[site]] <- current_vif
  
  print(current_vif$kept_predictors)
  print(car::vif(current_vif$final_fit))
  
  current_themes <- assign_themes(current_vif$kept_predictors, themes_keywords)
  theme_results[[site]] <- current_themes
  
  current_winners <- select_best_by_theme(
    data = current_df,
    response = species,
    themes = current_themes
  )
  
  winner_results[[site]] <- current_winners
}

save_vif_outputs <- function(predictor_results, species_name, save_dir) {
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  vif_steps_all <- list()
  vif_kept_all <- list()
  vif_dropped_all <- list()
  
  for (site in names(predictor_results)) {
    res <- predictor_results[[site]]
    
    # Kept predictors
    kept_tbl <- tibble(
      site = site,
      species = species_name,
      predictor = res$kept_predictors
    )
    vif_kept_all[[site]] <- kept_tbl
    
    # Dropped predictors
    dropped_tbl <- tibble(
      site = site,
      species = species_name,
      predictor = res$dropped_predictors
    )
    vif_dropped_all[[site]] <- dropped_tbl
    
    # Step log
    step_tbl <- purrr::map_dfr(res$steps, function(x) {
      tibble(
        site = site,
        species = species_name,
        step = x$step,
        action = x$action,
        dropped = ifelse(is.null(x$dropped), NA_character_, x$dropped),
        worst = ifelse(is.null(x$worst), NA_character_, x$worst),
        max_vif = ifelse(is.null(x$max_vif), NA_real_, x$max_vif)
      )
    })
    
    vif_steps_all[[site]] <- step_tbl
    
    # Save per-site txt summary
    txt_lines <- c(
      paste0("Site: ", site),
      paste0("Species: ", species_name),
      "",
      paste0("Kept predictors: ", paste(res$kept_predictors, collapse = ", ")),
      paste0("Dropped predictors: ", paste(res$dropped_predictors, collapse = ", ")),
      "",
      paste0("Final formula: ", deparse(res$final_formula))
    )
    
    writeLines(
      txt_lines,
      con = file.path(save_dir, paste0("VIF_", species_name, "_", site, "_summary.txt"))
    )
  }
  
  write.csv(bind_rows(vif_steps_all),
            file = file.path(save_dir, paste0("VIF_", species_name, "_steps_all_sites.csv")),
            row.names = FALSE)
  
  write.csv(bind_rows(vif_kept_all),
            file = file.path(save_dir, paste0("VIF_", species_name, "_kept_predictors_all_sites.csv")),
            row.names = FALSE)
  
  write.csv(bind_rows(vif_dropped_all),
            file = file.path(save_dir, paste0("VIF_", species_name, "_dropped_predictors_all_sites.csv")),
            row.names = FALSE)
}

save_vif_outputs(
  predictor_results = predictor_results,
  species_name = species,
  save_dir = precheck_save_dir
)

save_theme_assignments <- function(theme_results, species_name, save_dir) {
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  theme_tbl_all <- list()
  
  for (site in names(theme_results)) {
    current_themes <- theme_results[[site]]
    
    if (length(current_themes) == 0) next
    
    current_tbl <- purrr::imap_dfr(current_themes, function(vars, theme_name) {
      tibble(
        site = site,
        species = species_name,
        theme = theme_name,
        predictor = vars
      )
    })
    
    theme_tbl_all[[site]] <- current_tbl
  }
  
  theme_tbl_all <- bind_rows(theme_tbl_all)
  
  write.csv(
    theme_tbl_all,
    file = file.path(save_dir, paste0("Themes_", species_name, "_assigned_predictors_all_sites.csv")),
    row.names = FALSE
  )
}

save_theme_assignments(
  theme_results = theme_results,
  species_name = species,
  save_dir = precheck_save_dir
)

save_theme_selection_outputs <- function(winner_results, species_name, save_dir) {
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_candidate_results <- list()
  all_winners <- list()
  
  for (site in names(winner_results)) {
    res <- winner_results[[site]]
    
    if (!is.null(res$results) && nrow(res$results) > 0) {
      candidate_tbl <- res$results %>%
        mutate(site = site, species = species_name) %>%
        dplyr::relocate(site, species)
      
      all_candidate_results[[site]] <- candidate_tbl
      
      write.csv(
        candidate_tbl,
        file = file.path(save_dir, paste0("ThemeSelection_", species_name, "_", site, "_candidate_results.csv")),
        row.names = FALSE
      )
    }
    
    if (length(res$winners) > 0) {
      winner_tbl <- tibble(
        site = site,
        species = species_name,
        winner = res$winners
      )
      
      all_winners[[site]] <- winner_tbl
    }
  }
  
  if (length(all_candidate_results) > 0) {
    write.csv(
      bind_rows(all_candidate_results),
      file = file.path(save_dir, paste0("ThemeSelection_", species_name, "_all_sites_candidate_results.csv")),
      row.names = FALSE
    )
  }
  
  if (length(all_winners) > 0) {
    write.csv(
      bind_rows(all_winners),
      file = file.path(save_dir, paste0("ThemeSelection_", species_name, "_all_sites_winners.csv")),
      row.names = FALSE
    )
  }
}

save_theme_selection_outputs(
  winner_results = winner_results,
  species_name = species,
  save_dir = precheck_save_dir
)
# =========================================================
# STEP 8: BUILD GAMS
# =========================================================

plotGam <- function(gam) {
  return(plot(gam, trans = plogis, shift = coef(gam)[1], scheme = 2, seWithMean = TRUE))
}
 
plotGam1 <- function(gam) {
  return(plot(gam, trans = plogis, shift = coef(gam)[1], seWithMean = TRUE, scheme = 2, pages = 1))
}

build_gam_formula <- function(response, predictors, k = 4) {
  smooth_terms <- paste0("s(", predictors, ", k=", k, ")")
  as.formula(paste(response, "~", paste(smooth_terms, collapse = " + ")))
}

auto_gam <- function(data, response, predictors,
                     family = tw(link = "log", a = 1.1, b = 1.9),
                     k = 4,
                     p_thresh = 0.05,
                     verbose = TRUE) {
  
  current_preds <- predictors
  
  repeat {
    form <- build_gam_formula(response, current_preds, k)
    mod <- gam(form, data = data, family = family, method = "REML", select = TRUE)
    
    summ <- summary(mod)
    pvals <- summ$s.pv
    names(pvals) <- rownames(summ$s.table)
    
    worst_p <- max(pvals)
    
    if (verbose) {
      message("Max p-value = ", round(worst_p, 4))
    }
    
    if (worst_p <= p_thresh || length(current_preds) == 1) {
      break
    }
    
    worst_term <- names(which.max(pvals))
    worst_var <- gsub("s\\(|\\)", "", worst_term)
    
    if (verbose) {
      message("Dropping: ", worst_var)
    }
    
    current_preds <- setdiff(current_preds, worst_var)
  }
  
  return(mod)
}

final_models <- list()

for (site in names(binned_deseasoned)) {
  final_models[[site]] <- auto_gam(
    data = binned_deseasoned[[site]],
    response = species,
    predictors = winner_results[[site]]$winners
  )
}

# =========================================================
# STEP 9: VISUALIZE GAMS
# =========================================================
gam_save_dir <- file.path(plot_save_dir, "final_gams")
dir.create(gam_save_dir, recursive = TRUE, showWarnings = FALSE)

save_gam_outputs <- function(gam, sp, species_name, save_dir) {
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  base_name <- paste0("GAM_", species_name, "_", sp)
  
  summ <- summary(gam)
  
  # Full text summary
  summary_txt <- capture.output(print(summ))
  writeLines(
    summary_txt,
    con = file.path(save_dir, paste0(base_name, "_summary.txt"))
  )
  
  # Smooth-term table
  if (!is.null(summ$s.table)) {
    smooth_tbl <- as.data.frame(summ$s.table)
    smooth_tbl$term <- rownames(summ$s.table)
    smooth_tbl <- dplyr::relocate(smooth_tbl, term)
    
    write.csv(
      smooth_tbl,
      file = file.path(save_dir, paste0(base_name, "_smooth_terms.csv")),
      row.names = FALSE
    )
  }
  
  # Parametric coefficients table
  if (!is.null(summ$p.table)) {
    param_tbl <- as.data.frame(summ$p.table)
    param_tbl$term <- rownames(summ$p.table)
    param_tbl <- dplyr::relocate(param_tbl, term)
    
    write.csv(
      param_tbl,
      file = file.path(save_dir, paste0(base_name, "_parametric_terms.csv")),
      row.names = FALSE
    )
  }
  
  # One-row model metrics table
  metrics_tbl <- data.frame(
    site = sp,
    species = species_name,
    n = nrow(gam$model),
    dev_expl = summ$dev.expl,
    adj_r_sq = if (!is.null(summ$r.sq)) summ$r.sq else NA_real_,
    scale = summ$scale,
    n_smooth_terms = length(gam$smooth)
  )
  
  write.csv(
    metrics_tbl,
    file = file.path(save_dir, paste0(base_name, "_metrics.csv")),
    row.names = FALSE
  )
}

visualizeGAM <- function(gam, sp,
                         species_name = species,
                         save_plot = TRUE,
                         save_outputs = TRUE,
                         save_dir = gam_save_dir,
                         width = 10, height = 8, dpi = 300) {
  
  predictors <- gam$smooth |>
    sapply(function(x) x$term) |>
    unlist()
  
  plot_info <- smooth_estimates(gam) |> add_confint()
  
  dont_shift <- names(plot_info) %in%
    c(".smooth", ".type", ".by", ".se", predictors)
  
  plot_info <- plot_info |>
    gratia:::shift_values(
      i = dont_shift,
      h = coef(gam)[1],
      FUN = "+"
    ) |>
    transform_fun(fun = plogis)
  
  summ <- summary(gam)
  deviance <- round(summ$dev.expl * 100, 2)
  p_values <- setNames(summ$s.pv, rownames(summ$s.table))
  
  all_plots <- list()
  
  for (p in predictors) {
    current_plot <- dplyr::filter(
      plot_info,
      .smooth == paste0("s(", p, ")")
    )
    
    current_p_val <- p_values[[paste0("s(", p, ")")]]
    current_p_val <- max(current_p_val, 1e-6)
    
    current_plot$label <- paste0("p = ", signif(current_p_val, 3))
    
    x_vals <- current_plot[[p]]
    x_lim <- range(x_vals, na.rm = TRUE)
    
    plot <- ggplot(current_plot) +
      geom_ribbon(
        aes(x = .data[[p]], ymin = .lower_ci, ymax = .upper_ci),
        alpha = 0.2
      ) +
      geom_line(
        aes(x = .data[[p]], y = .estimate),
        linewidth = 1
      ) +
      geom_rug(
        data = gam$model,
        aes(x = .data[[p]]),
        sides = "b"
      ) +
      labs(
        y = "Partial effect",
        x = nameVar(p),
        subtitle = paste0(nameVar(p), " (", current_plot$label[1], ")")
      ) +
      theme_bw() +
      ylim(0, 1) +
      xlim(x_lim)
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  final_plot <- patchwork::wrap_plots(all_plots) +
    patchwork::plot_annotation(
      title = paste0(
        species_name, " at ", sp,
        " (", deviance, "% deviance explained)"
      )
    )
  
  print(final_plot)
  
  if (isTRUE(save_plot) || isTRUE(save_outputs)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
  }
  
  if (isTRUE(save_plot)) {
    file_name <- paste0("GAM_", species_name, "_", sp, ".png")
    out_file <- file.path(save_dir, file_name)
    
    ggplot2::ggsave(
      filename = out_file,
      plot = final_plot,
      width = width,
      height = height,
      dpi = dpi
    )
    
    message("Plot saved to: ", out_file)
  }
  
  if (isTRUE(save_outputs)) {
    save_gam_outputs(
      gam = gam,
      sp = sp,
      species_name = species_name,
      save_dir = save_dir
    )
  }
  
  return(final_plot)
}

nameVar <- function(var) {
  labels <- c(
    julian_day = "Julian Day",
    SSH = "Sea Surface Height (m)",
    SSH_anom = "De-seasoned Sea Surface Height (m)",
    FSLE = "FSLE Magnitude",
    FSLE_anom = "De-seasoned FSLE Magnitude",
    fsle_orient = "FSLE Orientation",
    fsle_orient_anom = "De-seasoned FSLE Orientation",
    mixed_layer = "Mixed Layer Depth (m)",
    mixed_layer_anom = "De-seasoned Mixed Layer Depth (m)",
    ice_conc = "Sea Ice Concentration",
    ice_conc_anom = "De-seasoned Sea Ice Concentration",
    ice_thickness = "Sea Ice Thickness (m)",
    ice_thickness_anom = "De-seasoned Sea Ice Thickness (m)"
  )
  
  if (var %in% names(labels)) return(labels[[var]])
  
  depth_suffix_label <- function(d) {
    if (d == "0") return("Sea Surface")
    paste0("@ ", d, "m")
  }
  
  if (grepl("_\\d+mon$", var)) {
    base <- sub("_(\\d+)mon$", "", var)
    lagN <- sub("^.*_(\\d+)mon$", "\\1", var)
    
    base_label <- switch(
      base,
      chla = "Chlorophyll (mg/m³)",
      productivity = "Net Primary Production (mg/m³/day C)",
      temperature = "Temperature (°C)",
      salinity = "Salinity (psu)",
      o2 = "Oxygen (mmol/m³)",
      EKE = "Eddy Kinetic Energy",
      base
    )
    return(paste0(lagN, " Month Lag: ", base_label))
  }
  
  if (grepl("^(temp|temperature|salinity|o2|chla|productivity|EKE|SSH|FSLE|mixed_layer|fsle_orient)_(sd|mad)_\\d+$", var)) {
    base <- sub("_(sd|mad)_\\d+$", "", var)
    stat <- sub("^.*_(sd|mad)_\\d+$", "\\1", var)
    depth <- sub("^.*_(sd|mad)_(\\d+)$", "\\2", var)
    
    base_label <- switch(
      base,
      temp = "Temperature (°C)",
      temperature = "Temperature (°C)",
      salinity = "Salinity (psu)",
      o2 = "Oxygen (mmol/m³)",
      chla = "Chlorophyll (mg/m³)",
      productivity = "Net Primary Production (mg/m³/day C)",
      EKE = "Eddy Kinetic Energy",
      SSH = "Sea Surface Height (m)",
      FSLE = "FSLE Magnitude",
      mixed_layer = "Mixed Layer Depth (m)",
      fsle_orient = "FSLE Orientation",
      base
    )
    
    stat_label <- if (stat == "sd") "SD" else "MAD"
    return(paste0(stat_label, " of ", base_label, " ", depth_suffix_label(depth)))
  }
  
  if (grepl("_(sd|mad)$", var)) {
    base <- sub("_(sd|mad)$", "", var)
    stat <- sub("^.*_(sd|mad)$", "\\1", var)
    
    base_label <- switch(
      base,
      fsle = "FSLE Magnitude",
      FSLE = "FSLE Magnitude",
      ssh = "Sea Surface Height (m)",
      SSH = "Sea Surface Height (m)",
      mixed_layer = "Mixed Layer Depth (m)",
      fsle_orient = "FSLE Orientation",
      ice_conc = "Sea Ice Concentration",
      ice_thickness = "Sea Ice Thickness (m)",
      base
    )
    
    stat_label <- if (stat == "sd") "SD" else "MAD"
    return(paste0(stat_label, " of ", base_label))
  }
  
  if (grepl("_\\d+$", var)) {
    base <- sub("_(\\d+)$", "", var)
    depth <- sub("^.*_(\\d+)$", "\\1", var)
    
    base_label <- switch(
      base,
      temperature = "Temperature (°C)",
      temp = "Temperature (°C)",
      salinity = "Salinity (psu)",
      o2 = "Oxygen (mmol/m³)",
      chla = "Chlorophyll (mg/m³)",
      productivity = "Net Primary Production (mg/m³/day C)",
      EKE = "Eddy Kinetic Energy",
      base
    )
    
    return(paste0(base_label, " ", depth_suffix_label(depth)))
  }
  
  if (grepl("_anom$", var)) {
    base <- sub("_anom$", "", var)
    return(paste0("De-seasoned ", nameVar(base)))
  }
  
  var
}
gam_plots_final <- list()

for (site in names(final_models)) {
  gam_plots_final[[site]] <- visualizeGAM(
    gam = final_models[[site]],
    sp = site,
    species_name = species,
    save_plot = TRUE,
    save_outputs = TRUE,
    save_dir = gam_save_dir
  )
}

# =========================================================
# STEP 10: MODEL CHECKING
# =========================================================

diag_save_dir <- file.path(plot_save_dir, "diagnostics")
dir.create(diag_save_dir, recursive = TRUE, showWarnings = FALSE)

acf_df <- function(x, lag.max = 200) {
  a <- acf(x, plot = FALSE, lag.max = lag.max, na.action = na.pass)
  tibble(lag = as.numeric(a$lag), acf = as.numeric(a$acf))
}

run_gam_check_base <- function(mod,
                               model_name = "GAM",
                               do_gam_check = TRUE,
                               do_concurvity = TRUE,
                               save_txt = FALSE,
                               txt_file = NULL) {
  
  out_lines <- character(0)
  add_line <- function(...) {
    txt <- paste0(...)
    out_lines <<- c(out_lines, txt)
    cat(txt, "\n")
  }
  
  add_line("==============================")
  add_line("Model: ", model_name)
  add_line("==============================")
  
  if (isTRUE(do_gam_check)) {
    add_line("")
    add_line("--- mgcv::gam.check() ---")
    
    old_par <- par(no.readonly = TRUE)
    gam_check_txt <- capture.output(
      suppressWarnings(gam.check(mod))
    )
    par(old_par)
    
    cat(paste(gam_check_txt, collapse = "\n"), "\n")
    out_lines <- c(out_lines, gam_check_txt)
  }
  
  if (isTRUE(do_concurvity)) {
    add_line("")
    add_line("--- mgcv::concurvity(full = TRUE) ---")
    
    conc <- mgcv::concurvity(mod, full = TRUE)
    
    conc_est <- if (is.list(conc) && !is.null(conc$estimate)) {
      conc$estimate
    } else {
      conc
    }
    
    conc_txt <- capture.output(print(conc_est))
    cat(paste(conc_txt, collapse = "\n"), "\n")
    out_lines <- c(out_lines, conc_txt)
    
    conc_tbl <- as.data.frame(conc_est)
    conc_tbl$term <- rownames(conc_est)
    conc_tbl <- dplyr::relocate(conc_tbl, term)
    
    add_line("")
    add_line("Sorted by worst concurvity:")
    
    if ("worst" %in% names(conc_tbl)) {
      sorted_txt <- capture.output(print(conc_tbl %>% arrange(desc(worst))))
    } else {
      sorted_txt <- capture.output(print(conc_tbl))
    }
    
    cat(paste(sorted_txt, collapse = "\n"), "\n")
    out_lines <- c(out_lines, sorted_txt)
  }
  
  if (isTRUE(save_txt) && !is.null(txt_file)) {
    writeLines(out_lines, txt_file)
  }
  
  invisible(out_lines)
}

plot_gam_diagnostics <- function(mod,
                                 model_name = "GAM",
                                 date = NULL,
                                 lag.max = 200,
                                 do_appraise = TRUE,
                                 save_plots = TRUE,
                                 save_dir = diag_save_dir,
                                 width = 10,
                                 height = 8,
                                 dpi = 300) {
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  p_appraise <- NULL
  p_acf <- NULL
  p_time <- NULL
  
  if (isTRUE(do_appraise)) {
    p_appraise <- appraise(mod) +
      plot_annotation(title = paste0(model_name, " | gratia::appraise()"))
    print(p_appraise)
    
    if (isTRUE(save_plots)) {
      ggplot2::ggsave(
        filename = file.path(save_dir, paste0(model_name, "_appraise.png")),
        plot = p_appraise,
        width = width,
        height = height,
        dpi = dpi
      )
    }
  }
  
  r <- residuals(mod, type = "deviance")
  
  if (!is.null(date)) {
    ord <- order(as.Date(date))
    r_use <- r[ord]
  } else {
    r_use <- r
  }
  
  acf_dat <- acf_df(r_use, lag.max = lag.max)
  
  p_acf <- ggplot(acf_dat, aes(lag, acf)) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(xend = lag, yend = 0)) +
    labs(
      title = paste0(model_name, " | Deviance residual ACF"),
      x = "Lag",
      y = "ACF"
    ) +
    theme_bw()
  
  print(p_acf)
  
  if (isTRUE(save_plots)) {
    ggplot2::ggsave(
      filename = file.path(save_dir, paste0(model_name, "_acf.png")),
      plot = p_acf,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  if (!is.null(date)) {
    df_time <- tibble(date = as.Date(date), resid = r) %>% arrange(date)
    
    p_time <- ggplot(df_time, aes(date, resid)) +
      geom_hline(yintercept = 0) +
      geom_line() +
      labs(
        title = paste0(model_name, " | Deviance residuals vs time"),
        x = NULL,
        y = "Deviance residual"
      ) +
      theme_bw()
    
    print(p_time)
    
    if (isTRUE(save_plots)) {
      ggplot2::ggsave(
        filename = file.path(save_dir, paste0(model_name, "_resid_time.png")),
        plot = p_time,
        width = width,
        height = height,
        dpi = dpi
      )
    }
    
    if (isTRUE(save_plots)) {
      write.csv(
        df_time,
        file = file.path(save_dir, paste0(model_name, "_residuals.csv")),
        row.names = FALSE
      )
    }
  }
  
  if (isTRUE(save_plots)) {
    write.csv(
      acf_dat,
      file = file.path(save_dir, paste0(model_name, "_acf_values.csv")),
      row.names = FALSE
    )
  }
  
  invisible(list(
    appraise_plot = p_appraise,
    acf_plot = p_acf,
    time_plot = p_time,
    acf_data = acf_dat
  ))
}

diagnostic_results <- list()

for (site in names(final_models)) {
  model_name <- paste0(site, "_final_", species)
  
  txt_file <- file.path(diag_save_dir, paste0(model_name, "_checks.txt"))
  
  run_gam_check_base(
    mod = final_models[[site]],
    model_name = model_name,
    do_gam_check = TRUE,
    do_concurvity = TRUE,
    save_txt = TRUE,
    txt_file = txt_file
  )
  
  diagnostic_results[[site]] <- plot_gam_diagnostics(
    mod = final_models[[site]],
    model_name = model_name,
    date = binned_deseasoned[[site]]$bin_start,
    lag.max = 200,
    do_appraise = TRUE,
    save_plots = TRUE,
    save_dir = diag_save_dir
  )
}
