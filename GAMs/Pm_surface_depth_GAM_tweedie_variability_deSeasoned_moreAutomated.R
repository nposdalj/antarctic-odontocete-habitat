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
# MODELS THAT HAVE BEEN IMPROVED: CI
# MODELS IN INITIAL PASSTHROUGH: EI, KGI

# ------------- Step 0: Choose Species ----------------
plot_save_dir <- "F:/Antarctica/GAMs"
# Modeling Gm for all sites, 40 km radius environmental data
species <- c('Pm') # options: BW29, BW37, Oo, Pm, Gm
# BW29 = Southern bottlenose whale, BW37 = Gray's and strap-toothed whales
# Oo = Killer whale, Pm = Sperm Whale, Gm = Long-finned pilot whale
# Note: not enough data to model BW58 at any site
sites <- c('EI','KGI','CI')
#ice_threshold <- 0.10

themes_keywords <- list(
  mesoscale = c("FSLE","SSH","ssh","EKE"),
  front_orientation = c("fsle_orient"),
  productivity = c("productivity","chla"),
  oxygen = c("o2"),
  salinity = c("salinity"),
  stratification = c("mixed_layer"),
  temperature = c("temp","temperature")
)

# Function to write out full name of a species/site code
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
  } else {
    fullname <- "Sperm Whale"
  }
  return(fullname)
}
# ------------- Step 1: Load Data -----------------
allData <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/allData.csv")
allData <- allData %>% subset(select=-X)
allData$date <- as.Date(allData$date, "%Y-%m-%d")
# Filter by species relevant data
# Only adding standard deviations of surface variables, feel free to change that if needed
depths <- c(0, 375, 1665) # Gm depths

# species is assumed to be a string like "Gm"
# (same as in your code where you used get(species))
all_species <- c("BW29","BW37","BW58","Gm","Pm","Oo")
drop_species <- setdiff(all_species, species)

# keep columns that end with _0, _16, or _1655
keep_depth_regex <- paste0("_(", paste(depths, collapse="|"), ")$")

sp_specific <- allData %>%
  dplyr::select(-dplyr::any_of(drop_species)) %>%
  dplyr::select(
    date, Site, julian_day, dplyr::all_of(species),
    dplyr::matches(keep_depth_regex),
    AAO:last_col()
  )

# Intermediate step to find ACF for seasonal model only, need to go back and adjust Trisha's ACF function
allData_EI <- allData %>% 
  filter(Site == "EI")
BlockMod_EI <- gam(Pm ~ s(SSH,k=4,sp=0.1) + s(chla_0,k=4,sp=0.1) + s(salinity_0,k=4,sp=0.1) + s(o2_375,k=4,sp=0.1),
                   family = tw(link = "log", a = 1.1, b = 1.9), data = allData_EI, method = "REML")
ACF = acf(residuals(BlockMod_EI), lag.max = 1500) 
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval_EI = ACFidx[1]

allData_CI <- allData %>% 
  filter(Site == "CI")
BlockMod_CI <- gam(Pm ~ s(SSH,k=4,sp=0.1) + s(chla_0,k=4,sp=0.1) + s(salinity_375,k=4) + s(o2_375,k=4) + s(chla_375,k=4),
                   family = tw(link = "log", a = 1.1, b = 1.9), data = allData_CI, method = "REML")
ACF = acf(residuals(BlockMod_CI), lag.max = 1500) 
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval_CI = ACFidx[1]

# ------------- Step 2: Average by ACF ------------
acf_table <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Autocorrelation/acf_table.csv")
acfVal <- function(site) {
  row_idx <- which(acf_table$site == site) # Row index for the site
  acf_val <- acf_table[row_idx,species][[1]]
  return(acf_val)
}

binByACF <- function(data, site, bin_days, species, depths,
                     id_cols = c("date","Site"),
                     group_cols = c("bin_start","Site"),
                     exclude_cols = c("X","X.1")) {
  
  # ---- basic checks ----
  stopifnot(is.data.frame(data))
  stopifnot(length(species) == 1)
  
  # Filter by site and make bin_start
  sp_filtered <- data %>%
    filter(.data$Site == site) %>%
    mutate(bin_start = floor_date(.data$date, unit = paste(bin_days, "days")))
  
  # Helper to safely create mean summaries only for columns that exist
  mean_col <- function(col) rlang::expr(mean(.data[[!!col]], na.rm = TRUE))
  
  # ---- build column sets dynamically ----
  
  # 1) Depth columns to keep/summarize: anything ending in _{depth}
  depth_regex <- paste0("_(", paste(depths, collapse = "|"), ")$")
  
  depth_cols <- names(sp_filtered)[grepl(depth_regex, names(sp_filtered))]
  
  # 2) Species column
  species_cols <- intersect(species, names(sp_filtered))  # should be 1 if present
  
  # 3) "Non-depth" predictor columns:
  #    take everything EXCEPT: ids, other species labels, depth cols, and excluded cols
  all_species <- c("BW29","BW37","BW58","Gm","Pm","Oo")
  other_species <- setdiff(all_species, species)
  
  non_depth_cols <- setdiff(
    names(sp_filtered),
    c(id_cols, other_species, exclude_cols, depth_cols)
  )
  
  # If you only want *numeric* non-depth columns summarized by mean, filter them:
  non_depth_cols <- non_depth_cols[sapply(sp_filtered[non_depth_cols], is.numeric)]
  
  # (Optional) ensure julian_day always included if present
  if ("julian_day" %in% names(sp_filtered) && !"julian_day" %in% non_depth_cols) {
    non_depth_cols <- c("julian_day", non_depth_cols)
  }
  
  # ---- build summarise expressions ----
  
  # Species: mean of raw values per bin (as you requested)
  species_expr <- set_names(
    lapply(species_cols, \(x) rlang::expr(mean(.data[[x]], na.rm = TRUE))),
    species_cols
  )
  
  # Non-depth: mean for each
  non_depth_expr <- set_names(
    lapply(non_depth_cols, mean_col),
    non_depth_cols
  )
  
  # Depth cols: mean for each
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

EI_binned <- binByACF(
  data     = sp_specific,
  site     = "EI",
  bin_days = ACFval_EI,
  species  = species,
  depths   = depths
)

CI_binned <- binByACF(
  data     = sp_specific,
  site     = "CI",
  bin_days = ACFval_CI,
  species  = species,
  depths   = depths
)


# ------------- Step 3: Plot Presence Timeseries --------------
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

# Remove days with less than a certain percentage of sea ice concentration
get_ice_thresholds <- function(df, response = "Pm",
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

thr_CI <- get_ice_thresholds(CI_binned, response = "Pm", q = 0.95)
thr_EI <- get_ice_thresholds(EI_binned, response = "Pm", q = 0.95)

# apply to one dataframe
EI_binned <- EI_binned %>%
  filter(ice_conc < thr_EI$thr)

CI_binned <- CI_binned %>%
  filter(ice_conc < thr_CI$thr)

# Deseason data
vars_to_anom <- c(
  # Large-scale / physical
  "SSH",
  "mixed_layer",
  "FSLE",
  "fsle_orient",
  
  # Surface (0 m)
  "temperature_0",
  "salinity_0",
  "o2_0",
  "chla_0",
  "productivity_0",
  "EKE_0",
  
  # Shallow (375 m)
  "temperature_375",
  "salinity_375",
  "o2_375",
  "chla_375",
  "productivity_375",
  "EKE_375",
  
  # Deep (1655 m)
  "temperature_1655",
  "salinity_1655",
  "o2_1655",
  "chla_1655",
  "productivity_1655",
  "EKE_1655"
)

make_monthly_anoms_binned <- function(
    df,
    vars,
    date_col = "bin_start",
    site_col = "Site",
    standardize = FALSE
) {
  stopifnot(date_col %in% names(df), site_col %in% names(df))
  
  # Coerce date if needed
  if (!inherits(df[[date_col]], c("Date", "POSIXt"))) {
    df[[date_col]] <- as.Date(df[[date_col]])
    if (any(is.na(df[[date_col]]))) {
      stop(sprintf("Column '%s' could not be coerced to Date.", date_col))
    }
  }
  
  # Keep only vars that exist
  vars <- intersect(vars, names(df))
  if (length(vars) == 0) stop("None of the requested vars are present in the data.")
  
  df_out <- df |>
    dplyr::mutate(.month = lubridate::month(.data[[date_col]])) |>
    dplyr::group_by(.data[[site_col]], .month) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(vars),
        ~ .x - mean(.x, na.rm = TRUE),
        .names = "{.col}_anom"
      )
    )
  
  if (isTRUE(standardize)) {
    df_out <- df_out |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(vars),
          ~ {
            s <- stats::sd(.x, na.rm = TRUE)
            if (all(is.na(.x)) || is.na(s) || s == 0) 0 else (.x - mean(.x, na.rm = TRUE)) / s
          },
          .names = "{.col}_z"
        )
      )
  }
  
  df_out |>
    dplyr::ungroup() |>
    dplyr::select(-.month)
}

EI_binned_deseasoned <- make_monthly_anoms_binned(EI_binned, vars_to_anom, date_col = "bin_start", standardize = FALSE)
CI_binned_deseasoned  <- make_monthly_anoms_binned(CI_binned,  vars_to_anom, date_col = "bin_start", standardize = FALSE)

# -------------- Step 4: VIF for Correlation -------------------
vif_stepwise_drop <- function(data, response, predictors,
                              vif_thresh = 10,
                              family = gaussian(),
                              max_iter = 200,
                              verbose = TRUE) {
  # keep only predictors that exist
  preds <- intersect(predictors, names(data))
  if (length(preds) < 2) stop("Need at least 2 predictors present in `data`.")
  
  # helper: safe VIF extraction
  get_vif <- function(fit) {
    v <- tryCatch(car::vif(fit), error = function(e) NA)
    if (all(is.na(v))) return(v)
    # car::vif can return a matrix for multi-df terms; convert to named numeric
    if (is.matrix(v)) {
      # common columns: GVIF, Df, GVIF^(1/(2*Df))
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
  
  # fallback: pick a column to drop if model/vif fails
  fallback_drop <- function(df, preds) {
    X <- df[, preds, drop = FALSE]
    X <- X[, sapply(X, is.numeric), drop = FALSE]
    if (ncol(X) < 2) return(preds[1])
    cm <- suppressWarnings(stats::cor(X, use = "pairwise.complete.obs"))
    cm[lower.tri(cm, diag = TRUE)] <- NA
    max_pair <- which(abs(cm) == max(abs(cm), na.rm = TRUE), arr.ind = TRUE)
    if (nrow(max_pair) == 0) return(colnames(X)[1])
    # drop the variable with higher mean absolute correlation
    i <- max_pair[1, 1]; j <- max_pair[1, 2]
    a <- colnames(X)[i]; b <- colnames(X)[j]
    mean_abs_cor <- function(var) mean(abs(stats::cor(X[[var]], X, use="pairwise.complete.obs")), na.rm = TRUE)
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
    
    # If VIF couldn't be computed (or some NA), fallback-drop
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
    
    # drop the worst offender
    preds <- setdiff(preds, max_var)
    log[[iter]]$action <- "drop (max vif)"
    log[[iter]]$dropped <- max_var
  }
  
  # final fit
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
                        drop_depths = c(375),
                        keep_lags = c(3, 6),
                        chla_prod_depth_max = 50,
                        
                        # raw core vars allowed
                        core_allow = c("FSLE","SSH","mixed_layer","fsle_orient"),
                        
                        # keep anoms ONLY for these (raw still included)
                        keep_anom_core = c("FSLE","SSH","mixed_layer"),
                        keep_anom_depth_bases = c("temperature","salinity","o2"),
                        keep_anom_depths = c(0, 375, 1655),
                        
                        # keep EKE_mad_* (but still drop other *_mad_* if they exist)
                        keep_EKE_mad = TRUE
) {
  
  nms <- names(df)
  
  # --- Start: drop IDs/response/known non-predictors ---
  pred <- setdiff(
    nms,
    c("bin_start","Site", species,
      "julian_day","AAO",
      "ice_conc","ice_thickness","ice_diff","ice_regime",
      "bathymetry")
  )
  
  # --- Drop STL + "*_anomaly" products (but NOT our monthly *_anom) ---
  pred <- pred[!grepl("(_anomaly$|_stl$)", pred)]
  
  # --- Drop ALL sd variables anywhere in the name (temp_sd_0, o2_sd_1655, ssh_sd, etc.) ---
  pred <- pred[!grepl("_sd", pred, ignore.case = TRUE)]
  
  # --- MAD handling: keep only EKE_mad_* if requested; drop other *_mad_* families ---
  if (isTRUE(keep_EKE_mad)) {
    keep_mad <- pred[grepl("^EKE_mad_\\d+$", pred)]
    drop_mad <- pred[grepl("_mad", pred, ignore.case = TRUE) & !grepl("^EKE_mad_\\d+$", pred)]
    pred <- setdiff(pred, drop_mad)
  } else {
    pred <- pred[!grepl("_mad", pred, ignore.case = TRUE)]
  }
  
  # --- Lags: keep only 3mon/6mon for productivity/chla ---
  lag_vars_all <- pred[grepl("_\\d+mon$", pred)]
  lag_vars_all <- lag_vars_all[grepl("^(productivity|chla)_", lag_vars_all)]
  lag_nums <- as.integer(sub(".*_(\\d+)mon$", "\\1", lag_vars_all))
  lag_keep_vars <- lag_vars_all[lag_nums %in% keep_lags]
  
  # --- Depth handling (raw depth vars only; exclude *_anom depth vars) ---
  depth_vars_all <- pred[grepl("_\\d+$", pred)]
  depth_vars_all <- setdiff(depth_vars_all, pred[grepl("_anom$", pred)]) # exclude deseasoned cols
  
  depth_nums <- as.integer(sub(".*_(\\d+)$", "\\1", depth_vars_all))
  
  # drop specified depths (e.g., 375 m)
  keep_depth <- !(depth_nums %in% drop_depths)
  
  # remove chla/productivity for depths > chla_prod_depth_max
  is_chla_prod <- grepl("^(chla|productivity)_", depth_vars_all)
  keep_chla_prod <- !(is_chla_prod & depth_nums > chla_prod_depth_max)
  
  depth_keep_vars <- depth_vars_all[keep_depth & keep_chla_prod]
  
  # --- Core (non-depth, non-lag, non-anom) vars ---
  anom_vars_all <- pred[grepl("_anom$", pred)]
  core_vars <- setdiff(pred, c(depth_vars_all, lag_vars_all, anom_vars_all))
  core_keep <- intersect(core_vars, core_allow)
  
  # --- Keep ONLY selected *_anom variables (raw stays included) ---
  anom_keep <- character(0)
  
  # core anomalies like SSH_anom
  anom_keep <- c(anom_keep, paste0(keep_anom_core, "_anom"))
  
  # depth anomalies like temperature_0_anom, salinity_1655_anom, o2_0_anom
  depth_anom_candidates <- as.vector(outer(keep_anom_depth_bases, keep_anom_depths, paste, sep = "_"))
  anom_keep <- c(anom_keep, paste0(depth_anom_candidates, "_anom"))
  
  # only keep anomalies that exist in df
  anom_keep <- intersect(anom_keep, pred)
  
  # --- Final predictor set ---
  final <- sort(unique(c(core_keep, depth_keep_vars, lag_keep_vars, anom_keep)))
  
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
    
    # prevent fsle_orient from entering mesoscale theme
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
                                 family = tw(link="log", a=1.1, b=1.9),
                                 prefer_anoms = TRUE,
                                 anom_suffix = "_anom$",
                                 date_col = "bin_start",
                                 site_col = "Site",
                                 plot = TRUE,
                                 plot_dir = NULL,
                                 plot_prefix = "theme_compare",
                                 aggregate_sites = TRUE,
                                 na.rm = TRUE) {
  
  stopifnot(is.data.frame(data),
            is.character(response), length(response) == 1,
            is.list(themes),
            date_col %in% names(data),
            response %in% names(data))
  
  winners <- character(0)
  all_results <- data.frame()
  
  is_anom <- function(x) grepl(anom_suffix, x, ignore.case = TRUE)
  base_name <- function(x) sub(anom_suffix, "", x, ignore.case = TRUE)
  
  # plotting deps
  if (isTRUE(plot)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Need ggplot2 for plotting.")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("Need dplyr for plotting.")
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("Need tidyr for plotting.")
    if (!requireNamespace("scales", quietly = TRUE)) stop("Need scales for plotting.")
    if (!requireNamespace("patchwork", quietly = TRUE)) stop("Need patchwork for plotting.")
    if (!requireNamespace("lubridate", quietly = TRUE)) stop("Need lubridate for plotting.")
  }
  
  # helper: safe date coercion
  .coerce_date <- function(x) {
    if (inherits(x, c("Date", "POSIXt"))) return(as.Date(x))
    as.Date(x)
  }
  
  # helper: make theme comparison plot (whale bar on top, env lines below)
  .plot_theme_vars <- function(df, theme, vars, response, date_col, site_col,
                               aggregate_sites = TRUE, na.rm = TRUE) {
    
    df <- df
    df[[date_col]] <- .coerce_date(df[[date_col]])
    
    # optionally aggregate across sites per bin_start
    if (isTRUE(aggregate_sites) && site_col %in% names(df)) {
      df_sum <- df |>
        dplyr::group_by(.data[[date_col]]) |>
        dplyr::summarise(
          !!response := if (na.rm) sum(.data[[response]], na.rm = TRUE) else sum(.data[[response]]),
          dplyr::across(dplyr::all_of(vars), ~ if (na.rm) mean(.x, na.rm = TRUE) else mean(.x)),
          .groups = "drop"
        )} else {
      # if not aggregating sites, just use the df as-is (may look messy with duplicates)
      df_sum <- df |>
        dplyr::select(dplyr::all_of(c(date_col, response, vars)))
    }
    
    # top: whale response as bars
    p_top <- ggplot2::ggplot(df_sum, ggplot2::aes(x = .data[[date_col]], y = .data[[response]])) +
      ggplot2::geom_col() +
      ggplot2::labs(
        title = paste0("Theme: ", theme),
        y = response,
        x = NULL
      ) +
      ggplot2::theme_bw()
    
    # bottom: env vars as overlapping lines (z-scored so they share a scale)
    df_long <- df_sum |>
      dplyr::select(dplyr::all_of(c(date_col, vars))) |>
      tidyr::pivot_longer(cols = dplyr::all_of(vars),
                          names_to = "variable",
                          values_to = "value")
    
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
    
    p_bottom <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[date_col]], y = .data$value_z, color = .data$variable)) +
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
    
    # Optional: if both raw and *_anom exist, keep only *_anom
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
    
    # ---- NEW: plot whale bars + env lines for this theme ----
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

#EI
EI_pred <- build_preds(EI_binned_deseasoned, species = species)

remove_zero_var <- function(df, vars) {
  vars[sapply(df[vars], function(x) sd(x, na.rm = TRUE) > 0)]
}
EI_pred <- remove_zero_var(EI_binned_deseasoned, EI_pred)

res_EI <- vif_stepwise_drop(
  data       = EI_binned_deseasoned,
  response   = species,          # e.g., "Gm"
  predictors = EI_pred,
  vif_thresh = 5,
  family     = gaussian(),       # like you were doing
  verbose    = TRUE
)

res_EI$kept_predictors
car::vif(res_EI$final_fit)

themes_EI <- assign_themes(res_EI$kept_predictors, themes_keywords)
EI_winners <- select_best_by_theme(
  data = EI_binned_deseasoned,
  response = "Pm",
  themes = themes_EI
)

#CI
CI_pred <- build_preds(CI_binned_deseasoned, species = species)

remove_zero_var <- function(df, vars) {
  vars[sapply(df[vars], function(x) sd(x, na.rm = TRUE) > 0)]
}
CI_pred <- remove_zero_var(CI_binned_deseasoned, CI_pred)

res_CI <- vif_stepwise_drop(
  data       = CI_binned_deseasoned,
  response   = species,          # e.g., "Gm"
  predictors = CI_pred,
  vif_thresh = 5,
  family     = gaussian(),       # like you were doing
  verbose    = TRUE
)

res_CI$kept_predictors
car::vif(res_CI$final_fit)
themes_CI <- assign_themes(res_CI$kept_predictors, themes_keywords)
CI_winners <- select_best_by_theme(
  data = CI_binned_deseasoned,
  response = "Pm",
  themes = themes_CI
)

# -------------- Step 5: Build GAMs ------------------------
# Function to visualize GAMs on a probability scale with the proper confidence interval
# Run this for each iteration of the model to plot smooth terms
plotGam <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],scheme=2,seWithMean=TRUE))
}
# Run this if all plots should be in one figure
plotGam1 <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],seWithMean=TRUE,scheme=2,pages=1))
}

 #Auto-build GAMs
 build_gam_formula <- function(response, predictors, k = 4) {
   smooth_terms <- paste0("s(", predictors, ", k=", k, ")")
   as.formula(paste(response, "~", paste(smooth_terms, collapse = " + ")))
 }

#Backward elimination
auto_gam <- function(data, response, predictors,
                     family = tw(link="log", a=1.1, b=1.9),
                     k = 4,
                     p_thresh = 0.05,
                     verbose = TRUE) {

  current_preds <- predictors

  repeat {

    form <- build_gam_formula(response, current_preds, k)

    mod <- gam(form, data = data, family = family, method = "REML",select = TRUE)

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

EI_final <- auto_gam(
  data = EI_binned_deseasoned,
  response = "Pm",
  predictors = EI_winners$winners
)

CI_final <- auto_gam(
  data = CI_binned_deseasoned,,
  response = "Pm",
  predictors = CI_winners$winners
)

fit_gam_inference <- function(data, response, predictors,
                              family = tw(link="log", a=1.1, b=1.9),
                              k = 4) {
  
  form <- build_gam_formula(response, predictors, k)
  
  gam(
    form,
    data = data,
    family = family,
    method = "REML",
    select = TRUE   # shrinkage penalty
  )
}

EI_final2  <- fit_gam_inference(EI_binned_deseasoned,  "Pm", EI_winners$winners)
CI_final2  <- fit_gam_inference(CI_binned_deseasoned,  "Pm", CI_winners$winners)

# ------------------ Step 6: Visualize GAMs -------------------
# Function to create a cleaner visualization of a GAM model
visualizeGAM <- function(gam, sp,
                         species_name = species,
                         save_plot = TRUE,
                         save_dir = plot_save_dir,
                         width = 10, height = 8, dpi = 300) {
  
  # Automatically extract smooth variable names
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
        aes(x = .data[[p]],
            ymin = .lower_ci,
            ymax = .upper_ci),
        alpha = 0.2
      ) +
      geom_line(
        aes(x = .data[[p]],
            y = .estimate),
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
  
  if (isTRUE(save_plot)) {
    
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    
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
  
  return(final_plot)
}

# Function to generate axis names from given variable names
nameVar <- function(var) {
  
  # ---- 1) Explicit labels you already curated ----
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
    
    temperature_0 = "Sea Surface Temperature (°C)",
    temperature_0_anom = "De-seasoned Sea Surface Temperature (°C)",
    temperature_375 = "Temperature @ 375m (°C)",
    temperature_375_anom = "De-seasoned Temperature @ 375m (°C)",
    temperature_1655 = "Temperature @ 1655m (°C)",
    temperature_1655_anom = "De-seasoned Temperature @ 1655m (°C)",
    
    salinity_0 = "Sea Surface Salinity (psu)",
    salinity_0_anom = "De-seasoned Sea Surface Salinity (psu)",
    salinity_375 = "Salinity @ 375m (psu)",
    salinity_375_anom = "De-seasoned Salinity @ 375m (psu)",
    salinity_1655 = "Salinity @ 1655m (psu)",
    salinity_1655_anom = "De-seasoned Salinity @ 1655m (psu)",
    
    EKE_0 = "Eddy Kinetic Energy (0m)",
    EKE_0_anom = "De-seasoned Eddy Kinetic Energy (0m)",
    EKE_375 = "Eddy Kinetic Energy @ 375m",
    EKE_375_anom = "De-seasoned Eddy Kinetic Energy @ 375m",
    EKE_1655 = "Eddy Kinetic Energy @ 1655m",
    EKE_1655_anom = "De-seasoned Eddy Kinetic Energy @ 1655m",
    
    EKE_mad_0 = "Eddy Kinetic Energy Variability (0m)",
    EKE_mad_0_anom = "De-seasoned Eddy Kinetic Energy Variability (0m)",
    
    chla_0 = "Chlorophyll (mg/m³)",
    chla_0_anom = "De-seasoned Chlorophyll (mg/m³)",
    chla_375 = "Chlorophyll @ 375m (mg/m³)",
    chla_375_anom = "De-seasoned Chlorophyll @ 375m (mg/m³)",
    
    o2_0 = "Oxygen (mmol/m³)",
    o2_0_anom = "De-seasoned Oxygen (mmol/m³)",
    o2_375 = "Oxygen @ 375m (mmol/m³)",
    o2_375_anom = "De-seasoned Oxygen @ 375m (mmol/m³)",
    o2_1655 = "Oxygen @ 1655m (mmol/m³)",
    o2_1655_anom = "De-seasoned Oxygen @ 1655m (mmol/m³)",
    
    productivity_0 = "Net Primary Production (mg/m³/day C)",
    productivity_0_anom = "De-seasoned Net Primary Production (mg/m³/day C)",
    productivity_375 = "Net Primary Production (mg/m³/day C) @ 375m",
    productivity_375_anom = "De-seasoned Net Primary Production (mg/m³/day C) @ 375m"
  )
  
  if (var %in% names(labels)) return(labels[[var]])
  
  # ---- 2) Helpers for pattern labels ----
  depth_suffix_label <- function(d) {
    if (d == "0") return("Sea Surface")
    paste0("@ ", d, "m")
  }
  
  # ---- 3) Handle lag variables like chla_3mon, productivity_6mon, o2_12mon ----
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
  
  # ---- 4) Handle *_sd and *_mad WITH depth (temp_sd_0, o2_sd_1655, EKE_mad_375, etc.) ----
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
  
  # ---- 5) Handle *_sd and *_mad WITHOUT depth (fsle_sd, ssh_sd, mixed_layer_sd, fsle_orient_sd) ----
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
      base
    )
    
    stat_label <- if (stat == "sd") "SD" else "MAD"
    return(paste0(stat_label, " of ", base_label))
  }
  
  # ---- 6) Handle plain depth variables you didn't explicitly list (o2_1655, chla_0, etc.) ----
  if (grepl("_(0|375|1655)$", var)) {
    base <- sub("_(0|375|1655)$", "", var)
    depth <- sub("^.*_(0|375|1655)$", "\\1", var)
    
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
  
  # ---- 7) Anomaly suffix fallback ----
  if (grepl("_anom$", var)) {
    base <- sub("_anom$", "", var)
    return(paste0("De-seasoned ", nameVar(base)))
  }
  
  # ---- fallback ----
  var
}

# Generating visualizations for each site's final model
EI_plots <- visualizeGAM(EI_final,'EI')
CI_plots <- visualizeGAM(CI_final,'CI')

# ------------------ Step 7: GAM Model Checking -------------------
# Helper: quick ACF plot as ggplot
acf_df <- function(x, lag.max = 200) {
  a <- acf(x, plot = FALSE, lag.max = lag.max, na.action = na.pass)
  tibble(lag = as.numeric(a$lag), acf = as.numeric(a$acf))
}

check_gam <- function(mod,
                      model_name = "GAM",
                      date = NULL,            # optionally pass the bin_start vector for residual ACF ordering
                      lag.max = 200,
                      do_gam_check = TRUE,
                      do_concurvity = TRUE,
                      do_appraise = TRUE,
                      print_tables = TRUE) {
  
  cat("\n==============================\n")
  cat("Model:", model_name, "\n")
  cat("==============================\n")
  
  # --- 1) Basic mgcv checks (k-index, residual patterns, QQ, etc.) ---
  if (isTRUE(do_gam_check)) {
    cat("\n--- mgcv::gam.check() ---\n")
    # This prints the k-index table + 4 diagnostic panels
    suppressWarnings(gam.check(mod))
  }
  
  # --- 2) gratia appraise (nice ggplots) ---
  p_appraise <- NULL
  if (isTRUE(do_appraise)) {
    # appraise() returns a patchwork object of diagnostics
    # (residuals vs fitted, QQ, histogram, response vs fitted)
    p_appraise <- appraise(mod) +
      plot_annotation(title = paste0(model_name, " | gratia::appraise()"))
    print(p_appraise)
  }
  
  # --- 3) Residual autocorrelation (important for time series) ---
  # Order residuals by date if provided; otherwise use model order
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
    labs(title = paste0(model_name, " | Deviance residual ACF"),
         x = "Lag", y = "ACF") +
    theme_bw()
  print(p_acf)
  
  # --- 4) Concurvity (GAM version of multicollinearity) ---
  if (isTRUE(do_concurvity)) {
    
    cat("\n--- mgcv::concurvity(full = TRUE) ---\n")
    
    conc <- mgcv::concurvity(mod, full = TRUE)
    
    # concurvity sometimes returns a matrix, sometimes a list
    conc_est <- if (is.list(conc) && !is.null(conc$estimate)) {
      conc$estimate
    } else {
      conc
    }
    
    print(conc_est)
    
    conc_tbl <- as.data.frame(conc_est)
    conc_tbl$term <- rownames(conc_est)
    conc_tbl <- dplyr::relocate(conc_tbl, term)
    
    cat("\nSorted by worst concurvity:\n")
    
    if ("worst" %in% names(conc_tbl)) {
      conc_tbl %>% arrange(desc(worst)) %>% print()
    } else {
      print(conc_tbl)
    }
  }
  
  # --- 5) Optional: a clean residuals vs time plot (often super informative) ---
  if (!is.null(date)) {
    df_time <- tibble(date = as.Date(date), resid = r) %>% arrange(date)
    p_time <- ggplot(df_time, aes(date, resid)) +
      geom_hline(yintercept = 0) +
      geom_line() +
      labs(title = paste0(model_name, " | Deviance residuals vs time"),
           x = NULL, y = "Deviance residual") +
      theme_bw()
    print(p_time)
  }
  
  invisible(list(appraise_plot = p_appraise, acf_plot = p_acf))
}

# ---- Run checks for each site model ----
# Use the *inference* fits (final2) since those are what you’ll interpret/report
check_gam(EI_final,  model_name = "EI_final (Pm)",  date = EI_binned_deseasoned$bin_start,  lag.max = 200)
check_gam(CI_final,  model_name = "CI_final (Pm)",  date = CI_binned_deseasoned$bin_start,  lag.max = 200)
