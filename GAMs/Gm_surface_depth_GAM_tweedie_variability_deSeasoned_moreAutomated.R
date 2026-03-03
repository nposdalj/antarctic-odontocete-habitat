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
# Modeling Gm for all sites, 40 km radius environmental data
species <- c('Gm') # options: BW29, BW37, Oo, Pm, Gm
# BW29 = Southern bottlenose whale, BW37 = Gray's and strap-toothed whales
# Oo = Killer whale, Pm = Sperm Whale, Gm = Long-finned pilot whale
# Note: not enough data to model BW58 at any site
sites <- c('EI','KGI','CI')
ice_threshold <- 0.15

themes_keywords <- list(
  mesoscale = c("FSLE","fsle","SSH","ssh","EKE"),
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
depths <- c(0, 16, 635) # Gm depths

# species is assumed to be a string like "Gm"
# (same as in your code where you used get(species))
all_species <- c("BW29","BW37","BW58","Gm","Pm","Oo")
drop_species <- setdiff(all_species, species)

# keep columns that end with _0, _16, or _635
keep_depth_regex <- paste0("_(", paste(depths, collapse="|"), ")$")

sp_specific <- allData %>%
  select(-any_of(drop_species)) %>%
  select(
    date, Site, julian_day, all_of(species),
    matches(keep_depth_regex),    # <-- keeps temperature_0, temp_sd_0, o2_sd_635, etc.
    AAO:last_col()                # <-- keeps all other variables after AAO
  )

# Intermediate step to find ACF for seasonal model only, need to go back and adjust Trisha's ACF function
allData_EI <- allData %>% 
  filter(Site == "EI")
BlockMod_EI <- gam(Gm ~ s(FSLE,k=4) + s(mixed_layer,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) +
                     s(EKE_mad_0,k=4,sp=0.1) + s(salinity_16,k=4,sp=0.1) + s(o2_635,k=4,sp=0.1),
                   family = tw(link = "log", a = 1.1, b = 1.9), data = allData_EI, method = "REML")
ACF = acf(residuals(BlockMod_EI), lag.max = 1500) 
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval_EI = ACFidx[1]

allData_KGI <- allData %>% 
  filter(Site == "KGI")
BlockMod_KGI <- gam(Gm ~ s(FSLE,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) +  s(EKE_mad_0,k=4) + s(salinity_635,k=4) +
                      s(o2_635,k=4) + s(productivity_16, k = 4), 
                    family = tw(link = "log", a = 1.1, b = 1.9), data = allData_KGI, method = "REML")
ACF = acf(residuals(BlockMod_KGI), lag.max = 1500) 
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval_KGI = ACFidx[1]

allData_CI <- allData %>% 
  filter(Site == "CI")
BlockMod_CI <- gam(Gm ~ s(FSLE,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) +
                     s(EKE_mad_0,k=4) + s(productivity_0,k=4) +
                     s(salinity_635,k=4) + s(o2_635,k=4) + s(chla_16,k=4),
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

KGI_binned <- binByACF(
  data     = sp_specific,
  site     = "KGI",
  bin_days = ACFval_KGI,
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
binned_plot <- grid.arrange(binnedTimeseries(EI_binned,'EI',EI_acf), binnedTimeseries(KGI_binned,'KGI',KGI_acf), 
                            binnedTimeseries(CI_binned,'CI',CI_acf), nrow=3, 
                            top = paste('ACF Binned Species Presence for ', name(species), sep=''))

# Remove days with less than a certain percentage of sea ice concentration
# apply to one dataframe
EI_binned <- EI_binned %>%
  filter(ice_conc < ice_threshold)

# do the same for other sites
KGI_binned <- KGI_binned %>%
  filter(ice_conc < ice_threshold)

CI_binned <- CI_binned %>%
  filter(ice_conc < ice_threshold)

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
  
  # Shallow (16 m)
  "temperature_16",
  "salinity_16",
  "o2_16",
  "chla_16",
  "productivity_16",
  "EKE_16",
  
  # Deep (635 m)
  "temperature_635",
  "salinity_635",
  "o2_635",
  "chla_635",
  "productivity_635",
  "EKE_635"
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
KGI_binned_deseasoned <- make_monthly_anoms_binned(KGI_binned, vars_to_anom, date_col = "bin_start", standardize = FALSE)
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

#Build predictor variables
build_preds <- function(df,
                        species,
                        drop_depths = c(16),
                        keep_lags = c(3, 6),
                        chla_prod_depth_max = 50,
                        keep_sd_for = c("FSLE","SSH","mixed_layer","fsle_orient"),
                        keep_depth_sd_prefixes = c("temp_sd", "EKE_mad") # add "salinity_sd" if desired
) {
  nms <- names(df)
  
  # start: drop IDs, response, and ice/AAO/julian/bathy
  pred <- setdiff(
    nms,
    c("bin_start","Site", species,
      "julian_day","AAO",
      "ice_conc","ice_thickness","ice_diff","ice_regime",
      "bathymetry")
  )
  
  # drop anomalies / STL products
  pred <- pred[!grepl("(_anom$|_anomaly$|_stl$)", pred)]
  
  # ---- depth handling ----
  # anything ending in _<number>
  depth_vars <- pred[grepl("_\\d+$", pred)]
  depth_nums <- as.integer(sub(".*_(\\d+)$", "\\1", depth_vars))
  
  # drop depth=16 (and any other drop_depths)
  keep_depth <- !(depth_nums %in% drop_depths)
  
  # remove chla/productivity for depths > chla_prod_depth_max
  is_chla_prod <- grepl("^(chla|productivity)_", depth_vars)
  keep_chla_prod <- !(is_chla_prod & depth_nums > chla_prod_depth_max)
  
  depth_keep_vars <- depth_vars[keep_depth & keep_chla_prod]
  
  # ---- lags handling ----
  # keep only 3mon and 6mon lags (and drop other lags)
  lag_vars <- pred[grepl("_\\d+mon$", pred)]
  
  # keep only productivity_ and chla_ lags
  lag_vars <- lag_vars[grepl("^(productivity|chla)_", lag_vars)]
  
  lag_nums <- as.integer(sub(".*_(\\d+)mon$", "\\1", lag_vars))
  lag_keep_vars <- lag_vars[lag_nums %in% keep_lags]
  
  # ---- non-depth, non-lag "core" vars ----
  core_vars <- pred[!(pred %in% depth_vars) & !(pred %in% lag_vars)]
  
  # keep only selected core vars (plus their *_sd if listed)
  # core env vars you probably want at EI:
  core_allow <- c("FSLE","SSH","mixed_layer","fsle_orient")
  core_keep <- intersect(core_vars, core_allow)
  
  # SD vars: keep only those explicitly requested
  sd_keep <- character(0)
  
  # sd for named core variables (e.g., fsle_sd, ssh_sd, mixed_layer_sd, fsle_orient_sd)
  for (v in keep_sd_for) {
    # map FSLE -> fsle_sd; SSH -> ssh_sd; mixed_layer -> mixed_layer_sd; fsle_orient -> fsle_orient_sd
    if (v == "FSLE") cand <- "fsle_sd"
    else if (v == "SSH") cand <- "ssh_sd"
    else if (v == "mixed_layer") cand <- "mixed_layer_sd"
    else if (v == "fsle_orient") cand <- "fsle_orient_sd"
    else cand <- paste0(v, "_sd")
    if (cand %in% pred) sd_keep <- c(sd_keep, cand)
  }
  
  # depth SD/MAD families (e.g., temp_sd_0, EKE_mad_0, temp_sd_635, EKE_mad_635, etc.)
  for (pref in keep_depth_sd_prefixes) {
    sd_keep <- c(sd_keep, pred[grepl(paste0("^", pref, "_\\d+$"), pred)])
  }
  
  # final predictor set
  final <- sort(unique(c(core_keep, sd_keep, depth_keep_vars, lag_keep_vars)))
  
  final
}

assign_themes <- function(kept_predictors, themes_keywords) {
  
  theme_list <- list()
  
  for (theme in names(themes_keywords)) {
    
    keywords <- themes_keywords[[theme]]
    
    matched <- kept_predictors[
      sapply(kept_predictors, function(var) {
        any(sapply(keywords, function(k)
          grepl(k, var, ignore.case = TRUE)))
      })
    ]
    
    if (length(matched) > 0) {
      theme_list[[theme]] <- matched
    }
  }
  
  return(theme_list)
}

select_best_by_theme <- function(data, response, themes,
                                 family = tw(link="log", a=1.1, b=1.9)) {
  
  winners <- c()
  
  for (theme in names(themes)) {
    
    vars <- themes[[theme]]
    vars <- vars[vars %in% names(data)]  # keep only existing
    
    if (length(vars) == 0) next
    
    models <- list()
    aics <- c()
    devs <- c()
    
    for (v in vars) {
      mod <- gam(
        as.formula(paste(response, "~ s(", v, ", k=4)")),
        data = data,
        family = family,
        method = "REML"
      )
      
      models[[v]] <- mod
      aics[v] <- AIC(mod)
      devs[v] <- summary(mod)$dev.expl
    }
    
    best <- names(which.min(aics))
    
    cat("\nTheme:", theme,
        "\n  Winner:", best,
        "\n  AIC:", round(aics[best],2),
        "\n  DevExpl:", round(devs[best],3), "\n")
    
    winners <- c(winners, best)
  }
  
  return(winners)
}

#EI
EI_pred <- build_preds(EI_binned, species = species)

remove_zero_var <- function(df, vars) {
  vars[sapply(df[vars], function(x) sd(x, na.rm = TRUE) > 0)]
}
EI_pred <- remove_zero_var(EI_binned, EI_pred)

res_EI <- vif_stepwise_drop(
  data       = EI_binned,
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
  data = EI_binned,
  response = "Gm",
  themes = themes_EI
)

#KGI
KGI_pred <- build_preds(KGI_binned, species = species)

remove_zero_var <- function(df, vars) {
  vars[sapply(df[vars], function(x) sd(x, na.rm = TRUE) > 0)]
}
KGI_pred <- remove_zero_var(KGI_binned, KGI_pred)

res_KGI <- vif_stepwise_drop(
  data       = KGI_binned,
  response   = species,          # e.g., "Gm"
  predictors = KGI_pred,
  vif_thresh = 10,
  family     = gaussian(),       # like you were doing
  verbose    = TRUE
)

res_KGI$kept_predictors
car::vif(res_KGI$final_fit)

#CI
CI_pred <- build_preds(CI_binned, species = species)

remove_zero_var <- function(df, vars) {
  vars[sapply(df[vars], function(x) sd(x, na.rm = TRUE) > 0)]
}
CI_pred <- remove_zero_var(CI_binned, CI_pred)

res_CI <- vif_stepwise_drop(
  data       = CI_binned,
  response   = species,          # e.g., "Gm"
  predictors = CI_pred,
  vif_thresh = 10,
  family     = gaussian(),       # like you were doing
  verbose    = TRUE
)

res_CI$kept_predictors
car::vif(res_CI$final_fit)

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
    
    mod <- gam(form, data = data, family = family, method = "REML")
    
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
  data = EI_binned,
  response = "Gm",
  predictors = res_EI$kept_predictors
)


KGI_final <- auto_gam(
  data = KGI_binned,
  response = "Gm",
  predictors = res_KGI$kept_predictors
)

CI_final <- auto_gam(
  data = CI_binned,
  response = "Gm",
  predictors = res_CI$kept_predictors
)

# ------------------ Step 6: Visualize GAMs -------------------
# Function to create a cleaner visualization of a GAM model
visualizeGAM <- function(gam, sp) {
  
  # Automatically extract smooth variable names
  predictors <- gam$smooth |> 
    sapply(function(x) x$term) |> 
    unlist()
  
  plot_info <- smooth_estimates(gam) |> add_confint()
  
  dont_shift <- names(plot_info) %in% 
    c('.smooth','.type','.by','.se', predictors)
  
  plot_info <- plot_info |>
    gratia:::shift_values(i = dont_shift,
                          h = coef(gam)[1],
                          FUN = '+') |>
    transform_fun(fun = plogis)
  
  summ <- summary(gam)
  deviance <- round(summ$dev.expl * 100, 2)
  p_values <- setNames(summ$s.pv, rownames(summ$s.table))
  
  all_plots <- list()
  
  for(p in predictors) {
    
    current_plot <- dplyr::filter(
      plot_info,
      .smooth == paste0("s(", p, ")")
    )
    
    current_p_val <- p_values[[paste0("s(", p, ")")]]
    current_p_val <- max(current_p_val, 1e-6)
    
    current_plot$label <- paste0("p = ",
                                 signif(current_p_val, 3))
    
    x_vals <- current_plot[[p]]
    x_lim <- range(x_vals, na.rm = TRUE)
    
    plot <- ggplot(current_plot) +
      geom_ribbon(aes(x = .data[[p]],
                      ymin = .lower_ci,
                      ymax = .upper_ci),
                  alpha = 0.2) +
      geom_line(aes(x = .data[[p]],
                    y = .estimate),
                linewidth = 1) +
      geom_rug(data = gam$model,
               aes(x = .data[[p]]),
               sides = "b") +
      labs(y = "Partial effect",
           x = nameVar(p)) +
      theme_bw() +
      ylim(0,1) +
      xlim(x_lim)
    
    all_plots[[length(all_plots)+1]] <- plot
  }
  
  final_plot <- patchwork::wrap_plots(all_plots) +
    patchwork::plot_annotation(
      title = paste0("Long-finned Pilot Whale at ",
                     sp,
                     " (",
                     deviance,
                     "% deviance explained)")
    )
  
  print(final_plot)
  return(final_plot)
}

# Function to generate axis names from given variable names
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
    
    temperature_0 = "Sea Surface Temperature (°C)",
    temperature_0_anom = "De-seasoned Sea Surface Temperature (°C)",
    temperature_16 = "Temperature @ 16m (°C)",
    temperature_16_anom = "De-seasoned Temperature @ 16m (°C)",
    temperature_635 = "Temperature @ 635m (°C)",
    temperature_635_anom = "De-seasoned Temperature @ 635m (°C)",
    
    salinity_0 = "Sea Surface Salinity (psu)",
    salinity_0_anom = "De-seasoned Sea Surface Salinity (psu)",
    salinity_16 = "Salinity @ 16m (psu)",
    salinity_16_anom = "De-seasoned Salinity @ 16m (psu)",
    salinity_635 = "Salinity @ 635m (psu)",
    salinity_635_anom = "De-seasoned Salinity @ 635m (psu)",
    
    EKE_0 = "Eddy Kinetic Energy (0m)",
    EKE_0_anom = "De-seasoned Eddy Kinetic Energy (0m)",
    EKE_16 = "Eddy Kinetic Energy @ 16m",
    EKE_16_anom = "De-seasoned Eddy Kinetic Energy @ 16m",
    EKE_635 = "Eddy Kinetic Energy @ 635m",
    EKE_635_anom = "De-seasoned Eddy Kinetic Energy @ 635m",
    
    EKE_mad_0 = "Eddy Kinetic Energy Variability (0m)",
    EKE_mad_0_anom = "De-seasoned Eddy Kinetic Energy Variability (0m)",
    
    chla_0 = "Chlorophyll (mg/m³)",
    chla_0_anom = "De-seasoned Chlorophyll (mg/m³)",
    chla_16 = "Chlorophyll @ 16m (mg/m³)",
    chla_16_anom = "De-seasoned Chlorophyll @ 16m (mg/m³)",
    
    o2_0 = "Oxygen (mmol/m³)",
    o2_0_anom = "De-seasoned Oxygen (mmol/m³)",
    o2_16 = "Oxygen @ 16m (mmol/m³)",
    o2_16_anom = "De-seasoned Oxygen @ 16m (mmol/m³)",
    o2_635 = "Oxygen @ 635m (mmol/m³)",
    o2_635_anom = "De-seasoned Oxygen @ 635m (mmol/m³)",
    
    productivity_0 = "Net Primary Production (mg/m³/day C)",
    productivity_0_anom = "De-seasoned Net Primary Production (mg/m³/day C)",
    productivity_16 = "Net Primary Production @ 16m",
    productivity_16_anom = "De-seasoned Net Primary Production @ 16m"
  )
  
  if (var %in% names(labels)) {
    return(labels[[var]])
  } else {
    return(var)  # fallback if not defined
  }
}

# Generating visualizations for each site's final model
KGI_plots <- visualizeGAM(KGI_final,'KGI')
EI_plots <- visualizeGAM(EI_final,'EI')
CI_plots <- visualizeGAM(CI_final,'CI')
