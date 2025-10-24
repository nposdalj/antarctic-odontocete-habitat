library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)
library(itsadug)

# ------------- Step 0: Choose Species ----------------
# Modeling Gm for all sites, 40 km radius environmental data
species <- c('Oo') # options: BW29, BW37, Oo, Pm, Gm
# BW29 = Southern bottlenose whale, BW37 = Gray's and strap-toothed whales
# Oo = Killer whale, Pm = Sperm Whale, Gm = Long-finned pilot whale
# Note: not enough data to model BW58 at any site
sites <- c('EI','KGI','CI')

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
sp_specific <- allData  %>% subset(select=-c(Gm,BW37,BW29,BW58,Pm)) %>%
  subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                  temperature_0,salinity_0,EKE_0,temperature_16,salinity_16,EKE_16,
                  temperature_635,salinity_635,EKE_635, chla_0,o2_0,productivity_0,chla_16,
                  o2_16,productivity_16, chla_635,o2_635,productivity_635,
                  ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                  chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime,fsle_orient, chla_1mon,
                  chla_2mon, chla_3mon, chla_4mon, chla_5mon, chla_6mon, productivity_1mon,
                  productivity_2mon, productivity_3mon, productivity_4mon, 
                  productivity_5mon, productivity_6mon))

# Intermediate step to find ACF for seasonal model only, need to go back and adjust Trisha's ACF function
allData_EI <- allData %>% 
  filter(Site == "EI")
allData_EI$Binary <- ifelse(allData_EI$Gm > 0, 1, 0)
BlockMod_EI <- gam(Oo ~ s(ice_conc,k=4),
                   family = tw(link = "log", a = 1.1, b = 1.9), data = allData_EI, method = "REML")
ACF = acf(residuals(BlockMod_EI), lag.max = 1500) 
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval_EI = ACFidx[1]

allData_KGI <- allData %>% 
  filter(Site == "KGI")
allData_KGI$Binary <- ifelse(allData_KGI$Gm > 0, 1, 0)
BlockMod_KGI <- gam(Oo ~ s(ice_conc,k=4),
                    family = tw(link = "log", a = 1.1, b = 1.9), data = allData_KGI, method = "REML")
ACF = acf(residuals(BlockMod_KGI), lag.max = 1500) 
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval_KGI = ACFidx[1]

allData_CI <- allData %>% 
  filter(Site == "CI")
allData_CI$Binary <- ifelse(allData_CI$Gm > 0, 1, 0)
BlockMod_CI <- gam(Oo ~ s(ice_conc,k=4),
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

binByACF <- function(site, bin) {
  # Filtering by site and adding bin_start date
  sp_filtered <- sp_specific %>% filter(Site==site) %>%
    mutate(bin_start = floor_date(date, unit = paste(bin, 'days')))
  
  # Function for taking mean with correct syntax for a column
  mean_col <- function(col) expr(mean(!!sym(col), na.rm = TRUE))
  
  # Expression to summarize species presence
  # replace this:
  # species_expr <- set_names(list(expr(as.integer(any(!!sym(species) > 0)))), species)
  
  # with this: mean of the raw species values per bin
  species_expr <- set_names(
    list(expr(mean(!!sym(species), na.rm = TRUE))),
    species
  )
  
  # Taking the average of environnmental variables at surface and at 768 m depth
  summarize_cols <- list(julian_day = mean_col("julian_day"), AAO = mean_col("AAO"),
                         SSH = mean_col("SSH"), FSLE = mean_col("FSLE"), mixed_layer = mean_col("mixed_layer"),
                         ice_conc = mean_col("ice_conc"), ice_thickness = mean_col("ice_thickness"),
                         ice_diff = mean_col("ice_diff"),
                         temperature_0 = mean_col("temperature_0"), salinity_0 = mean_col("salinity_0"),
                         EKE_mad_0 = mean_col("EKE_mad_0"), chla_0 = mean_col('chla_0'),
                         o2_0 = mean_col('o2_0'), productivity_0 = mean_col('productivity_0'),
                         SSH_sd = mean_col('ssh_sd'), FSLE_sd = mean_col('fsle_sd'), 
                         mixed_layer_sd = mean_col('mixed_layer_sd'), temp_0_sd = mean_col('temp_sd_0'),
                         chla_0_sd = mean_col('chla_sd_0'), productivity_0_sd = mean_col('productivity_sd_0'),
                         EKE_0_mad = mean_col('EKE_mad_0'), o2_0_sd = mean_col('o2_sd_0'),
                         salinity_0_sd = mean_col('salinity_sd_0'),temperature_635 = mean_col('temperature_635'),
                         salinity_635 = mean_col('salinity_635'), EKE_635 = mean_col('EKE_635'),
                         o2_635 = mean_col('o2_635'), fsle_orient= mean_col('fsle_orient'), 
                         temperature_16 = mean_col('temperature_16'),
                         salinity_16 = mean_col('salinity_16'), EKE_16 = mean_col('EKE_16'),
                         o2_16 = mean_col('o2_16'))
  
  # Adding shallow dive depth for applicable variables
  # Customize this list to add standard deviations of variables at depth (if needed)
  summarize_cols[[paste0("temperature_", depths[2])]] <- mean_col(paste0("temperature_", depths[2]))
  summarize_cols[[paste0("salinity_", depths[2])]] <- mean_col(paste0("salinity_", depths[2]))
  summarize_cols[[paste0("EKE_", depths[2])]] <- mean_col(paste0("EKE_", depths[2]))
  summarize_cols[[paste0("chla_", depths[2])]] <- mean_col(paste0("chla_", depths[2]))
  summarize_cols[[paste0("o2_", depths[2])]] <- mean_col(paste0("o2_", depths[2]))
  summarize_cols[[paste0("productivity_", depths[2])]] <- mean_col(paste0("productivity_", depths[2]))
  # Adding deep dive depth for species with multiple depths
  if (species != "BW29") {
    summarize_cols[[paste0("temperature_", depths[3])]] <- mean_col(paste0("temperature_", depths[3]))
    summarize_cols[[paste0("salinity_", depths[3])]] <- mean_col(paste0("salinity_", depths[3]))
    summarize_cols[[paste0("EKE_", depths[3])]] <- mean_col(paste0("EKE_", depths[3]))
    summarize_cols[[paste0("chla_", depths[3])]] <- mean_col(paste0("chla_", depths[3]))
    summarize_cols[[paste0("o2_", depths[3])]] <- mean_col(paste0("o2_", depths[3]))
    summarize_cols[[paste0("productivity_", depths[3])]] <- mean_col(paste0("productivity_", depths[3]))
  }
  
  # Joining together dataframes and grouping to make final binned data
  all_summaries <- c(species_expr, summarize_cols)
  sp_binned <- sp_filtered %>%
    group_by(bin_start, Site) %>%
    summarise(!!!all_summaries, .groups = "drop")
  
  # Adding ice_regime column
  for (x in 1:nrow(sp_binned)) {
    if(sp_binned[x,'ice_diff'] == 0) {
      sp_binned[x,'ice_regime'] <- 'none'
    } else if(sp_binned[x,'ice_diff'] <= -0.01) {
      sp_binned[x,'ice_regime'] <- 'decreasing'
    } else if(sp_binned[x,'ice_diff'] >= 0.01) {
      sp_binned[x,'ice_regime'] <- 'increasing'
    } else
      sp_binned[x,'ice_regime'] <- 'stable'
  }
  
  return(sp_binned)
}
#EI_acf <- acfVal('EI')
EI_binned <- binByACF('EI',ACFval_EI)
#KGI_acf <- acfVal('KGI')
KGI_binned <- binByACF('KGI',ACFval_KGI)
#CI_acf <- acfVal('CI')
CI_binned <- binByACF('CI',ACFval_CI)

# ------------- Step 3: Plot Presence Timeseries --------------
binnedTimeseries <- function(data,site, bin) { # Function to create a timeseries plot
  # Making timeseries 
  ggplot(data = data, mapping = aes(x = bin_start, y = get(species))) + geom_col(width = 1, color = "slateblue") +
    scale_x_date(date_labels = "%b %Y")+
    labs(subtitle = name(site), y = NULL, x = NULL) + 
    theme(plot.subtitle = element_text(size = 9, face = "bold"), 
          plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), units = "line"))
  
}
binned_plot <- grid.arrange(binnedTimeseries(EI_binned,'EI',ACFval_EI),
                            binnedTimeseries(KGI_binned,'KGI',ACFval_KGI),
                            binnedTimeseries(CI_binned,'CI',ACFval_CI), nrow=3, 
                            top = paste('ACF Binned Species Presence for ', name(species), sep=''))


# -------------- Step 4: VIF for Correlation -------------------
# Not including any variables at depth for initial predictors.
# Also, not modeling AAO (varies on yearly timescales) and ice thickness
# Only including julian day for KGI, since it has almost 1 year of data

#VIF FUNCTION


# ELEPHANT ISLAND
EI_pred <- c('ice_conc', 'ice_diff', 'ice_thickness')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
#....several steps later
EI_pred <- c('ice_conc', 'ice_diff')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# ice_conc ice_diff 
# 1.051221 1.051221 

# KING GEORGE ISLAND
KGI_pred <- c('ice_conc', 'ice_diff', 'ice_thickness')
mod_formula <- reformulate(KGI_pred, response = species)
KGI_vif <- lm(mod_formula, data = KGI_binned)
vif(KGI_vif)
#...several steps later
KGI_pred <- c('ice_conc', 'ice_diff')
mod_formula <- reformulate(KGI_pred, response = species)
KGI_vif <- lm(mod_formula, data = KGI_binned)
vif(KGI_vif)
# ice_conc ice_diff 
# 1.004708 1.004708 

# CLARENCE ISLAND
CI_pred <- c('ice_conc', 'ice_diff', 'ice_thickness')
mod_formula <- reformulate(CI_pred, response = species)
CI_vif <- lm(mod_formula, data = CI_binned)
vif(CI_vif)
# ice_conc      ice_diff ice_thickness 
# 3.270406      1.005895      3.262777 

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

# -------------------- Step 5a: Elephant Island GAM ------------------------------
#Full model
EI_gam <- gam(Oo ~ s(ice_conc,k=3,sp=0.1) + s(ice_diff,k=4,sp=0.1),
              family = tw(link = "log", a = 1.1, b = 1.9), data = EI_binned, method = "REML")

# Remove ice_diff
EI_final <- gam(Oo ~ s(ice_conc,k=3,sp=0.1),
                family = tw(link = "log", a = 1.1, b = 1.9), data = EI_binned, method = "REML")


# -------------------- Step 5b: King George Island GAM -------------------------
# Full GAM
KGI_gam <- gam(Oo ~ s(ice_conc,k=4,sp=0.1) + s(ice_diff,k=4,sp=0.1),
               family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

# Remove ice_diff
KGI_final <- gam(Oo ~ s(ice_conc,k=3),
                 family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")


# --------------------- Step 5c: Clarence Island GAM ------------------------------------
# Final model
CI_gam <- gam(Oo ~ s(ice_conc,k=3,sp=0.1) + s(ice_diff,k=4,sp=0.1) + s(ice_thickness,k=4,sp=0.1),
              family = tw(link = "log", a = 1.1, b = 1.9),  data = CI_binned, method = "REML")

# Remove ice_diff
CI_gam <- gam(Oo ~ s(ice_conc,k=3,sp=0.1) + s(ice_thickness,k=3,sp=0.1),
              family = tw(link = "log", a = 1.1, b = 1.9),  data = CI_binned, method = "REML")

# Remove ice thickness
CI_final <- gam(Oo ~ s(ice_conc,k=3,sp=0.1),
                family = tw(link = "log", a = 1.1, b = 1.9),  data = CI_binned, method = "REML")


# ------------------ Step 6: Visualize GAMs -------------------
# Function to create a cleaner visualization of a GAM model
visualizeGAM <- function(gam, predictors, sp) {
  # Accessing and preparing data
  # extracting data to use for a ggplot of the GAM
  plot_info <- smooth_estimates(gam)
  
  # adding confidence interval
  plot_info <- plot_info %>% add_confint()
  
  # shifting by intercept and transforming to logistic (probability as y-axis)
  dont_shift <- names(plot_info) %in% c('.smooth','.type','.by','.se',predictors)
  plot_info <- plot_info %>%
    gratia:::shift_values(i = dont_shift, h = coef(gam)[1], FUN = '+') %>%
    transform_fun(fun = plogis)
  
  # Extracting deviance explained and p-values from the model
  summary <- summary(gam)
  deviance <- round(summary$dev.expl * 100, 3)
  p_values <- setNames(summary$s.pv, rownames(summary$s.table))
  
  # Creating list to store plots in 
  all_plots <- list()
  
  # Looping through all terms in GAM model
  for(p in predictors) {
    # Keeping only the current predictor
    current_plot <- filter(plot_info, .smooth == paste0('s(',p,')'))
    
    # Extracting p-value
    current_p_val <- p_values[[paste0("s(", p, ")")]]
    if(current_p_val == 0) {
      current_p_val <- 1 * 10^-6
    }
    current_plot$label <- paste0("p-value = ", round(current_p_val,6))
    
    # Setting position for p-value label
    current_plot$label_x <- max(current_plot[,p],na.rm=TRUE) - 
      (0.2*((max(current_plot[,p],na.rm=TRUE)) - (max(current_plot[,p],na.rm=TRUE))))
    current_plot$label_y <- 0.98
    
    # Setting limits for x axis
    x_vals <- pull(current_plot, p)
    x_lim <- range(x_vals, na.rm = TRUE)
    
    # Creating plot
    plot <- ggplot() +
      # Confidence interval ribbon plot
      geom_ribbon(data=current_plot,
                  aes(ymin = .lower_ci, ymax = .upper_ci, y = .estimate, x = .data[[p]]), alpha = 0.2) +
      
      # Line plot for smooth function
      geom_line(data=current_plot,aes(y = .estimate, x = .data[[p]]), linewidth = 1) +
      
      # Rug plot to show observed predictor values
      geom_rug(data = gam$model, aes(x = .data[[p]]), sides = "b", length = grid::unit(0.03, "npc")) +
      
      # Adding p-value label
      geom_label(data=current_plot,
                 aes(x = label_x, y = label_y, label = label), stat = "unique",
                 size = 3.5, alpha = 0.6, label.padding = unit(0.4, "lines"),hjust=1.2,vjust=1) +
      
      # Styling and cropping plot
      labs(y = "Partial effect", x = nameVar(p)) + theme_bw() + ylim(0,1) + xlim(x_lim) + 
      scale_x_continuous(expand = c(0, 0))
    
    # Adding plot object to list of plots
    all_plots[[length(all_plots)+1]] <- plot
  }
  
  # Setting rows and columns for final display
  if(length(predictors) == 1) {
    row <- 1
    col <- 1
  } else if(length(predictors) == 2) {
    row <- 1
    col <- 2
  } else if(length(predictors) == 3 || length(predictors) == 4) {
    row <- 2
    col <- 2
  } else if(length(predictors) == 5 || length(predictors) == 6) {
    row <- 2
    col <- 3
  } else {
    row <- 3
    col <- 3
  }
  
  # Aggregating all the plots into one figure
  final_plot <- wrap_plots(all_plots, nrow = row, ncol = col, guides = "collect") &
    plot_annotation(title = paste0("Killer Whale at ",
                                   sp," Presence (",deviance,"% Deviance Explained)"))
  
  print(final_plot)
  return(final_plot)
}

# Function to generate axis names from given variable names
nameVar <- function(var) {
  if(paste(var) == 'julian_day') {
    return("Julian Day")
  } else if(paste(var) == 'SSH') {
    return('Sea Surface Height (m)')
  } else if(paste(var) == 'AAO') {
    return('Antarctic Oscillation Index')
  } else if(paste(var) == 'FSLE,fsle_orient') {
    return('FSLE Magnitude')
  } else if(paste(var) == 'mixed_layer') {
    return('Mixed Layer Depth (m)')
  } else if(paste(var) == 'ice_conc') {
    return('Sea Ice Concentration')
  } else if(paste(var) == 'ice_thickness') {
    return('Sea Ice Thickness')
  } else if(paste(var) == 'ice_diff') {
    return('Difference in Sea Ice Concentration')
  } else if(paste(var) == 'temperature_0') {
    return('Sea Surface Temperature (°C)')
  } else if(paste(var) == 'salinity_0') {
    return('Sea Surface Salinity (psu)')
  } else if(paste(var) == 'salinity_16') {
    return('Sea Surface Salinity @ 16m (psu)')
  } else if(paste(var) == 'EKE_0') {
    return('Eddy Kinetic Energy (cm^2/s^2)')
  } else if(paste(var) == 'EKE_16') {
    return('Eddy Kinetic Energy @ 16m (cm^2/s^2)')
  } else if(paste(var) == 'EKE_635') {
    return('Eddy Kinetic Energy @ 635m (cm^2/s^2)')
  } else if(paste(var) == 'chla_0') {
    return('Chlorophyll (mg/m^3)')
  } else if(paste(var) == 'chla_16') {
    return('Chlorophyll @ 16m (mg/m^3)')
  } else if(paste(var) == 'o2_0') {
    return('Oxygen (mmol/m^3)')
  } else if(paste(var) == 'o2_635') {
    return('Oxygen @ 635m (mmol/m^3)')
  } else if(paste(var) == 'productivity_0') {
    return('Net Primary Production (mg/m^3/day carbon)')
  } else if(paste(var) == 'productivity_16') {
    return('Net Primary Production @ 16m (mg/m^3/day carbon)')
  } 
}

# Generating visualizations for each site's final model
KGI_pred <- c('ice_conc')
KGI_plots <- visualizeGAM(KGI_final, KGI_pred, 'KGI')

EI_pred <- c('ice_conc')
EI_plots <- visualizeGAM(EI_final, EI_pred, 'EI')

CI_pred <- c('ice_conc')
CI_plots <- visualizeGAM(CI_final, CI_pred, 'CI')
