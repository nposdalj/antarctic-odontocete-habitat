library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)

# ------------- Step 0: Choose Species ----------------
# Modeling Pm for all sites it is present (EI & CI), 40 km radius environmental data
species <- c('Oo') # options: BW29, BW37, Oo, Pm, Gm
# BW29 = Southern bottlenose whale, BW37 = Gray's and strap-toothed whales
# Oo = Killer whale, Pm = Sperm Whale, Gm = Long-finned pilot whale
# Note: not enough data to model BW58 at any site
sites <- c('CI') # Only modeling CI because not enough data for Oo at other sites

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
allData <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/allData_40km.csv")
allData <- allData %>% subset(select=-X)
allData$date <- as.Date(allData$date, "%Y-%m-%d")
# Filter by species relevant data
# Only adding standard deviations of surface variables, feel free to change that if needed
if(species =='BW29') {
  depths <- c(0, 768)
  sp_specific <- allData %>% subset(select=-c(BW37,BW58,Oo,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_768,salinity_768,EKE_768,
                    chla_0,o2_0,productivity_0,chla_768,o2_768,productivity_768,
                    ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime))
} else if(species =='BW37') {
  depths <- c(0, 67, 920) 
  sp_specific <- allData %>% subset(select=-c(BW29,BW58,Oo,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_67,salinity_67,EKE_67,
                    temperature_920,salinity_920,EKE_920, chla_0,o2_0,productivity_0,chla_67,
                    o2_67,productivity_67, chla_920,o2_920,productivity_920,
                    ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime))
} else if(species =='Oo') {
  depths <- c(0, 11, 455) 
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_11,salinity_11,EKE_11,
                    temperature_455,salinity_455,EKE_455, chla_0,o2_0,productivity_0,chla_11,
                    o2_11,productivity_11, chla_455,o2_455,productivity_455,
                    ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime))
} else if(species =='Pm') {
  depths <- c(0, 375, 1665)
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Oo,Gm)) %>%
    subset(select=c(date,julian_day,Site,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_375,salinity_375,EKE_375,
                    temperature_1665,salinity_1665,EKE_1665, chla_0,o2_0,productivity_0,chla_375,
                    o2_375,productivity_375, chla_1665,o2_1665,productivity_1665,
                    ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime))
} else if(species =='Gm') {
  depths <- c(0, 16, 635) 
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Oo,Pm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_16,salinity_16,EKE_16,
                    temperature_635,salinity_635,EKE_635, chla_0,o2_0,productivity_0,chla_16,
                    o2_16,productivity_16, chla_635,o2_635,productivity_635,
                    ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime))
} else
  print('Species code not valid. Check inputs.')


# ------------- Step 2: Average by ACF ------------
acf_table <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Autocorrelation/acf_table.csv")
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
  species_expr <- set_names(list(expr(as.integer(any(!!sym(species) == 1)))), species)
  
  # Taking the average of environnmental variables
  summarize_cols <- list(julian_day = mean_col("julian_day"), AAO = mean_col("AAO"),
                         SSH = mean_col("SSH"), FSLE = mean_col("FSLE"), mixed_layer = mean_col("mixed_layer"),
                         ice_conc = mean_col("ice_conc"), ice_thickness = mean_col("ice_thickness"),
                         ice_diff = mean_col("ice_diff"),
                         temperature_0 = mean_col("temperature_0"), salinity_0 = mean_col("salinity_0"),
                         EKE_0 = mean_col("EKE_0"), chla_0 = mean_col('chla_0'),
                         o2_0 = mean_col('o2_0'), productivity_0 = mean_col('productivity_0'),
                         SSH_sd = mean_col('ssh_sd'), FSLE_sd = mean_col('fsle_sd'), 
                         mixed_layer_sd = mean_col('mixed_layer_sd'), temp_0_sd = mean_col('temp_sd_0'),
                         chla_0_sd = mean_col('chla_sd_0'), productivity_0_sd = mean_col('productivity_sd_0'),
                         EKE_0_mad = mean_col('EKE_mad_0'), o2_0_sd = mean_col('o2_sd_0'),
                         salinity_0_sd = mean_col('salinity_sd_0'))
  
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
  
  return(sp_binned)
}
CI_acf <- acfVal('CI')
CI_binned <- binByACF('CI',CI_acf)

# ------------- Step 3: Plot Presence Timeseries --------------
binnedTimeseries <- function(data,site, bin) { # Function to create a timeseries plot
  # Making timeseries 
  ggplot(data = data, mapping = aes(x = bin_start, y = get(species))) + geom_col(width = 1, color = "slateblue") +
    scale_x_date(date_labels = "%b %Y")+
    labs(subtitle = name(site), y = NULL, x = NULL) + 
    theme(plot.subtitle = element_text(size = 9, face = "bold"), 
          plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), units = "line"))
  
}
binned_plot <- grid.arrange(
binnedTimeseries(CI_binned,'CI',CI_acf), nrow=1, 
                            top = paste('ACF Binned Species Presence for ', name(species), sep=''))
# -------------- Step 4: VIF for Correlation -------------------
# Not including any variables at depth for initial predictors.
# Also, not modeling AAO, julian day (vary on yearly timescales) and ice thickness
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0 
# 5.233589       8.574598       8.747562      13.393871       1.752501      18.192115      62.949706 
# EKE_0           o2_0         chla_0 productivity_0 
# 2.347280      13.991321      18.614977      27.334713 

# dropping temperature
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0 
# 3.425412       6.597550       7.244992       5.790618       1.645800      15.809096       2.329536 
# o2_0         chla_0 productivity_0 
# 3.515797      10.111880      17.523329 

# dropping primary production
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0",'o2_0','chla_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE         SSH mixed_layer    ice_conc    ice_diff  salinity_0       EKE_0        o2_0      chla_0 
# 2.508181    6.804206    4.981695    4.000060    1.503635   13.424566    2.201136    3.998256    2.357812

# dropping salinity
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', 
             "EKE_0",'o2_0','chla_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE         SSH mixed_layer    ice_conc    ice_diff       EKE_0        o2_0      chla_0 
# 2.490610    1.988522    2.476330    3.929869    1.502093    2.089286    2.785802    1.734605 

# final preditors for CI: FSLE, SSH, mixed layer depth, sea ice concentration,
#   difference in ice concentration, EKE, oxygen, chlorophyll


# -------------- Step 5: Build GAMs ------------------------
# Function to visualize GAMs on a probability scale with the proper confidence interval
# Run this for each iteration of the model to plot smooth terms
plotGam <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],seWithMean=TRUE))
}
# Run this if all plots should be in one figure
plotGam1 <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],seWithMean=TRUE,pages=1))
}

# CLARENCE ISLAND
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing at 0.1
# initial predictors: FSLE, SSH, mixed layer depth, sea ice concentration,
#   difference in sea ice concentration, eddy kinetic energy, oxygen concentration, chlorophyll

# weighing 1s at 8 to make ratio of 0s to 1s roughly 1:1 (reduce zero inflation)
CI_binned$weights <- ifelse(CI_binned$Oo == 1, 8, 1)

# FSLE
CI_gam <- gam(Oo ~ s(FSLE,k=4,sp=0.1), family=binomial, weights=weights, data=CI_binned)
# Formula:
#   Oo ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2747     0.1698  -1.618    0.106
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE) 1.984  2.361  16.31 0.000805 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0673   Deviance explained = 8.52%
# UBRE = 1.2533  Scale est. = 1         n = 105

# sea surface height
CI_gam <- gam(Oo ~ s(SSH,k=4,sp=0.1), family=binomial, weights=weights, data=CI_binned)
# Formula:
#   Oo ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.08334    0.14945  -0.558    0.577
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(SSH) 1.976  2.388  0.905   0.579
# 
# R-sq.(adj) =  -0.00586   Deviance explained = 1.08%
# UBRE = 1.4317  Scale est. = 1         n = 105

# mixed layer
CI_gam <- gam(Oo ~ s(mixed_layer,k=4,sp=0.1), family=binomial, weights=weights, data=CI_binned)
# Formula:
#   Oo ~ s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1155     0.1515  -0.763    0.446
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(mixed_layer) 1.904  2.289  4.323   0.133
# 
# R-sq.(adj) =  0.0177   Deviance explained = 2.96%
# UBRE = 1.3853  Scale est. = 1         n = 105

# sea ice concentration
CI_gam <- gam(Oo ~ s(ice_conc,k=4,sp=0.1), family=binomial, weights=weights, data=CI_binned)
# Formula:
#   Oo ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.6514     0.1972  -3.304 0.000954 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(ice_conc) 2.015  2.374  47.99  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.291   Deviance explained = 27.9%
# UBRE = 0.78784  Scale est. = 1         n = 105

# difference in ice concentration
CI_gam <- gam(Oo ~ s(ice_diff,k=4,sp=0.1), family=binomial, weights=weights, data=CI_binned)
# Formula:
#   Oo ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1599     0.1545  -1.035    0.301
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(ice_diff) 1.68  2.077  7.541  0.0277 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0531   Deviance explained = 5.48%
# UBRE = 1.3204  Scale est. = 1         n = 105

# eddy kinetic energy
CI_gam <- gam(Oo ~ s(EKE_0,k=4,sp=0.1), family=binomial, weights=weights, data=CI_binned)
# Formula:
#   Oo ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.09908    0.15155  -0.654    0.513
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.848   2.26  2.242   0.351
# 
# R-sq.(adj) =  -0.00288   Deviance explained = 1.46%
# UBRE = 1.4203  Scale est. = 1         n = 105

# oxygen concentration
CI_gam <- gam(Oo ~ s(o2_0,k=4,sp=0.1), family=binomial, weights=weights, data=CI_binned)
# Formula:
#   Oo ~ s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2207     0.1619  -1.363    0.173
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(o2_0) 2.106   2.48  10.81 0.00472 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0671   Deviance explained = 7.48%
# UBRE = 1.2805  Scale est. = 1         n = 105

# significant variables: FSLE, ice concentration, difference in ice concentration, oxygen
# building final model with these variables

# FSLE and ice concentration
CI_gam <- gam(Oo ~ s(FSLE,k=4,sp=0.1) + s(ice_conc,k=4,sp=0.1), 
              family=binomial, weights=weights, data=CI_binned)
# AIC: 182.897
# summary:
# Formula:
#   Oo ~ s(FSLE, k = 4, sp = 0.1) + s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.7025     0.1999  -3.514 0.000442 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(FSLE)     1.879  2.244  1.217   0.609    
# s(ice_conc) 1.899  2.253 36.152  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.311   Deviance explained = 31.2%
# UBRE = 0.74188  Scale est. = 1         n = 105
# ...................................................
# FSLE no longer significant, drop it
# concurvity and plots look good (though presence increases with ice concentration??)

# dropping FSLE and adding ice concentration difference
CI_gam <- gam(Oo ~ s(ice_diff,k=4,sp=0.1) + s(ice_conc,k=4,sp=0.1), 
              family=binomial, weights=weights, data=CI_binned)
# AIC: 187.75
# summary:
# Formula:
#   Oo ~ s(ice_diff, k = 4, sp = 0.1) + s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.6647     0.1972   -3.37 0.000751 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(ice_diff) 1.623  1.996  1.121   0.566    
# s(ice_conc) 1.985  2.332 43.467  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.293   Deviance explained = 29.3%
# UBRE = 0.78643  Scale est. = 1         n = 105
# .......................................................
# ice difference no longer significant, dropping it
# concurvity and plots look good (same weird presence & ice trend)

# dropping ice difference, adding oxygen concentration
CI_gam <- gam(Oo ~ s(o2_0,k=4,sp=0.1) + s(ice_conc,k=4,sp=0.1), 
              family=binomial, weights=weights, data=CI_binned)
# AIC: 190.5685
# summary:
# Formula:
#   Oo ~ s(o2_0, k = 4, sp = 0.1) + s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.6662     0.2000  -3.331 0.000866 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(o2_0)     1.835  2.200  0.644   0.802    
# s(ice_conc) 1.909  2.263 40.352  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.282   Deviance explained = 28.2%
# UBRE = 0.81494  Scale est. = 1         n = 105
# ..................................................
# dropping oxygen because it is not significant

# final variable: ice concentration (AIC = 187.7228)

# experimenting with different knots and smoothing
# trying more knots
CI_gam <- gam(Oo ~ s(ice_conc,k=8,sp=0.1), family=binomial, weights=weights, data=CI_binned) # AIC = 187.6161
# AIC went down a little bit, but model seems to be noiser, sticking with 4 knots
# trying different smoothing
CI_gam <- gam(Oo ~ s(ice_conc,k=4,sp=0.01), family=binomial, weights=weights, data=CI_binned) # AIC = 188.1105
CI_gam <- gam(Oo ~ s(ice_conc,k=4,sp=1), family=binomial, weights=weights, data=CI_binned) # AIC = 187.033
# trend becomes clearer at sp = 1, AIC goes down a little bit

# final model (AIC = 187,033, 27.6% deviance explained)
CI_final <- gam(Oo ~ s(ice_conc,k=4,sp=1), family=binomial, weights=weights, data=CI_binned)


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
    plot_annotation(title = paste0("Killer Whale at ",sp," Presence (",deviance,"% Deviance Explained)"))
  
  print(final_plot)
  return(final_plot)
}

# Function to generate axis names from given variable names
nameVar <- function(var) {
  if(paste(var) == 'julian_day') {
    return("Julian Day")
  } else if(paste(var) == 'SSH') {
    return('Sea Surface Height (m)')
  } else if(paste(var) == 'FSLE') {
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
  } else if(paste(var) == 'EKE_0') {
    return('Eddy Kinetic Energy (cm^2/s^2)')
  } else if(paste(var) == 'chla_0') {
    return('Chlorophyll (mg/m^3)')
  } else if(paste(var) == 'o2_0') {
    return('Oxygen (mmol/m^3)')
  } else if(paste(var) == 'productivity_0') {
    return('Net Primary Production (mg/m^3/day carbon)')
  } 
}

# Generating visualizations for each site's final model
CI_pred <- c('ice_conc')
CI_plots <- visualizeGAM(CI_final, CI_pred, 'CI')