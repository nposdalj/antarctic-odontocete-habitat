library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)

# ------------- Step 0: Choose Species ----------------
# Modeling Pm for all sites it is present (EI & CI), 40 km radius environmental data
species <- c('BW29') # options: BW29, BW37, Oo, Pm, Gm
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
EI_acf <- acfVal('EI')
EI_binned <- binByACF('EI',EI_acf)
KGI_acf <- acfVal('KGI')
KGI_binned <- binByACF('KGI',KGI_acf)
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
binned_plot <- grid.arrange(binnedTimeseries(EI_binned,'EI',EI_acf), binnedTimeseries(KGI_binned,'KGI',KGI_acf),
                            binnedTimeseries(CI_binned,'CI',CI_acf), nrow=3, 
                            top = paste('ACF Binned Species Presence for ', name(species), sep=''))

# -------------- Step 4: VIF for Correlation -------------------
# Not including any variables at depth for initial predictors.
# Also, not modeling AAO (varies on yearly timescales) and ice thickness
# Including julian day only for KGI because it covers almost a whole year

# ELEPHANT ISLAND
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0          EKE_0 
# 2.833739       2.856879       5.121568       3.794281       1.457163       4.353255      29.742875       1.491466 
# o2_0         chla_0 productivity_0 
# 22.806866      13.042628      19.709775 

# dropping surface temperature
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0           o2_0 
# 2.738288       2.291862       4.758202       2.844682       1.339495       4.339637       1.431586      15.438108 
# chla_0 productivity_0 
# 11.539813      18.579052 

# dropping primary production
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0",'o2_0','chla_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE         SSH mixed_layer    ice_conc    ice_diff  salinity_0       EKE_0        o2_0      chla_0 
# 1.806160    1.903622    4.532330    2.653441    1.295439    3.786823    1.297755    7.003150   10.557158 

# dropping chlorophyll
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0",'o2_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE         SSH mixed_layer    ice_conc    ice_diff  salinity_0       EKE_0        o2_0 
# 1.597221    1.630715    4.489612    2.011245    1.263594    2.298703    1.295622    5.432405 

# dropping oxygen
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0")
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE         SSH mixed_layer    ice_conc    ice_diff  salinity_0       EKE_0 
# 1.520709    1.540471    2.094414    1.987500    1.101917    1.509491    1.269623 

# final predictors for EI: FSLE, SSH, mixed layer depth, ice concentration, difference in ice concentration,
#   salinity, EKE

# KING GEORGE ISLAND
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0','julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0          EKE_0 
# 1.819185       1.838033       8.357019       2.224855       1.428314       3.324515      54.996723       1.140717 
# o2_0         chla_0 productivity_0     julian_day 
# 19.128410       8.033942      14.587727      10.206261 

# dropping temperature
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
              "EKE_0",'o2_0','chla_0','productivity_0','julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0           o2_0 
# 1.731630       1.684189       4.657725       2.219804       1.311916       3.160807       1.139838       8.149174 
# chla_0 productivity_0     julian_day 
# 7.650235       8.919327       9.734475 

# dropping julian day
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
              "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0           o2_0 
# 1.686486       1.676030       4.610218       2.199005       1.292921       3.022264       1.081825       3.255846 
# chla_0 productivity_0 
# 7.244615       8.210635 

# dropping primary production
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
              "EKE_0",'o2_0','chla_0')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE         SSH mixed_layer    ice_conc    ice_diff  salinity_0       EKE_0        o2_0      chla_0 
# 1.481281    1.584158    5.072583    1.945712    1.182638    2.546771    1.096340    2.103921    3.578545

# dropping mixed layer depth
KGI_pred <- c("FSLE", "SSH", "ice_conc", 'ice_diff', "salinity_0", 
              "EKE_0",'o2_0','chla_0')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE        SSH   ice_conc   ice_diff salinity_0      EKE_0       o2_0     chla_0 
# 1.529578   1.470737   1.580225   1.136723   1.911472   1.094967   1.850450   2.494099 

# final predictors for KGI: FSLE, SSH, ice concentration, difference in ice concentration,
#   salinity, EKE, oxygen, chlorophyll

# CLARENCE ISLAND
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
# CI vif model is not converging when run as a logistic regression, so running as a linear regression instead
CI_vif <- glm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0          EKE_0 
# 3.015100       7.437624       4.207211       3.052102       1.649232       7.196531      21.534652       1.210504 
# o2_0         chla_0 productivity_0 
# 7.602843       8.405939      30.503936 

# dropping primary production
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# FSLE           SSH   mixed_layer      ice_conc      ice_diff    salinity_0 temperature_0         EKE_0          o2_0 
# 3.010532      7.424873      3.972272      2.473190      1.648487      6.621887      4.798855      1.210336      3.787841 
# chla_0 
# 2.568309 

# dropping SSH
CI_pred <- c("FSLE", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# FSLE   mixed_layer      ice_conc      ice_diff    salinity_0 temperature_0         EKE_0          o2_0        chla_0 
# 2.989799      3.971194      2.469405      1.487423      2.761430      4.437775      1.149541      3.481934      2.568269 

# final predictors for CI: FSLE, mixed layer, sea ice concentration, difference in sea ice concentration,
#   salinity, temperature, EKE, oxygen, chlorophyll

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

# -------------------- Step 5a: Elephant Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing at 0.1
# NOT weighing for BW29 at EI since majority of bins have presence
# List of predictors are stored in EI_pred include FSLE, SSH, mixed layer depth, ice concentration,
#   difference in ice concentration, salinity, eddy kinetic energy

# starting with FSLE
EI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3829     0.2458   1.558    0.119
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(FSLE) 1.646  2.019   0.43   0.783
# 
# R-sq.(adj) =  -0.00457   Deviance explained = 1.42%
# UBRE = 0.40804  Scale est. = 1         n = 69

# sea surface height
EI_gam <- gam(BW29 ~ s(SSH,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW29 ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3907     0.2502   1.562    0.118
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(SSH) 1.732  2.115  3.099   0.241
# 
# R-sq.(adj) =  0.028   Deviance explained = 3.89%
# UBRE = 0.37725  Scale est. = 1         n = 69

# mixed layer depth
EI_gam <- gam(BW29 ~ s(mixed_layer,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW29 ~ s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   0.5213     0.2880    1.81   0.0703 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(mixed_layer) 1.35  1.619  8.041  0.0111 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.132   Deviance explained = 13.3%
# UBRE = 0.23842  Scale est. = 1         n = 69

# sea ice concentration
EI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.4061     0.2523    1.61    0.107
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.414  1.702  1.688   0.346
# 
# R-sq.(adj) =  0.00869   Deviance explained = 2.41%
# UBRE = 0.38795  Scale est. = 1         n = 69

# difference in ice concentration
EI_gam <- gam(BW29 ~ s(ice_diff,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW29 ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3868     0.2464    1.57    0.116
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_diff) 1.469  1.779  0.787   0.573
# 
# R-sq.(adj) =  0.00815   Deviance explained = 2.41%
# UBRE = 0.38965  Scale est. = 1         n = 69

# salinity
EI_gam <- gam(BW29 ~ s(salinity_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW29 ~ s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3829     0.2469   1.551    0.121
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(salinity_0) 1.606  1.978  1.316   0.548
# 
# R-sq.(adj) =  0.00286   Deviance explained = 1.96%
# UBRE = 0.39972  Scale est. = 1         n = 69

# eddy kinetic energy
EI_gam <- gam(BW29 ~ s(EKE_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW29 ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3842     0.2463    1.56    0.119
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.632  1.991  0.636   0.711
# 
# R-sq.(adj) =  -0.00631   Deviance explained = 1.29%
# UBRE = 0.40943  Scale est. = 1         n = 69

# only significant variable was mixed layer depth, trying different # of basis functions and smoothing parameters
EI_gam <- gam(BW29 ~ s(mixed_layer,k=4,sp=0.1), family=binomial, data=EI_binned)
# AIC: 85.45081

# trying more knots
EI_gam <- gam(BW29 ~ s(mixed_layer,k=8,sp=0.1), family=binomial, data=EI_binned)
# AIC: 86.02437
# AIC lower, sticking with 4 knots

# trying different smoothing
EI_gam <- gam(BW29 ~ s(mixed_layer,k=4,sp=0.01), family=binomial, data=EI_binned)
# AIC: 85.80641
# trying higher smoothing
EI_gam <- gam(BW29 ~ s(mixed_layer,k=4,sp=1), family=binomial, data=EI_binned)
# AIC: 85.08565
# AIC went down, error bars smaller, deviance explained decreased by only 0.2%
# sticking with higher smoothing
# trying even higher smoothing
EI_gam <- gam(BW29 ~ s(mixed_layer,k=4,sp=10), family=binomial, data=EI_binned)
# AIC: 85.02764
# marginal difference in confidence interval, deviance explained does not change
# marginal difference in AIC
# preferring sp of 1 to not oversmooth trends

# final model for EI
EI_final <- gam(BW29 ~ s(mixed_layer,k=4,sp=1), family=binomial, data=EI_binned)

# -------------------- Step 5b: King George Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in EI_pred include FSLE, SSH, ice concentration, 
#   difference in ice concentration, salinity, eddy kinetic energy, oxygen concentration, chlorophyll

# weighing 1s by 6 to make the ratio of 1s to 0s near 1:1
KGI_binned$weights <- ifelse(KGI_binned$BW29 == 1, 6, 1)

# starting with FSLE
KGI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1), weights=weights, family=binomial, data=KGI_binned)
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.2905     0.1319  -2.202   0.0276 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE) 2.143  2.535  26.76 5.67e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0823   Deviance explained = 8.55%
# UBRE = 1.1401  Scale est. = 1         n = 181

# sea surface height
KGI_gam <- gam(BW29 ~ s(SSH,k=4,sp=0.1), weights=weights, family=binomial, data=KGI_binned)
# Formula:
#   BW29 ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.0505     0.2669  -3.936 8.29e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(SSH) 1.529  1.864  51.68  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.222   Deviance explained = 22.8%
# UBRE = 0.80584  Scale est. = 1         n = 181

# sea ice concentration
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1), weights=weights, family=binomial, data=KGI_binned)
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -2.5332     0.9483  -2.671  0.00755 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc) 1.025  1.049  8.311  0.0027 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.203   Deviance explained = 22.7%
# UBRE = 0.80132  Scale est. = 1         n = 181

# difference in sea ice concentration
KGI_gam <- gam(BW29 ~ s(ice_diff,k=4,sp=0.1), weights=weights, family=binomial, data=KGI_binned)
# Formula:
#   BW29 ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1833     0.1219  -1.504    0.133
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_diff) 1.73  2.112  14.25 0.00239 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0521   Deviance explained = 5.29%
# UBRE = 1.2108  Scale est. = 1         n = 181

# salinity
KGI_gam <- gam(BW29 ~ s(salinity_0,k=4,sp=0.1), weights=weights, family=binomial, data=KGI_binned)
# Formula:
#   BW29 ~ s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.5588     0.3488  -4.469 7.86e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(salinity_0) 1.788  2.136  47.85  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.311   Deviance explained = 32.3%
# UBRE = 0.58941  Scale est. = 1         n = 181

# eddy kinetic energy
KGI_gam <- gam(BW29 ~ s(EKE_0,k=4,sp=0.1), weights=weights, family=binomial, data=KGI_binned)
# Formula:
#   BW29 ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1655     0.1248  -1.326    0.185
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(EKE_0) 1.404  1.681  3.752  0.0647 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0187   Deviance explained =  2.6%
# UBRE = 1.2691  Scale est. = 1         n = 181

# oxygen concentration
KGI_gam <- gam(BW29 ~ s(o2_0,k=4,sp=0.1), weights=weights, family=binomial, data=KGI_binned)
# Formula:
#   BW29 ~ s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.4329     0.1512  -2.863   0.0042 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(o2_0) 1.834  2.227   28.2 1.74e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.124   Deviance explained = 12.6%
# UBRE =  1.043  Scale est. = 1         n = 181

# chlorophyll
KGI_gam <- gam(BW29 ~ s(chla_0,k=4,sp=0.1), weights=weights, family=binomial, data=KGI_binned)
# Formula:
#   BW29 ~ s(chla_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.6381     0.1880  -3.395 0.000687 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(chla_0) 1.568  1.871  41.37 3.48e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.22   Deviance explained = 19.6%
# UBRE = 0.88042  Scale est. = 1         n = 181

# significant predictors for KGI: FSLE, SSH, ice concentration, difference in ice concentration, salinity
#   eddy kinetic energy, oxygen concentration, chlorophyll
# Building final GAM by adding these variables one at a time (and correcting for significance/concurvity)

# starting with ice concentration and FSLE
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 317.2929
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -2.791      1.044  -2.674   0.0075 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc) 1.016  1.032  7.523  0.0044 **
#   s(FSLE)     2.033  2.420  9.895  0.0105 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.232   Deviance explained = 25.8%
# UBRE =  0.753  Scale est. = 1         n = 181
# ..................................................
# concurvity looks good
# FSLE plot has large error bars (if standard error of intercept is included), very little impact on probability
# might drop later on

# adding sea surface height
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) + s(SSH,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 290.2585
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(FSLE, k = 4, sp = 0.1) + 
#   s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -3.799      1.203  -3.159  0.00158 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc) 1.006  1.012  6.289 0.010014 *  
#   s(FSLE)     1.953  2.334 13.219 0.001591 ** 
#   s(SSH)      1.356  1.612 25.900 0.000999 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.312   Deviance explained = 32.9%
# UBRE = 0.60364  Scale est. = 1         n = 181
# ....................................................
# concurvity looks good
# plots for FSLE and SSH have large error bars (if error of intercept included)
#   also have low impact on final probability, might drop later on

# adding in ice difference
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) +
                 s(ice_diff,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 288.796
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(FSLE, k = 4, sp = 0.1) + 
#   s(SSH, k = 4, sp = 0.1) + s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -4.039      1.262  -3.202  0.00137 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(ice_conc) 1.004  1.009  6.620 0.00854 ** 
#   s(FSLE)     1.945  2.325 13.279 0.00151 ** 
#   s(SSH)      1.360  1.618 26.823 0.00082 ***
#   s(ice_diff) 1.145  1.271  3.183 0.14490    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.317   Deviance explained = 33.8%
# .....................................................
# concurvity looks good
# dropping difference in ice concentration because it is not significant

# removing ice concentration difference, adding in salinity
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) +
                 s(salinity_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 275.6125
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(FSLE, k = 4, sp = 0.1) + 
#   s(SSH, k = 4, sp = 0.1) + s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -2.543      1.084  -2.347   0.0189 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.048  1.093  0.568 0.468175    
# s(FSLE)       1.899  2.280 14.689 0.000852 ***
#   s(SSH)        1.357  1.614  0.834 0.573759    
# s(salinity_0) 1.539  1.809 16.583 0.002123 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.358   Deviance explained = 37.1%
# UBRE = 0.52272  Scale est. = 1         n = 181
# ...........................................................
# salinity and ice concentration have problematic concurvity
# dropping FSLE first because it is no longer significant

# removing FSLE
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) +s(SSH,k=4,sp=0.1) +
                 s(salinity_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 291.9322
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(SSH, k = 4, sp = 0.1) + 
#   s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -2.2091     0.9176  -2.408   0.0161 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc)   1.075  1.144  0.694 0.43738   
# s(SSH)        1.408  1.696  0.480 0.72751   
# s(salinity_0) 1.621  1.921 14.651 0.00158 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =    0.3   Deviance explained = 32.4%
# UBRE = 0.61289  Scale est. = 1         n = 181
# .........................................................
# salinity concurvity still problematic with ice concentration
# removing salinity, even though it is significant right now

# removing salinity
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) +s(SSH,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 305.6803
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -3.122      1.051  -2.971  0.00297 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc) 1.013  1.026  6.109 0.01030 * 
#   s(SSH)      1.391  1.663 22.442 0.00207 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.259   Deviance explained = 28.3%
# UBRE = 0.68884  Scale est. = 1         n = 181
# ................................................................
# concurvity looks good
# plots still a little iffy - SSH error bars big

# adding oxygen concentration
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) +s(SSH,k=4,sp=0.1) +
                 s(o2_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 276.4975
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(SSH, k = 4, sp = 0.1) + 
#   s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -2.0176     0.6437  -3.134  0.00172 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc) 1.126  1.238  7.963  0.00499 ** 
#   s(SSH)      1.324  1.567  1.539  0.46401    
# s(o2_0)     1.556  1.876 27.630 1.51e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.346   Deviance explained = 36.1%
# UBRE = 0.52761  Scale est. = 1         n = 181
# .....................................................
# SSH is no longer significant, removing it
# concurvity is okay

# removing sea surface height
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) +s(o2_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 274.9579
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.8594     0.5527  -3.364 0.000768 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc) 1.160  1.296  12.02 0.000707 ***
#   s(o2_0)     1.622  1.974  32.91  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.349   Deviance explained = 35.8%
# UBRE = 0.5191  Scale est. = 1         n = 181
# ...............................................................
# concurvity looks good, plots look good, AIC went down

# adding chlorophyll
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) +s(o2_0,k=4,sp=0.1) +
                 s(chla_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 243.2747
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(o2_0, k = 4, sp = 0.1) + 
#   s(chla_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -4.266      1.243  -3.431 0.000601 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc) 1.016  1.031  10.28 0.001017 ** 
#   s(o2_0)     1.271  1.470  12.80 0.000647 ***
#   s(chla_0)   1.061  1.114  17.25 3.21e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.444   Deviance explained = 43.7%
# UBRE = 0.34406  Scale est. = 1         n = 181
# ..........................................................
# even though AIC went down and all variables are significant, chlorophyll & oxygen have high concurvity
# model without chlorophyll: AIC = 274.9579, deviance explained = 35.8%
# trying model without oxygen, will choose better of the two

# removing oxygen
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(chla_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 260.9019
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(chla_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -5.183      1.604   -3.23  0.00124 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc) 1.002  1.004  10.37 0.001132 ** 
#   s(chla_0)   1.280  1.490  15.29 0.000681 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.396   Deviance explained =   39%
# UBRE = 0.44145  Scale est. = 1         n = 181
# ............................................................
# large error bars because of model coefficient
# while chlorophyll has smaller AIC and more deviance explained, its effect size is very small based on graph
#   effect becomes negligible once model coefficient bars are included (which is not true for oxygen)
# so, keeping oxygen in final model and dropping chlorophyll

# final variables: oxygen and ice concentration
# now, experimenting with different k values and smoothing parameters

# trying different knots, starting with ice concentration
KGI_gam <- gam(BW29 ~ s(ice_conc,k=8,sp=0.1) + s(o2_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 271.5074
# even though AIC went down, error bars grow too much to justify it
# trying different knots for oxygen
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(o2_0,k=8,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 271.6329
# keeping knots at 4 for both variables, more knots seem to be fitting to noise

# trying different sp for ice concentration
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.01) + s(o2_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 273.0405
# error bars become massive, not using
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=1) + s(o2_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 275.34 
# AIC went up, not using

# trying different sp for oxygen
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(o2_0,k=4,sp=0.01),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 272.345
# AIC slightly went down, but error bars too large
# trying higher sp
KGI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(o2_0,k=4,sp=1),
               weights=weights, family=binomial, data=KGI_binned)
# AIC: 281.5063
# AIC went up, keeping smoothing at 0.1 for oxygen

# final model for KGI: ice concentration and oxygen
# AIC is 274.9599, 35.8% explained
KGI_final <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1),
               weights=weights, family=binomial, data=KGI_binned)

# -------------------- Step 5c: Clarence Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in EI_pred include FSLE, mixed layer depth, ice concentration,
#   difference in ice concentration, salinity, temperature, eddy kinetic energy, oxygen concentration, chlorophyll

# weighing 1s by 2 to make the ratio of 1s to 0s near 1:1
CI_binned$weights <- ifelse(CI_binned$BW29 == 1, 2, 1)

# starting with FSLE
CI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.02311    0.28313  -0.082    0.935
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(FSLE) 1.546  1.889  11.65  0.0016 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.244   Deviance explained = 21.7%
# UBRE = 0.58242  Scale est. = 1         n = 51

# mixed layer depth
CI_gam <- gam(BW29 ~ s(mixed_layer,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW29 ~ s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1293     0.2499   0.517    0.605
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(mixed_layer) 1.727  2.109  4.862  0.0927 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0429   Deviance explained = 5.84%
# UBRE = 0.88909  Scale est. = 1         n = 51

# sea ice concentration
CI_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -1.2378     0.7079  -1.749   0.0804 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc) 1.02   1.04  8.253 0.00302 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.356   Deviance explained = 32.1%
# UBRE = 0.36444  Scale est. = 1         n = 51

# difference in sea ice concentration
CI_gam <- gam(BW29 ~ s(ice_diff,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW29 ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.05615    0.26027   0.216    0.829
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_diff) 1.407  1.687  6.095   0.113
# 
# R-sq.(adj) =  0.0883   Deviance explained = 9.12%
# UBRE = 0.81455  Scale est. = 1         n = 51

# salinity
CI_gam <- gam(BW29 ~ s(salinity_0,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW29 ~ s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.08268    0.25154   0.329    0.742
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(salinity_0) 1.713  2.101  8.468  0.0154 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.146   Deviance explained = 13.3%
# UBRE = 0.74678  Scale est. = 1         n = 51

# temperature
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW29 ~ s(temperature_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1321     0.3145   -0.42    0.674
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(temperature_0) 1.572  1.905  20.92 2.01e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.472   Deviance explained = 41.8%
# UBRE = 0.20274  Scale est. = 1         n = 51

# eddy kinetic energy
CI_gam <- gam(BW29 ~ s(EKE_0,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW29 ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1577     0.2483   0.635    0.525
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.511  1.848  4.766   0.133
# 
# R-sq.(adj) =  0.039   Deviance explained = 5.88%
# UBRE = 0.87982  Scale est. = 1         n = 51

# oxygen concentration
CI_gam <- gam(BW29 ~ s(o2_0,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW29 ~ s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.04717    0.25958   0.182    0.856
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(o2_0) 1.798  2.155  14.52 0.00207 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.237   Deviance explained = 20.3%
# UBRE = 0.61837  Scale est. = 1         n = 51

# chlorophyll
CI_gam <- gam(BW29 ~ s(chla_0,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW29 ~ s(chla_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1418     0.2967  -0.478    0.633
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(chla_0) 1.679  2.016  21.56 1.99e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.421   Deviance explained = 35.2%
# UBRE = 0.33128  Scale est. = 1         n = 51

# final significant variables: FSLE, ice concentration, salinity, temperature, oxygen concentration, chlorophyll
# now building final model with these variables

# FSLE and ice concentration
CI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1) + s(ice_conc,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)
# AIC: 317.2929
# summary:
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1) + s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -2.791      1.044  -2.674   0.0075 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(FSLE)     2.033  2.420  9.895  0.0105 * 
#   s(ice_conc) 1.016  1.032  7.523  0.0044 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.232   Deviance explained = 25.8%
# UBRE =  0.753  Scale est. = 1         n = 181
# .........................................................
# FSLE error bars are large, even though it is significant consider removing in future
# concurvity looks good

# adding salinity
CI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1) + s(ice_conc,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)
# AIC: 273.7946
# summary:
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1) + s(ice_conc, k = 4, sp = 0.1) + 
#   s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -2.2967     0.9077   -2.53   0.0114 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE)       1.924  2.318 14.496  0.00102 ** 
#   s(ice_conc)   1.080  1.153  0.524  0.51502    
# s(salinity_0) 1.589  1.873 35.363 7.85e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.363   Deviance explained =   37%
# UBRE = 0.51268  Scale est. = 1         n = 181
# ......................................................
# concurvity strong between salinity and ice concentration
# removing ice concentration

# removing ice concentration
CI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1) + s(salinity_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)
# AIC: 272.0136
# summary:
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1) + s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7792     0.3796  -4.686 2.78e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE)       1.930  2.327  14.91 0.000886 ***
#   s(salinity_0) 1.697  2.021  46.31  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.366   Deviance explained = 36.9%
# UBRE = 0.50284  Scale est. = 1         n = 181
# ....................................................
# plots look good, concurvity is good

# adding temperature
CI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1) + s(salinity_0,k=4,sp=0.1) +
                s(temperature_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)
# AIC: 263.4327
# summary:
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1) + s(salinity_0, k = 4, sp = 0.1) + 
#   s(temperature_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.1006     0.4453  -4.717  2.4e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE)          1.879  2.261  18.81 0.000157 ***
#   s(salinity_0)    1.525  1.802  46.28 4.58e-07 ***
#   s(temperature_0) 1.842  2.162  11.10 0.006597 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.396   Deviance explained = 39.8%
# UBRE = 0.45543  Scale est. = 1         n = 181
# ..................................................
# salinity and temperature have high (ish) concurvity
# plots look okay (trend for temperature is not super present)
# will remove temperature later if it drops in significance

# adding oxygen
CI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1) + s(salinity_0,k=4,sp=0.1) +
                s(temperature_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)
# AIC: 226.8324
# summary:
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1) + s(salinity_0, k = 4, sp = 0.1) + 
#   s(temperature_0, k = 4, sp = 0.1) + s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.4109     0.4352   -5.54 3.03e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE)          1.723  2.089  18.88 0.000105 ***
#   s(salinity_0)    1.574  1.857  41.62  < 2e-16 ***
#   s(temperature_0) 1.671  1.931  28.63 1.18e-05 ***
#   s(o2_0)          1.258  1.440  21.63 2.17e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.49   Deviance explained =   49%
# UBRE = 0.25322  Scale est. = 1         n = 181
# ..........................................................
# plots loook really good, improved temperature and FSLE trend
# potentially problematic concurvity between salinity & temp and salinity & oxygen
#   will address later

# adding chlorophyll
CI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1) + s(salinity_0,k=4,sp=0.1) +
                s(temperature_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1) + s(chla_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)
# AIC: 228.5841
# summary:
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1) + s(salinity_0, k = 4, sp = 0.1) + 
#   s(temperature_0, k = 4, sp = 0.1) + s(o2_0, k = 4, sp = 0.1) + 
#   s(chla_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.4164     0.5361  -4.507 6.57e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE)          1.720  2.085 18.710 0.000112 ***
#   s(salinity_0)    1.569  1.845 29.041 3.46e-06 ***
#   s(temperature_0) 1.604  1.864 18.300 0.002254 ** 
#   s(o2_0)          1.247  1.417 16.450 0.000213 ***
#   s(chla_0)        1.074  1.133  0.022 0.931852    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.487   Deviance explained = 49.1%
# UBRE = 0.2629  Scale est. = 1         n = 181
# ......................................................
# chlorophyll not significant, will drop it
# will also drop salinity, see how it impacts concurvity 

# dropping salinity and chlorophyll
CI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1) +
                s(temperature_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)
# AIC: 289.993
# summary:
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1) + s(temperature_0, k = 4, sp = 0.1) + 
#   s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.0235     0.1953  -5.241  1.6e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE)          1.957  2.342  17.18  0.00106 ** 
#   s(temperature_0) 2.017  2.385  26.37 2.93e-06 ***
#   s(o2_0)          1.644  1.993  20.56 6.71e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.341   Deviance explained = 33.6%
# UBRE = 0.60217  Scale est. = 1         n = 181
# ........................................................
# concurvity is good, plots look great
# replacing temperature with salinity to see if it's better

# replacing temperature with salinity
CI_gam <- gam(BW29 ~ s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)
# AIC: 263.6876
# summary:
# Formula:
#   BW29 ~ s(FSLE, k = 4, sp = 0.1) + s(salinity_0, k = 4, sp = 0.1) + 
#   s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.6193     0.3355  -4.827 1.39e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE)       1.873  2.259   9.97  0.00694 ** 
#   s(salinity_0) 1.739  2.052  30.66 1.36e-06 ***
#   s(o2_0)       1.520  1.848  10.51  0.01408 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.39   Deviance explained = 39.7%
# UBRE = 0.45684  Scale est. = 1         n = 181
# ........................................................
# while AIC is lower/explained is higher, plot trends have become weaker
# results comparable to keeping salinity and temperature, removing oxygen (line 1218)
# going to try keeping temperature, salinity, and oxygen while removing FSLE 
#   since its trend is smallest in magnitude

# replacing FSLE with temperature
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)
# AIC: 245.2002
# summary:
# Formula:
#   BW29 ~ s(temperature_0, k = 4, sp = 0.1) + s(salinity_0, k = 4, 
#                                                sp = 0.1) + s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.8390     0.3471  -5.298 1.17e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(temperature_0) 1.784  2.027  21.05 3.75e-05 ***
#   s(salinity_0)    1.694  2.019  43.74  < 2e-16 ***
#   s(o2_0)          1.305  1.521  21.05 2.24e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.44   Deviance explained = 43.9%
# UBRE = 0.3547  Scale est. = 1         n = 181
# .................................................
# plots look good, AIC looks good
# however, concurvity got worse

#  next step: determining which variables to keep/discard based on concurvity
# table of AIC and deviance explained if a variable is discarded (leaving FSLE + other two vars)
# variable excluded:         AIC:           deviance explained:
# oxygen                     263            39.8
# FSLE                       245            43.9
# temperature                263            39.7
# salinity                   289            33.6
# none excluded              226            49
#
# cutoff for concurvity is typically ~0.8, only pairwise concurvity above is temp & salinity (0.84)
# high concurvity typically means large confidence intervals, however confidence intervals
#   for model with highest concurvity are most reasonable
# "none excluded" model with highest concurvity also has lowest AIC and most deviance explained
# concurvity is not too high, AIC is low, and deviance explained is high ->
#   choosing to not exclude, final model will have FSLE, salinity, temperature, and oxygen

# now, experimenting with knots and smoothing for relevant variables
# starting with knots
# temperature
CI_gam <- gam(BW29 ~ s(temperature_0,k=8,sp=0.1) + s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 226.4493
# keeping at 4 for temp
# FSLE
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=8,sp=0.1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 228.1161
# keeping at 4 for FSLE
# salinity
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=8,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 218.6601
# even though AIC decreased, plot trend didn't significantly change (now more fit to noise?)
# keeping at 4 for salinity
# oxygen
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=8,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 225.397
# AIC went down but confidence interval went up, keeping oxygen k at 4

# experimenting with smoothing parameter
# temperature
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.01) + s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 226.9692
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=1) + s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 228.2801
# AIC didn't decrease, not changing sp
# FSLE
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=4,sp=0.01) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 227.8694
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=4,sp=1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 227.5535
# AIC went up, not changing sp for FSLE
# salinity
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.01) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 221.9909
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 233.1805
# not changing smoothing for salinity to 0.01 because even though AIC went down
#   too much of an increase in confidence interval for FSLE
# oxygen
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.01), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 226.1212
CI_gam <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=4,sp=1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=1), 
              weights=weights, family=binomial, data=KGI_binned) # AIC: 227.9339
# marginal decrease in AIC, large increase in confidence interval for oxygen in its higher range
# keeping smoothing for oxygen at 0.1

# final model for BW29 at CI has temperature, FSLE, salinity, and oxygen
# AIC: 226.8324, deviance explained: 49%
CI_final <- gam(BW29 ~ s(temperature_0,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) +
                s(salinity_0,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1), 
              weights=weights, family=binomial, data=KGI_binned)


# ------------------ Step 6: Visualize GAMs -------------------
# Function to create a cleaner visualization of a GAM model
visualizeGAM <- function(gam, predictors, siteName) {
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
    current_plot$label <- paste0("p-value = ", round(current_p_val,5))
    
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
    row <- 2
    col <- 1
  } else if(length(predictors) == 3 || length(predictors) == 4) {
    row <- 2
    col <- 2
  } else if(length(predictors) == 5 || length(predictors) == 6) {
    row <- 3
    col <- 2
  }
  
  # Aggregating all the plots into one figure
  final_plot <- wrap_plots(all_plots, nrow = row, ncol = col, guides = "collect") &
    plot_annotation(title = paste0(name(species), " at ", siteName,
                                   " (",deviance,"% Deviance Explained)"))
  
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
KGI_pred <- c('ice_conc','o2_0')
KGI_plots <- visualizeGAM(KGI_final, KGI_pred, 'KGI')

EI_pred <- c('mixed_layer')
EI_plots <- visualizeGAM(EI_final, EI_pred, 'EI')


CI_pred <- c('temperature_0','FSLE','salinity_0','o2_0')
CI_plots <- visualizeGAM(CI_final, CI_pred, 'CI')
