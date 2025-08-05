library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)

# ------------- Step 0: Choose Species ----------------
# Modeling Pm for all sites it is present (EI & CI), 40 km radius environmental data
species <- c('BW37') # options: BW29, BW37, Oo, Pm, Gm
# BW29 = Southern bottlenose whale, BW37 = Gray's and strap-toothed whales
# Oo = Killer whale, Pm = Sperm Whale, Gm = Long-finned pilot whale
# Note: not enough data to model BW58 at any site
sites <- c('EI','KGI','CI') # might not model CI, unclear if enough data

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
# 2.675251       2.753947       8.732690       3.633710       1.600542       4.390069      58.158197       1.909625 
# o2_0         chla_0 productivity_0 
# 40.316700      19.626474      20.767076 

# removing temperature
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
              "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0           o2_0 
# 2.682726       2.105412       6.631590       2.562893       1.521564       3.992775       1.565064      20.397754 
# chla_0 productivity_0 
# 14.416170      19.281382 

# removing oxygen concentration
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0", 'chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0         chla_0 
# 1.858249       1.934275       4.160423       2.313631       1.429246       3.732810       1.358885      13.517687 
# productivity_0 
# 8.850236 

# remove chlorophyll
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0", 'productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0 productivity_0 
# 1.869524       1.221080       3.015767       2.082012       1.352137       2.970050       1.370421       5.951237

# remove primary production
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0")
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE         SSH mixed_layer    ice_conc    ice_diff  salinity_0       EKE_0 
# 1.246420    1.217921    1.771057    1.847067    1.303808    1.382892    1.345480 

# final predictors for EI: FSLE, SSH, mixed layer depth, ice concentration, difference in ice, salinity, EKE

# KING GEORGE ISLAND
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0', 'julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
# logistic model did not converge, doing a regular linear model for VIF
KGI_vif <- glm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0          EKE_0 
# 5.204969      13.600947      12.158026      23.756336       2.412271      26.896812      42.386452       3.375282 
# o2_0         chla_0 productivity_0     julian_day 
# 14.477522      42.417088      47.234936       7.714122 

# removing primary production
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
              "temperature_0", "EKE_0",'o2_0','chla_0','julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# FSLE           SSH   mixed_layer      ice_conc      ice_diff    salinity_0 temperature_0         EKE_0          o2_0 
# 3.017585     13.597898     12.103759     18.533760      2.164271     22.225327     23.484862      2.834284      9.234252 
# chla_0    julian_day 
# 10.321175      6.188530 

# removing temperature
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
              "EKE_0",'o2_0','chla_0','julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# FSLE         SSH mixed_layer    ice_conc    ice_diff  salinity_0       EKE_0        o2_0      chla_0  julian_day 
# 2.116292   13.152098    8.740210   15.123865    1.869976   21.818307    2.587255    4.161131    9.392733    4.744134

# removing salinity
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',  
              "EKE_0",'o2_0','chla_0','julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# FSLE         SSH mixed_layer    ice_conc    ice_diff       EKE_0        o2_0      chla_0  julian_day 
# 2.115521   11.500794    5.241359   11.236468    1.865719    2.394771    3.933761    8.828312    3.981740 

# removing sea ice concentration
KGI_pred <- c("FSLE", "SSH", "mixed_layer", 'ice_diff',  
              "EKE_0",'o2_0','chla_0','julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)

# removing chlorophyll
KGI_pred <- c("FSLE", "SSH", "mixed_layer", 'ice_diff',  
              "EKE_0",'o2_0','julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# FSLE         SSH mixed_layer    ice_diff       EKE_0        o2_0  julian_day 
# 2.101461    3.325650    2.510340    1.844511    1.721748    1.988136    3.682233 

# final predictors for KGI: FSLE, SSH, mixed layer depth, difference in ice, EKE, oxygen, julian day

# CLARENCE ISLAND
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
# logistic model did not converge, doing a regular linear model for VIF
CI_vif <- glm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0          EKE_0 
# 2.978411       5.678823       4.314468       2.714885       1.562568       6.494601      15.370744       1.240458 
# o2_0         chla_0 productivity_0 
# 6.143649       6.089811      18.944454 

# remove primary production
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# FSLE           SSH   mixed_layer      ice_conc      ice_diff    salinity_0 temperature_0         EKE_0          o2_0 
# 2.966221      5.662829      4.092079      2.364949      1.562566      6.074576      4.594500      1.240456      3.834732 
# chla_0 
# 2.616437

# remove salinity
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',  
             "temperature_0", "EKE_0",'o2_0','chla_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# FSLE           SSH   mixed_layer      ice_conc      ice_diff temperature_0         EKE_0          o2_0        chla_0 
# 2.862442      2.699328      2.773518      2.305051      1.551299      4.269732      1.219936      3.789067      2.519791 

# predictors for CI: FSLE, SSH, mixed layer depth, ice concentration, difference in ice, temperature, EKE, oxygen, chlorophyll

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
# NOT weighing for BW37 at EI since roughly half of bins have presence
# List of predictors are stored in EI_pred include FSLE, SSH, mixed layer depth, ice concentration,
#   difference in ice concentration, salinity, eddy kinetic energy

# FSLE
EI_gam <- gam(BW37 ~ s(FSLE,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW37 ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2144     0.2942  -0.729    0.466
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(FSLE) 1.526  1.866  0.377   0.782
# 
# R-sq.(adj) =  -0.0119   Deviance explained = 1.58%
# UBRE = 0.46068  Scale est. = 1         n = 47

# sea surface height
EI_gam <- gam(BW37 ~ s(SSH,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW37 ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2141     0.2939  -0.728    0.466
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(SSH) 1.59  1.949  0.131   0.924
# 
# R-sq.(adj) =  -0.028   Deviance explained = 0.548%
# UBRE = 0.47762  Scale est. = 1         n = 47

# mixed layer depth
EI_gam <- gam(BW37 ~ s(mixed_layer,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW37 ~ s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2162     0.3024  -0.715    0.475
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(mixed_layer) 1.548  1.885   2.59   0.179
# 
# R-sq.(adj) =  0.0721   Deviance explained = 7.85%
# UBRE = 0.37543  Scale est. = 1         n = 47

# sea ice concentration
EI_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2182     0.2950   -0.74     0.46
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.341  1.594  0.526   0.767
# 
# R-sq.(adj) =  -0.018   Deviance explained = 0.997%
# UBRE = 0.46085  Scale est. = 1         n = 47

# sea ice difference
EI_gam <- gam(BW37 ~ s(ice_diff,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW37 ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2063     0.2983  -0.692    0.489
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_diff) 1.411  1.698  1.021   0.391
# 
# R-sq.(adj) =  0.0381   Deviance explained =  5.2%
# UBRE = 0.40599  Scale est. = 1         n = 47

# salinity
EI_gam <- gam(BW37 ~ s(salinity_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW37 ~ s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2289     0.3009  -0.761    0.447
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(salinity_0) 1.452  1.772  1.927   0.346
# 
# R-sq.(adj) =  0.00929   Deviance explained = 3.15%
# UBRE = 0.43597  Scale est. = 1         n = 47

# eddy kinetic energy
EI_gam <- gam(BW37 ~ s(EKE_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   BW37 ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.2137     0.2935  -0.728    0.467
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.584  1.921  0.084   0.961
# 
# R-sq.(adj) =  -0.0329   Deviance explained = 0.195%
# UBRE = 0.48225  Scale est. = 1         n = 47

# nothing is significant, look into this later

# -------------------- Step 5b: King George Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing at 0.1
# NOT weighing for BW37 at KGI since over half of bins have presence
# List of predictors are stored in EI_pred include FSLE, SSH, mixed layer depth, difference in ice,
#   eddy kinetic energy, oxygen, julian day

# FSLE
KGI_gam <- gam(BW37 ~ s(FSLE,k=4,sp=0.1), family=binomial, data=KGI_binned)
# Formula:
#   BW37 ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.7452     0.5437   1.371     0.17
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(FSLE) 1.247  1.449  3.771  0.0646 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.232   Deviance explained = 21.1%
# UBRE = 0.23099  Scale est. = 1         n = 24

# sea surface height
KGI_gam <- gam(BW37 ~ s(SSH,k=4,sp=0.1), family=binomial, data=KGI_binned)
# Formula:
#   BW37 ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.6164     0.6092   1.012    0.312
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(SSH) 1.18  1.337  6.719  0.0196 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.444   Deviance explained = 40.8%
# UBRE = -0.034892  Scale est. = 1         n = 24

# mixed layer depth
KGI_gam <- gam(BW37 ~ s(mixed_layer,k=4,sp=0.1), family=binomial, data=KGI_binned)
# Formula:
#   BW37 ~ s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.6181     0.5403   1.144    0.253
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(mixed_layer) 1.194  1.364  6.608  0.0238 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.329   Deviance explained = 30.1%
# UBRE = 0.10784  Scale est. = 1         n = 24

# difference in ice concentration
KGI_gam <- gam(BW37 ~ s(ice_diff,k=4,sp=0.1), family=binomial, data=KGI_binned)
# Formula:
#   BW37 ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.5181     0.4281    1.21    0.226
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_diff) 1.432  1.707  1.254   0.488
# 
# R-sq.(adj) =  0.0548   Deviance explained = 8.34%
# UBRE = 0.41541  Scale est. = 1         n = 24

# eddy kinetic energy
KGI_gam <- gam(BW37 ~ s(EKE_0,k=4,sp=0.1), family=binomial, data=KGI_binned)
# Formula:
#   BW37 ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.5345     0.4381    1.22    0.222
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.38  1.657   1.46   0.389
# 
# R-sq.(adj) =  0.00771   Deviance explained =  5.3%
# UBRE = 0.45129  Scale est. = 1         n = 24

# oxygen concentration
KGI_gam <- gam(BW37 ~ s(o2_0,k=4,sp=0.1), family=binomial, data=KGI_binned)
# Formula:
#   BW37 ~ s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.5680     0.4489   1.265    0.206
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(o2_0) 1.32  1.563  1.818   0.297
# 
# R-sq.(adj) =  0.024   Deviance explained = 6.56%
# UBRE = 0.42973  Scale est. = 1         n = 24

# julian day
KGI_gam <- gam(BW37 ~ s(julian_day,k=4,sp=0.1), family=binomial, data=KGI_binned)
# Formula:
#   BW37 ~ s(julian_day, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.5771     0.4508    1.28      0.2
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(julian_day) 1.388  1.666  3.035   0.289
# 
# R-sq.(adj) =  0.161   Deviance explained = 15.7%
# UBRE = 0.31444  Scale est. = 1         n = 24

# significant variables: SSH, mixed layer depth

# combining both into one model
KGI_gam <- gam(BW37 ~ s(SSH,k=4,sp=0.1) + s(mixed_layer,k=4,sp=0.1),
               family=binomial, data=KGI_binned)
# AIC: 20.35083
# summary:
# Formula:
#   BW37 ~ s(SSH, k = 4, sp = 0.1) + s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.8666     0.7369   1.176     0.24
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(SSH)         1.184  1.340  5.378  0.0487 *
#   s(mixed_layer) 1.141  1.267  3.715  0.1075  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.686   Deviance explained = 56.9%
# UBRE = -0.15205  Scale est. = 1         n = 24
# ................................................
# even though mixed layer depth currently not significant, it has minimal concurvity with SSH
#   it also contributes to lower AIC and higher deviance explained
# going to see if changing smoothing parameter makes it significant

# changing smoothing parameter for mixed layer depth
KGI_gam <- gam(BW37 ~ s(SSH,k=4,sp=0.1) + s(mixed_layer,k=4,sp=1),
               family=binomial, data=KGI_binned)
# AIC: 20.46859
# summary:
# Formula:
#   BW37 ~ s(SSH, k = 4, sp = 0.1) + s(mixed_layer, k = 4, sp = 1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.9113     0.7358   1.239    0.216
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(SSH)         1.185  1.342  5.263  0.0512 .
# s(mixed_layer) 1.016  1.032  3.327  0.0745 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.682   Deviance explained = 55.7%
# UBRE = -0.14714  Scale est. = 1         n = 24
# ..........................................
# lower smoothing
KGI_gam <- gam(BW37 ~ s(SSH,k=4,sp=0.1) + s(mixed_layer,k=4,sp=0.01),
               family=binomial, data=KGI_binned)
# Formula:
#   BW37 ~ s(SSH, k = 4, sp = 0.1) + s(mixed_layer, k = 4, sp = 0.01)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.6976     0.7560   0.923    0.356
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(SSH)         1.176  1.326  5.829  0.0392 *
#   s(mixed_layer) 1.645  2.036  3.289  0.2048  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =    0.7   Deviance explained = 61.2%
# UBRE = -0.16778  Scale est. = 1         n = 24
#...............................................
# higher smoothing for both variables
KGI_gam <- gam(BW37 ~ s(SSH,k=4,sp=1) + s(mixed_layer,k=4,sp=1),
               family=binomial, data=KGI_binned)
# Formula:
#   BW37 ~ s(SSH, k = 4, sp = 1) + s(mixed_layer, k = 4, sp = 1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.9448     0.7302   1.294    0.196
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(SSH)         1.023  1.045  4.949  0.0297 *
#   s(mixed_layer) 1.016  1.032  3.231  0.0787 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.681   Deviance explained = 55.1%
# UBRE = -0.15273  Scale est. = 1         n = 24
# ..............................................
# dropping mixed layer depth because still not significant

# model has just SSH, changing smoothing & basis functions
KGI_gam <- gam(BW37 ~ s(SSH,k=4,sp=1),
               family=binomial, data=KGI_binned) # AIC: 22.8688
KGI_gam <- gam(BW37 ~ s(SSH,k=4,sp=0.01),
               family=binomial, data=KGI_binned) # AIC: 24.37222
# changing smoothing to 1, even though AIC increases the partial effect plot becomes more interpretable
KGI_gam <- gam(BW37 ~ s(SSH,k=8,sp=1),
               family=binomial, data=KGI_binned) # AIC: 22.97687
# negligile impact on confidence interval, AIC, so keeping knots at 4

# final model for KGI (AIC = 22.8688, deviance explained 40.7%)
KGI_final <- gam(BW37 ~ s(SSH,k=4,sp=1),
               family=binomial, data=KGI_binned)

# -------------------- Step 5c: Clarence Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing at 0.1
# might not end up modeling CI because less presence (74 bins without presence, only 7 with)
# List of predictors are stored in EI_pred include FSLE, SSH, mixed layer depth, ice concenctraion, difference in ice,
#   temperature, eddy kinetic energy, oxygen, chlorophyll

# setting weights for CI, weighing presence at 10 to make ratio of 1s to 0s roughly 1:1
CI_binned$weights <- ifelse(CI_binned$BW37 == 1, 10, 1)

# FSLE
CI_gam <- gam(BW37 ~ s(FSLE,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW37 ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.4350     0.2093  -2.078   0.0377 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE) 1.888  2.263   20.9 6.15e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.127   Deviance explained = 14.2%
# UBRE = 1.1837  Scale est. = 1         n = 81

# sea surface height
CI_gam <- gam(BW37 ~ s(SSH,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW37 ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.5001     0.2220  -2.253   0.0243 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(SSH) 1.567   1.91  15.75 0.000223 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.159   Deviance explained = 16.7%
# UBRE = 1.1152  Scale est. = 1         n = 81

# mixed layer depth
CI_gam <- gam(BW37 ~ s(mixed_layer,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW37 ~ s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.4003     0.2164   -1.85   0.0643 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(mixed_layer) 1.484   1.79  14.83 0.00135 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0844   Deviance explained = 10.7%
# UBRE =   1.26  Scale est. = 1         n = 81

# difference in ice
CI_gam <- gam(BW37 ~ s(ice_diff,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW37 ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1329     0.1725   -0.77    0.441
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_diff) 1.422   1.72  3.617   0.113
# 
# R-sq.(adj) =  0.0407   Deviance explained = 4.72%
# UBRE = 1.4066  Scale est. = 1         n = 81

# eddy kinetic energy
CI_gam <- gam(BW37 ~ s(EKE_0,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW37 ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.07282    0.16787  -0.434    0.664
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.909  2.328  1.197   0.651
# 
# R-sq.(adj) =  -0.0119   Deviance explained = 1.03%
# UBRE = 1.5095  Scale est. = 1         n = 81

# oxygen concentration
CI_gam <- gam(BW37 ~ s(o2_0,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW37 ~ s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.4040     0.1974  -2.047   0.0407 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(o2_0) 1.965  2.337  22.62 2.11e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.16   Deviance explained = 16.7%
# UBRE = 1.1259  Scale est. = 1         n = 81

# temperature
CI_gam <- gam(BW37 ~ s(temperature_0,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW37 ~ s(temperature_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.6340     0.4096  -3.989 6.64e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(temperature_0) 1.728  2.074  39.69  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.486   Deviance explained = 47.2%
# UBRE = 0.36873  Scale est. = 1         n = 81

# sea ice concentration
CI_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -136.8      109.1  -1.254     0.21
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc)   1      1  1.589   0.207
# 
# R-sq.(adj) =  0.322   Deviance explained = 36.5%
# UBRE = 0.61241  Scale est. = 1         n = 81

# chlorophyll
CI_gam <- gam(BW37 ~ s(chla_0,k=4,sp=0.1), weights=weights, family=binomial, data=CI_binned)
# Formula:
#   BW37 ~ s(chla_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.8798     0.2614  -3.366 0.000762 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(chla_0) 1.767  2.154  36.74  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.345   Deviance explained = 32.9%
# UBRE = 0.72042  Scale est. = 1         n = 81

# significant predictors for CI: FSLE, SSH, mixed layer depth, oxygen concentration, temperature, chlorophyll
# adding these variables one at a time to develop final model

# FSLE and SSH
CI_gam <- gam(BW37 ~ s(FSLE,k=4,sp=0.1) + s(SSH,k=4,sp=0.1),
              weights=weights, family=binomial, data=CI_binned)
# AIC: 152.8525
# Formula:
# BW37 ~ s(FSLE, k = 4, sp = 0.1) + s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.7697     0.2472  -3.113  0.00185 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(FSLE) 1.869  2.226  15.66 0.00183 **
#   s(SSH)  1.512  1.837  17.09 0.00384 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.268   Deviance explained = 27.8%
# UBRE = 0.88707  Scale est. = 1         n = 81
# ..................................................
# concurvity and plots not problematic

# adding mixed layer depth
CI_gam <- gam(BW37 ~ s(FSLE,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) + s(mixed_layer,k=4,sp=0.1),
              weights=weights, family=binomial, data=CI_binned)
# AIC: 149.8564
# Formula:
#   BW37 ~ s(FSLE, k = 4, sp = 0.1) + s(SSH, k = 4, sp = 0.1) + s(mixed_layer, 
#                                                                 k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.8664     0.2591  -3.344 0.000825 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(FSLE)        1.857  2.209 14.460 0.00315 **
#   s(SSH)         1.525  1.855 12.367 0.00797 **
#   s(mixed_layer) 1.480  1.785  5.071 0.14657   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.289   Deviance explained = 30.8%
# UBRE = 0.85008  Scale est. = 1         n = 81
# .................................................
# plots & concurvities look alright
# removing mixed layer because not significant

# removing mixed layer, adding in oxygen concentration
CI_gam <- gam(BW37 ~ s(FSLE,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1),
              weights=weights, family=binomial, data=CI_binned)
# AIC: 142.4519
# summary:
# Formula:
#   BW37 ~ s(FSLE, k = 4, sp = 0.1) + s(SSH, k = 4, sp = 0.1) + s(o2_0, 
#                                                                 k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.0139     0.2884  -3.515  0.00044 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(FSLE) 1.781  2.113  9.619 0.01306 * 
#   s(SSH)  1.472  1.775 14.404 0.00995 **
#   s(o2_0) 1.759  2.096  7.655 0.01518 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.337   Deviance explained = 34.6%
# UBRE = 0.75867  Scale est. = 1         n = 81
# ................................................
# plots and concurvities look reasonable

# adding in temperature
CI_gam <- gam(BW37 ~ s(FSLE,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) + s(o2_0,k=4,sp=0.1) +
                s(temperature_0,k=4,sp=0.1),
              weights=weights, family=binomial, data=CI_binned)
# # AIC: 110.2529
# Formula:
#   BW37 ~ s(FSLE, k = 4, sp = 0.1) + s(SSH, k = 4, sp = 0.1) + s(o2_0, 
#                                                                 k = 4, sp = 0.1) + s(temperature_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7518     0.4366  -4.013    6e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE)          1.599  1.885  0.080    0.948    
# s(SSH)           1.390  1.648  4.270    0.159    
# s(o2_0)          1.396  1.639  0.133    0.870    
# s(temperature_0) 1.619  1.886 18.910 5.11e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.52   Deviance explained = 51.8%
# UBRE = 0.36115  Scale est. = 1         n = 81
# ................................................
# temperature has near problematic concurvity with oxygen
# will try removing oxygen first

# removing oxygen
CI_gam <- gam(BW37 ~ s(FSLE,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) +
                s(temperature_0,k=4,sp=0.1),
              weights=weights, family=binomial, data=CI_binned)
# AIC: 109.326
# Formula:
#   BW37 ~ s(FSLE, k = 4, sp = 0.1) + s(SSH, k = 4, sp = 0.1) + s(temperature_0, 
#                                                                 k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7679     0.4252  -4.158 3.21e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(FSLE)          1.661  1.974  0.327   0.846    
# s(SSH)           1.402  1.674  5.142   0.137    
# s(temperature_0) 1.649  1.944 23.038 6.6e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.517   Deviance explained = 50.9%
# UBRE = 0.3497  Scale est. = 1         n = 81
# ....................................................
# keeping oxygen out of model because of lower AIC without it
# FSLE has least significant effect, removing it next

# removing FSLE
CI_gam <- gam(BW37 ~ s(SSH,k=4,sp=0.1) +
                s(temperature_0,k=4,sp=0.1),
              weights=weights, family=binomial, data=CI_binned)
# AIC: 106.3866
# Formula:
#   BW37 ~ s(SSH, k = 4, sp = 0.1) + s(temperature_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7373     0.4131  -4.206  2.6e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(SSH)           1.418  1.701  5.262   0.132    
# s(temperature_0) 1.703  2.034 34.212  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.528   Deviance explained = 50.8%
# UBRE = 0.31342  Scale est. = 1         n = 81
# ............................................
# model got better without keeping FSLE
# if remove SSH, AIC will go up (110.8671 based on model of just temp at line 789)
# will keep SSH and temperature for now, see if adding another variable changes

# adding chlorophyll
CI_gam <- gam(BW37 ~ s(SSH,k=4,sp=0.1) + s(chla_0,k=4,sp=0.1) + 
                s(temperature_0,k=4,sp=0.1),
              weights=weights, family=binomial, data=CI_binned)
# AIC: 106.5035
# Formula:
#   BW37 ~ s(SSH, k = 4, sp = 0.1) + s(chla_0, k = 4, sp = 0.1) + 
#   s(temperature_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.8303     0.4845  -3.777 0.000159 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(SSH)           1.391  1.662  3.599 0.211041    
# s(chla_0)        1.359  1.608  1.021 0.345442    
# s(temperature_0) 1.612  1.897 16.050 0.000188 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.531   Deviance explained =   52%
# UBRE = 0.31486  Scale est. = 1         n = 81
# ..............................................
# concurvity is not problematic
# will drop SSH first to see if chlorophyll becomes significant

# dropping SSH
CI_gam <- gam(BW37 ~ s(chla_0,k=4,sp=0.1) + 
                s(temperature_0,k=4,sp=0.1),
              weights=weights, family=binomial, data=CI_binned)
# AIC: 108.7752
# Formula:
#   BW37 ~ s(chla_0, k = 4, sp = 0.1) + s(temperature_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.6738     0.4298  -3.894 9.86e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(chla_0)        1.387  1.659  1.639     0.23    
# s(temperature_0) 1.643  1.947 22.478 8.74e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.505   Deviance explained = 49.5%
# UBRE = 0.3429  Scale est. = 1         n = 81
# ................................................
# chlorophyll still not significant, dropping it from model

# dropping chlorophyll, only keeping temperature
# refer to line 789 for output
# now, experimenting with different knots/smoothing
CI_gam <- gam(BW37 ~ s(temperature_0,k=4,sp=0.01),
              weights=weights, family=binomial, data=CI_binned) # AIC: 100.9191
CI_gam <- gam(BW37 ~ s(temperature_0,k=4,sp=1),
              weights=weights, family=binomial, data=CI_binned) # AIC: 140.6829
# changin sp to 1 since trend is stronger and AIC is much lower
# checking knots
CI_gam <- gam(BW37 ~ s(temperature_0,k=8,sp=1),
              weights=weights, family=binomial, data=CI_binned) # AIC: 123.3664
# AIC increased, model is not as good

# final model for CI (AIC = 100.9191, 52.6% explained)
CI_final <- gam(BW37 ~ s(temperature_0,k=4,sp=0.01),
              weights=weights, family=binomial, data=CI_binned)

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
    plot_annotation(title = paste0("Long-Finned Pilot Whale at ", siteName,
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
KGI_pred <- c('SSH')
KGI_plots <- visualizeGAM(KGI_final, KGI_pred, 'KGI')

CI_pred <- c('temperature_0')
CI_plots <- visualizeGAM(CI_final, CI_pred, 'CI')
