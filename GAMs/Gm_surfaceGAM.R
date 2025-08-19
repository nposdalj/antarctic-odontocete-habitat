library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)
# MODELS THAT HAVE BEEN IMPROVED: CI
# MODELS IN INITIAL PASSTHROUGH: EI, KGI

# ------------- Step 0: Choose Species ----------------
# Modeling Gm for all sites, 40 km radius environmental data
species <- c('Gm') # options: BW29, BW37, Oo, Pm, Gm
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
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime,fsle_orient))
} else if(species =='BW37') {
  depths <- c(0, 67, 920) 
  sp_specific <- allData %>% subset(select=-c(BW29,BW58,Oo,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_67,salinity_67,EKE_67,
                    temperature_920,salinity_920,EKE_920, chla_0,o2_0,productivity_0,chla_67,
                    o2_67,productivity_67, chla_920,o2_920,productivity_920,
                    ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime,fsle_orient))
} else if(species =='Oo') {
  depths <- c(0, 11, 455) 
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_11,salinity_11,EKE_11,
                    temperature_455,salinity_455,EKE_455, chla_0,o2_0,productivity_0,chla_11,
                    o2_11,productivity_11, chla_455,o2_455,productivity_455,
                    ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime,fsle_orient))
} else if(species =='Pm') {
  depths <- c(0, 375, 1665)
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Oo,Gm)) %>%
    subset(select=c(date,julian_day,Site,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_375,salinity_375,EKE_375,
                    temperature_1665,salinity_1665,EKE_1665, chla_0,o2_0,productivity_0,chla_375,
                    o2_375,productivity_375, chla_1665,o2_1665,productivity_1665,
                    ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime,fsle_orient))
} else if(species =='Gm') {
  depths <- c(0, 16, 635) 
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Oo,Pm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_16,salinity_16,EKE_16,
                    temperature_635,salinity_635,EKE_635, chla_0,o2_0,productivity_0,chla_16,
                    o2_16,productivity_16, chla_635,o2_635,productivity_635,
                    ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                    chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime,fsle_orient))
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
                         salinity_0_sd = mean_col('salinity_sd_0'), fsle_orient = mean_col('fsle_orient'))
  
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
EI_acf <- acfVal('EI')
EI_binned <- binByACF('EI',EI_acf)
KGI_acf <- acfVal('KGI')
KGI_binned <- binByACF('KGI',KGI_acf)
CI_acf <- acfVal('CI')
CI_binned <- binByACF('CI',CI_acf)

# ACF binned (2 days) has 154 0s, 27 1s (14% 1s)
KGI_3day <- binByACF('KGI',3) # 101 0s, 22 1s (18% 1s)
KGI_4day <- binByACF('KGI',4) # 74 0s, 19 1s (20% 1s)
KGI_5day <- binByACF('KGI',5) # 60 0s, 17 1s (22% 1s)
KGI_6day <- binByACF('KGI',6) # 50 0s, 15 1s (23% 1s)
KGI_week <- binByACF('KGI',7) # 44 0s, 14 1s (24% 1s)


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
# Only including julian day for KGI, since it has almost 1 year of data

# ELEPHANT ISLAND
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
                     "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
# EI vif model is not converging when run as a logistic regression, so running as a linear regression instead
EI_vif <- glm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0  temperature_0          EKE_0 
# 3.917328       2.796444       7.806685       5.448504       3.456836       7.245168      65.689317       2.110835 
# o2_0         chla_0 productivity_0 
# 45.948050      23.479121      32.495525 

# dropping surface temperature
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
                     "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0          EKE_0           o2_0 
# 3.635886       2.146968       7.606461       3.349641       3.140330       6.595141       1.882706      26.562145 
# chla_0 productivity_0 
# 20.869381      32.205804 

# dropping primary production
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
             "EKE_0",'chla_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# FSLE         SSH mixed_layer    ice_conc fsle_orient  salinity_0       EKE_0      chla_0 
# 1.702349    1.803669    4.519144    2.940721    2.375007    4.464684    1.501540   11.582560 

# dropping chlorophyll
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
             "EKE_0")
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# FSLE         SSH mixed_layer    ice_conc fsle_orient  salinity_0       EKE_0 
# 1.381906    1.296292    2.166030    2.286215    2.315872    1.817996    1.450946 

# Final predictors for Elephant Island: FSLE, SSH, mixed layer depth, sea ice concentraiton,
#   fsle orientation, salinity, EKE

# KING GEORGE ISLAND
# adding julian day into initial formula
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
                     "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0', 'julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), family = binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0  temperature_0          EKE_0 
# 1.406687       2.433127       7.188120       2.958964       1.598854       4.607197      40.902517       1.140619 
# o2_0         chla_0 productivity_0     julian_day 
# 14.993033       9.906773      19.276807       3.296272 

# removing temp_0
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
                    "EKE_0",'o2_0','chla_0','productivity_0', 'julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), family = binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0
# 1.454792       2.028160       4.876820       2.908102       1.357751       4.092606       1.072813
# o2_0         chla_0 productivity_0     julian_day
# 2.755725       7.688033       6.222557       3.105939

# removing chlorophyll
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
              "EKE_0",'o2_0','productivity_0', 'julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), family = binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0          EKE_0           o2_0 
# 1.336210       2.189414       4.166698       2.339382       1.278746       4.122134       1.094359       2.521988 
# productivity_0     julian_day 
# 2.461906       2.533958 

# Final predictors for KGI: FSLE, SSH, mixed layer depth, sea ice concentration, fsle orientation
#  salinity, EKE, oxygen, net primary production, julian day

# CLARENCE ISLAND
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
                     "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0','EKE_0_mad')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), family = binomial, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0  temperature_0 
# 3.519468       6.640285       6.328253       3.032206       3.026913      12.325767      37.402842 
# EKE_0           o2_0         chla_0 productivity_0      EKE_0_mad 
# 1.511918       6.699801      13.700005      20.432526       5.350101 

# removing SST
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
                     "EKE_0",'o2_0','chla_0','productivity_0','EKE_0_mad')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), family = binomial, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0          EKE_0 
# 3.125057       6.710290       5.933576       2.805004       2.680816      10.635756       1.509862 
# o2_0         chla_0 productivity_0      EKE_0_mad 
# 3.282087      13.285292       7.511489       4.783141 

# removing chlorophyll
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
             "EKE_0",'o2_0','productivity_0','EKE_0_mad')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), family = binomial, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0          EKE_0 
# 2.565896       4.642922       4.105160       2.062533       1.964091       7.405345       1.406003 
# o2_0 productivity_0      EKE_0_mad 
# 2.571408       3.053382       3.212377

# removing salinity
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient',  
             "EKE_0",'o2_0','productivity_0','EKE_0_mad')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), family = binomial, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient          EKE_0           o2_0 
# 2.255048       2.373907       2.645606       1.943394       1.978269       1.378912       2.405171 
# productivity_0      EKE_0_mad 
# 2.936078       2.798508 

# Final Clarence Island predictors: FSLE, SSH, mixed layer thickness, sea ice concentration,
#   FSLE orientation, EKE, oxygen concentration, net primary production

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
# starting by building GAMs one predictor at a time to find significant variables
# Similar process: weighing to reach roughly 1:1 1s to 0s and setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in EI_pred include FSLE, SSH, mixed layer, ice concentration, 
#   ice concentration difference, EKE, salinity

# not adding weights for EI because ratio of 0s to 1s is already close to 1:1 (13:17)
# starting with FSLE magnitude and orientation
EI_gam <- gam(Gm ~ s(FSLE,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Gm ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.2756     0.3732   0.738     0.46
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(FSLE) 1.43  1.731  0.609   0.569
# 
# R-sq.(adj) =  0.0247   Deviance explained = 5.24%
# UBRE = 0.45869  Scale est. = 1         n = 30

# sea surface height
EI_gam <- gam(Gm ~ s(SSH,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Gm ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.2696     0.3695    0.73    0.466
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(SSH) 1.488  1.809  0.129   0.916
# 
# R-sq.(adj) =  -0.0357   Deviance explained = 1.27%
# UBRE = 0.51703  Scale est. = 1         n = 30
# ................................................

# mixed layer depth
EI_gam <- gam(Gm ~ s(mixed_layer,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Gm ~ s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.6672     0.5397   1.236    0.216
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(mixed_layer) 1.044  1.087  5.054  0.0236 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.301   Deviance explained = 29.7%
# UBRE = 0.097976  Scale est. = 1         n = 30

# ice concentration
EI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.4991     0.4903   1.018    0.309
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.167  1.311  2.644   0.114
# 
# R-sq.(adj) =  0.144   Deviance explained = 14.8%
# UBRE = 0.31068  Scale est. = 1         n = 30

# salinity (LOG TRANSFORM)
EI_gam <- gam(Gm ~ s(salinity_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary
# Formula:
#   Gm ~ s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.2560     0.4018   0.637    0.524
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(salinity_0) 1.249  1.457  3.624   0.113
# 
# R-sq.(adj) =  0.0977   Deviance explained = 11.4%
# UBRE = 0.36219  Scale est. = 1         n = 30

# eddy kinetic energy
EI_gam <- gam(Gm ~ s(EKE_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Gm ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.2722     0.3871   0.703    0.482
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.415  1.712  2.722   0.229
# 
# R-sq.(adj) =  0.112   Deviance explained = 11.2%
# UBRE = 0.37621  Scale est. = 1         n = 30

# The only significant predictor is mixed layer depth
# playing with different smoothing parameters, basis functions to find best model
EI_gam <- gam(Gm ~ s(mixed_layer,k=4,sp=0.1), family=binomial, data=EI_binned)
# AIC: 32.93928
# gam.check does not indicate that k should be increased, trying it anyways:
EI_gam <- gam(Gm ~ s(mixed_layer,k=8,sp=0.1), family=binomial, data=EI_binned)
# AIC: 33.1593

# AIC went up, reverting to k = 4
# trying sp = 0.01
EI_gam <- gam(Gm ~ s(mixed_layer,k=8,sp=0.01), family=binomial, data=EI_binned)
# AIC: 34.43422
# AIC went up, trying sp = 1
EI_gam <- gam(Gm ~ s(mixed_layer,k=8,sp=1), family=binomial, data=EI_binned)
# AIC: 32.91472
# negligible reduction in AIC, sticking with sp = 0.1

# final model: mixed layer thickness as only significant predictor
# Could be a result of less data (~30 data points vs 150+ for CI and KGI)
#   could also have to do with weighing - EI did not need to have detections weighed at all to reach 1:1 ratio
EI_final <- gam(Gm ~ s(mixed_layer,k=4,sp=0.1), family=binomial, data=EI_binned)

# -------------------- Step 5b: King George Island GAM -------------------------
# initial gam based on list of KGI predictors
KGI_gam <- gam(Gm ~ s(FSLE,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) + s(ice_conc,k=4) + 
                 s(ice_diff,k=4) + s(salinity_0,k=4) +
                 s(EKE_0,k=4) + s(o2_0,k=4) + s(productivity_0,k=4) + s(julian_day,k=4), 
               data = KGI_binned, family=binomial,method='REML')
# summary()
# 
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(FSLE, k = 4) + s(SSH, k = 4) + s(mixed_layer, k = 4) + 
#   s(ice_conc, k = 4) + s(ice_diff, k = 4) + s(salinity_0, k = 4) + 
#   s(EKE_0, k = 4) + s(o2_0, k = 4) + s(productivity_0, k = 4) + 
#   s(julian_day, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -27.92      22.88   -1.22    0.222
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(FSLE)           2.114  2.491  1.003 0.55932   
# s(SSH)            1.000  1.000  1.182 0.27698   
# s(mixed_layer)    1.000  1.000  0.999 0.31765   
# s(ice_conc)       1.854  1.979  1.201 0.54380   
# s(ice_diff)       1.000  1.000  1.967 0.16078   
# s(salinity_0)     1.000  1.000  0.004 0.95116   
# s(EKE_0)          1.620  1.904  0.888 0.60863   
# s(o2_0)           1.000  1.000  3.732 0.05338 . 
# s(productivity_0) 1.000  1.000  8.726 0.00314 **
#   s(julian_day)     1.000  1.000  0.821 0.36501   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.446   Deviance explained = 46.8%
# -REML = 40.187  Scale est. = 1         n = 181
# ....................................................................
# concurvity()
# para   s(FSLE)    s(SSH) s(mixed_layer) s(ice_conc) s(ice_diff) s(salinity_0)  s(EKE_0)
# worst    2.335937e-22 0.5864634 0.8963933      0.8580711   0.9781698   0.6751140     0.9866093 0.7140100
# observed 2.335937e-22 0.3744678 0.8620985      0.8393429   0.9401030   0.2933448     0.9572558 0.4814734
# estimate 2.335937e-22 0.5344532 0.8692273      0.8262123   0.9578391   0.3229273     0.9541493 0.4209709
# s(o2_0) s(productivity_0) s(julian_day)
# worst    0.9716414         0.9703266     0.9920461
# observed 0.9527751         0.9398072     0.8546254
# estimate 0.9541944         0.9029360     0.8914166
# .....................................................................
# gam.check()
# 
# Method: REML   Optimizer: outer newton
# full convergence after 19 iterations.
# Gradient range [-1.697423e-05,6.362028e-06]
# (score 40.1868 & scale 1).
# Hessian positive definite, eigenvalue range [1.443786e-06,0.3172399].
# Model rank =  31 / 31 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                     k'  edf k-index p-value    
# s(FSLE)           3.00 2.11    0.96    0.34    
# s(SSH)            3.00 1.00    0.92    0.17    
# s(mixed_layer)    3.00 1.00    0.97    0.33    
# s(ice_conc)       3.00 1.85    0.84    0.02 *  
#   s(ice_diff)       3.00 1.00    0.79  <2e-16 ***
#   s(salinity_0)     3.00 1.00    0.96    0.34    
# s(EKE_0)          3.00 1.62    1.05    0.78    
# s(o2_0)           3.00 1.00    0.83    0.01 ** 
#   s(productivity_0) 3.00 1.00    0.91    0.17    
# s(julian_day)     3.00 1.00    0.80  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Building GAM one predictor at a time because multivariate GAM above had difficult to interpret results
KGI_gam <- gam(Gm ~ s(FSLE,k=4), family = binomial, method = "REML", data = KGI_binned)
# plot does not show significance (horizontal line test), but does show shape
# summary(KGI_gam)
# 
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(FSLE, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7618     0.2116  -8.326   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(FSLE) 1.945  2.337  0.907   0.573
# 
# R-sq.(adj) =  0.0028   Deviance explained = 1.82%
# -REML = 77.199  Scale est. = 1         n = 181
# gam.check()
# Method: REML   Optimizer: outer newton
# full convergence after 3 iterations.
# Gradient range [-3.110683e-06,-3.110683e-06]
# (score 77.1986 & scale 1).
# Hessian positive definite, eigenvalue range [0.1435551,0.1435551].
# Model rank =  4 / 4 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#           k'  edf k-index p-value
# s(FSLE) 3.00 1.94    0.94    0.3

KGI_gam <- gam(Gm ~ s(SSH,k=4), family = binomial, method = "REML", data = KGI_binned)
# plot shows significant relationship (higher SSH, higher probability) 
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(SSH, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.5621     0.4052  -6.323 2.57e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(SSH)   1      1  15.61 7.77e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.14   Deviance explained = 19.3%
# -REML = 61.996  Scale est. = 1         n = 181
# ..................................................
# AIC(KGI_gam) = 127.1259
# Trying linear fit because edf is 1
KGI_gam <- gam(Gm ~ SSH, family = binomial, method = 'REML', data=KGI_binned)
# AIC() = 127.1255
# Very marginal difference, but going with linear fit for SSH

KGI_gam <- gam(Gm ~ s(mixed_layer,k=4), family = binomial, method = "REML", data = KGI_binned)
# large error bar at very end, but reasonable errors for mixed layers between 15 to 80 m
# Doesn't pass horizontal line test
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.9392     0.2633  -7.365 1.78e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(mixed_layer) 2.483    2.8   4.59   0.208
# 
# R-sq.(adj) =  0.0236   Deviance explained = 5.14%
# -REML = 75.615  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(ice_conc,k=4), family = binomial, method = "REML", data = KGI_binned)
# not great: full-range error bars for concentrations of ~0.25 and above (trendline makes sense however)
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -23.33      19.35  -1.205    0.228
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.853  1.981  4.139    0.15
# 
# R-sq.(adj) =  0.125   Deviance explained =   19%
# -REML = 62.293  Scale est. = 1         n = 181


KGI_gam <- gam(Gm ~ s(ice_diff,k=4), family = binomial, method = "REML", data = KGI_binned)
# error bars big towards beginning, but trend looks alright
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(ice_diff, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7902     0.2197  -8.149 3.68e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_diff) 1.864   2.22  2.152   0.417
# 
# R-sq.(adj) =  0.00503   Deviance explained = 2.11%
# -REML =  76.63  Scale est. = 1         n = 181
# gam.check()
# Method: REML   Optimizer: outer newton
# full convergence after 4 iterations.
# Gradient range [-1.467041e-05,-1.467041e-05]
# (score 76.6299 & scale 1).
# Hessian positive definite, eigenvalue range [0.2613356,0.2613356].
# Model rank =  4 / 4 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#               k'  edf k-index p-value    
# s(ice_diff) 3.00 1.86    0.68  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


KGI_gam <- gam(Gm ~ s(salinity_0,k=4), family = binomial, method = "REML", data = KGI_binned)
# borderline fails horizontal line test, but clear trend in salinity
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(salinity_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.4923     0.4385  -5.684 1.31e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(salinity_0) 2.162  2.535  11.02 0.00862 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0965   Deviance explained = 16.9%
# -REML = 65.595  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(EKE_0,k=4), family = binomial, method = "REML", data = KGI_binned)
# minor trend, good amount of variability
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(EKE_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7431     0.2089  -8.346   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.32  1.563  0.213   0.886
# 
# R-sq.(adj) =  -0.00509   Deviance explained = 0.314%
# -REML =  77.49  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(o2_0,k=4), family = binomial, method = "REML", data = KGI_binned)
# reasonable trend in plot, reasonable error bars
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.1809     0.3051  -7.149 8.75e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(o2_0) 2.499  2.822  12.91 0.00462 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0955   Deviance explained = 13.4%
# -REML = 69.103  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(productivity_0,k=4), family = binomial, method = "REML", data = KGI_binned)
# Clear trend in presence based off plot
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(productivity_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.0638     0.2869  -7.194 6.27e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(productivity_0)   1      1  9.116 0.00253 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0647   Deviance explained = 9.12%
# -REML = 70.042  Scale est. = 1         n = 181
# .....................................................
# AIC() = 142.5871
# Trying linear fit because edf is 1
KGI_gam <- gam(Gm ~ productivity_0, family = binomial, method = 'REML', data = KGI_binned)
# AIC() = 142,5868
# Marginal difference, but going with linear fit for NPP

KGI_gam <- gam(Gm ~ s(julian_day,k=4), family = binomial, method = "REML", data = KGI_binned)
# Clear seasonal trend, showing more spring/early summer presence
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(julian_day, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.5203     0.4072  -6.189 6.04e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(julian_day) 2.748  2.954  21.52 8.86e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.236   Deviance explained = 24.9%
# -REML = 61.037  Scale est. = 1         n = 181

# Variables that were significant (p < 0.05): SSH (linear), salinity, oxygen, net primary production (linear), julian day
# Variables that were insignifacnt: FSLE, mixed layer depth, EKE, ice concentration, ice concentration difference
# Still including ice concentration (p=0.15) and/or ice concentration difference (p=0.417) because ecologically important

# Next step: one by one, add significant variables to model
# Starting with linear terms
KGI_gam <- gam(Gm ~ SSH + productivity_0, family = binomial, method = 'REML', data = KGI_binned)
# AIC: 113.4565
# summary:
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ SSH + productivity_0
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    63.34362   15.71880   4.030 5.58e-05 ***
#   SSH            39.41786    9.69153   4.067 4.76e-05 ***
#   productivity_0 -0.14846    0.04641  -3.199  0.00138 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# R-sq.(adj) =  0.302   Deviance explained = 29.5%
# -REML = 53.206  Scale est. = 1         n = 181

# Adding sea ice concentration, building rest of model around it
KGI_gam <- gam(Gm ~ SSH + productivity_0 + s(ice_conc,k=4), family = binomial, method = 'REML', data = KGI_binned)
# AIC: 101.0181
# summary:
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ SSH + productivity_0 + s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -1.99801   30.53945  -0.065 0.947837    
# SSH            13.43412   13.05840   1.029 0.303587    
# productivity_0 -0.17350    0.04998  -3.472 0.000517 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.825  1.972   1.17   0.553
# 
# R-sq.(adj) =  0.381   Deviance explained = 40.3%
# -REML = 44.221  Scale est. = 1         n = 181
# ............................................................
# Adding weights to make the ratio of 1s to 0s near 1:1
KGI_binned$weights <- ifelse(KGI_binned$Gm == 1, 5, 1)
KGI_gam <- gam(Gm ~ SSH + productivity_0 + s(ice_conc,k=4), family = binomial, method = 'REML', 
               weights=weights,data = KGI_binned)
# AIC: 227.155 (jump is okay because weights changed)
# summary
# Formula:
#   Gm ~ SSH + productivity_0 + s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -12.84468   22.15035  -0.580    0.562    
# SSH              7.47595    7.70685   0.970    0.332    
# productivity_0  -0.15815    0.02648  -5.971 2.35e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.895  1.993  2.148   0.326
# 
# R-sq.(adj) =  0.463   Deviance explained = 45.5%
# -REML = 109.99  Scale est. = 1         n = 181

# Dropping SSH because it is no longer significant
KGI_gam <- gam(Gm ~ productivity_0 + s(ice_conc,k=4), family = binomial, method = 'REML', 
               weights=weights,data = KGI_binned)
# AIC: 226.7411
# summary
# Formula:
#   Gm ~ productivity_0 + s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -23.92483   17.81480  -1.343    0.179    
# productivity_0  -0.17054    0.02355  -7.241 4.46e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.894  1.993  2.425   0.279
# 
# R-sq.(adj) =  0.465   Deviance explained = 45.2%
# -REML = 113.53  Scale est. = 1         n = 181

# Now adding in more significant smooth terms
KGI_gam <- gam(Gm ~ productivity_0 + s(ice_conc,k=4) + s(o2_0,k=4), family = binomial, method = 'REML', 
               weights=weights,data = KGI_binned)
# AIC: 227.3305
# summary:
# Formula:
#   Gm ~ productivity_0 + s(ice_conc, k = 4) + s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -23.02691   17.50413  -1.316    0.188    
# productivity_0  -0.15742    0.02759  -5.705 1.16e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.878  1.988  2.579   0.249
# s(o2_0)     1.862  2.262  3.283   0.231
# 
# R-sq.(adj) =  0.474   Deviance explained = 46.2%
# -REML = 113.07  Scale est. = 1         n = 181
# Removing sea ice concentration (will add back at end) to see if it is confounding effect of other smooths
KGI_gam <- gam(Gm ~ productivity_0 + s(o2_0,k=4), family = binomial, method = 'REML', 
               weights=weights,data = KGI_binned)
# summary:
# Formula:
#   Gm ~ productivity_0 + s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     0.17990    0.22454   0.801    0.423    
# productivity_0 -0.09595    0.02157  -4.449 8.63e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(o2_0) 2.501  2.828  29.45 2.11e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.239   Deviance explained = 21.4%
# -REML = 164.13  Scale est. = 1         n = 181

# ice concentration is important, deviance explained dropped by around ~20% without it
# Continuing to build GAM with significant terms from earlier, will add sea ice back at end
# Adding salinity now
KGI_gam <- gam(Gm ~ productivity_0 + s(o2_0,k=4) + s(salinity_0,k=4), family = binomial, method = 'REML', 
               weights=weights,data = KGI_binned)
# AIC: 230.8216
# summary:
# Formula:
#   Gm ~ productivity_0 + s(o2_0, k = 4) + s(salinity_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -0.62127    0.49547  -1.254     0.21    
# productivity_0 -0.17156    0.02319  -7.399 1.37e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(o2_0)       1.000  1.000   7.41  0.00649 ** 
#   s(salinity_0) 2.761  2.952  26.05 1.02e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.465   Deviance explained = 45.2%
# -REML = 117.45  Scale est. = 1         n = 181
# .............................................
# None of concurvity or gam.check outputs were concerning
# Even though AIC slightly increased, keeping it in the model for now
#    - It has a significant effect & increases % deviance explained


# Next step: adding in julian day
KGI_gam <- gam(Gm ~ productivity_0 + s(o2_0,k=4) + s(salinity_0,k=4) + s(julian_day,k=4), 
               family = binomial, method = 'REML', weights=weights,data = KGI_binned)
# AIC: 227.8585
# summary:
# Formula:
#   Gm ~ productivity_0 + s(o2_0, k = 4) + s(salinity_0, k = 4) + 
#   s(julian_day, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -0.95468    0.57359  -1.664    0.096 .  
# productivity_0 -0.17739    0.02394  -7.410 1.26e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(o2_0)       1.000  1.000 11.539 0.000682 ***
#   s(salinity_0) 2.805  2.967 18.385 0.000393 ***
#   s(julian_day) 1.000  1.000  4.098 0.042946 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.473   Deviance explained = 46.4%
# -REML =  115.3  Scale est. = 1         n = 181
# .....................................................
# Salinity and julian day have high concurvity (0.91), will deal with this after adding back sea ice

# Adding sea ice concentration back in
KGI_gam <- gam(Gm ~ productivity_0 + s(o2_0,k=4) + s(salinity_0,k=4) + s(julian_day,k=4) + s(ice_conc,k=4), 
               family = binomial, method = 'REML', weights=weights,data = KGI_binned)
# AIC: 225.6064
# summary:
# Formula:
#   Gm ~ productivity_0 + s(o2_0, k = 4) + s(salinity_0, k = 4) + 
#   s(julian_day, k = 4) + s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -17.60052   14.15379  -1.244    0.214    
# productivity_0  -0.17401    0.02544  -6.841 7.87e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(o2_0)       1.000  1.000  6.783 0.00921 **
#   s(salinity_0) 1.000  1.000  1.449 0.22879   
# s(julian_day) 1.000  1.000  2.734 0.09826 . 
# s(ice_conc)   1.853  1.982  2.033 0.32873   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.475   Deviance explained =   47%
# -REML = 110.44  Scale est. = 1         n = 181
# ....................................................
# Most problematic concurvity is between ice concentration and salinity, going to remove salinity

# Next step: removing salinity
KGI_gam <- gam(Gm ~ productivity_0 + s(o2_0,k=4) + s(julian_day,k=4) + s(ice_conc,k=4), 
               family = binomial, method = 'REML', weights=weights,data = KGI_binned)
# AIC: 224.9739
# summary:
# Formula:
#   Gm ~ productivity_0 + s(o2_0, k = 4) + s(julian_day, k = 4) + 
#   s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -18.57047   14.65539  -1.267    0.205    
# productivity_0  -0.18234    0.02497  -7.303 2.81e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(o2_0)       1.000  1.000  6.032  0.0140 *
#   s(julian_day) 1.000  1.000  5.147  0.0233 *
#   s(ice_conc)   1.856  1.983  3.180  0.1728  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.477   Deviance explained = 46.7%
# -REML = 111.21  Scale est. = 1         n = 181
# ............................................
# Some issues with julian day & sea ice concentration concurvity, but these are just in worst concurvity
#     - ignoring for now because observed and estimate concurvity are very reasonable (ranging from 0.15 to 0.25)
# However, gam.check shows significant values for all three variables (but especially for julian day)
#     - going to increase k for julian day

# Next step: trying different k for julian day
KGI_gam <- gam(Gm ~ productivity_0 + s(o2_0,k=4) + s(julian_day,k=8) + s(ice_conc,k=4), 
               family = binomial, method = 'REML', weights=weights,data = KGI_binned)
# gam.check:
# Method: REML   Optimizer: outer newton
# full convergence after 9 iterations.
# Gradient range [-6.806735e-06,2.415695e-06]
# (score 111.207 & scale 1).
# Hessian positive definite, eigenvalue range [6.638096e-06,0.3405381].
# Model rank =  15 / 15 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                 k'  edf k-index p-value    
# s(o2_0)       3.00 1.00    0.80   0.005 ** 
#   s(julian_day) 7.00 1.00    0.72  <2e-16 ***
#   s(ice_conc)   3.00 1.86    0.78   0.010 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# .................................................
# gam.check output above stayed the same after doubling k (as did AIC and summary)
# likely culprit: missing predictors
#   - hypothesis: variables at depth better explain species presence (surface is limitation for this model)
# plots imply we should return to k=4 for julian day

# final model:
KGI_final <- gam(Gm ~ productivity_0 + s(o2_0,k=4) + s(julian_day,k=4) + s(ice_conc,k=4), 
               family = binomial, method = 'REML', weights=weights,data = KGI_binned)
# sea ice concentration still not significant, but its p-value is now down to 0.17
#    - if i remove sea ice concentration, deviance explain goes own to 42.3% and AIC increases to 244.5775
# keeping weights in as well, without them p-values of all predictors skyrocket
# my guess: this model will become significant for sea ice concentration 
#   once we figure out why its standard error is so high
#     - at that point, will further develop model to see if other predictors are also significant

# Next strategy: mess with smoothing parameter of final model
# Final model smoothing parameters were ~3000 for oxygen and julian day, ~10^-5 for ice concentration
# Starting with sp of 0.1 for all
KGI_gam <- gam(Gm ~ productivity_0 + s(o2_0,k=4,sp=0.1) + s(julian_day,k=4,sp=0.1) + s(ice_conc,k=4,sp=0.1), 
  family = binomial, weights=weights,data = KGI_binned)
# AIC: 232.7961
# summary:
# Formula:
#   Gm ~ productivity_0 + s(o2_0, k = 4, sp = 0.1) + s(julian_day, 
#                                                      k = 4, sp = 0.1) + s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -0.85940    0.73009  -1.177    0.239    
# productivity_0 -0.20379    0.04505  -4.524 6.08e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(o2_0)       1.440  1.756 13.462 0.00929 **
#   s(julian_day) 1.255  1.453  4.758 0.04443 * 
#   s(ice_conc)   1.064  1.124  8.886 0.00204 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.46   Deviance explained = 44.6%
# UBRE = 0.28617  Scale est. = 1         n = 181
# .................................................
# Model now looks a lot more reasonable
# Starting from beginning using manual sp: multiple model with all variables, find significant ones, etc.

# Single models, testing each parameter for significance
KGI_gam <- gam(Gm ~ s(FSLE,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Formula:
#   Gm ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1771     0.1207  -1.468    0.142
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(FSLE) 2.163  2.547  4.048   0.104
# 
# R-sq.(adj) =  0.0223   Deviance explained = 2.85%
# UBRE = 1.1786  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(SSH,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Formula:
#   Gm ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.1727     0.2871  -4.084 4.43e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(SSH) 1.553  1.884   54.1  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.248   Deviance explained = 25.2%
# UBRE = 0.67826  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(mixed_layer,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Formula:
#   Gm ~ s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -0.226      0.125  -1.808   0.0706 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(mixed_layer) 1.927   2.34  8.938 0.00944 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0463   Deviance explained = 4.92%
# UBRE = 1.1303  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.1210     0.3485  -3.216   0.0013 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(ice_conc) 1.383   1.64   15.7 0.00026 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.162   Deviance explained =   18%
# UBRE = 0.83502  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(ice_diff,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Formula:
#   Gm ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1805     0.1217  -1.483    0.138
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(ice_diff) 1.743   2.09  7.046  0.0384 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0202   Deviance explained = 2.51%
# UBRE = 1.1816  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(salinity_0,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Formula:
#   Gm ~ s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.9088     0.2154  -4.219 2.46e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(salinity_0) 1.945  2.344  46.01  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.221   Deviance explained = 22.9%
# UBRE = 0.73389  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(EKE_0,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1480     0.1192  -1.242    0.214
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.788  2.181  2.064    0.37
# 
# R-sq.(adj) =  5.62e-05   Deviance explained = 0.877%
# UBRE =  1.218  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(o2_0,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Formula:
#   Gm ~ s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.4206     0.1386  -3.034  0.00242 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(o2_0) 1.943  2.343  37.33  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.165   Deviance explained = 14.7%
# UBRE = 0.91437  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(productivity_0,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Formula:
#   Gm ~ s(productivity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.4397     0.1486   -2.96  0.00308 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(productivity_0) 1.791  2.177  36.63  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.125   Deviance explained = 11.9%
# UBRE = 0.97455  Scale est. = 1         n = 181

KGI_gam <- gam(Gm ~ s(julian_day,k=4,sp=0.1),family = binomial, weights=weights,data = KGI_binned) 
# Formula:
#   Gm ~ s(julian_day, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.7167     0.1766  -4.058 4.95e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(julian_day) 2.089  2.473  49.65  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.296   Deviance explained = 24.8%
# UBRE = 0.69251  Scale est. = 1         n = 181

# Significant predictors are: SSH, mixed layer depth, ice concentration, ice concentration difference,
# salinity, oxygen, productivity, julian day (all with edf greater than 1, no linear relationships)
# Building final GAM one predictor at a time as before, checking for concurvity each time

# Starting with sea ice concentration and mixed_layer
KGI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(mixed_layer,sp=0.1,k=4),
               family = binomial, weights=weights,data = KGI_binned) 
# AIC: 304.3143
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1) + s(mixed_layer, sp = 0.1, 
#                                         k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.7568     0.5503  -3.192  0.00141 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)    1.161  1.297  12.83 0.000432 ***
#   s(mixed_layer) 1.748  2.151  22.60 1.12e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.241   Deviance explained = 25.8%
# UBRE = 0.68129  Scale est. = 1         n = 181
# .........................................
# concurvity is not problematic, keeping both variables in

# Adding difference in ice concentration
KGI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(mixed_layer,sp=0.1,k=4) + s(ice_diff,sp=0.1,k=4),
               family = binomial, weights=weights,data = KGI_binned)
# AIC: 305.2008
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1) + s(mixed_layer, sp = 0.1, 
#                                         k = 4) + s(ice_diff, sp = 0.1, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.7799     0.5608  -3.174   0.0015 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)    1.151  1.278 12.523 0.000504 ***
#   s(mixed_layer) 1.715  2.113 19.554 4.53e-05 ***
#   s(ice_diff)    1.291  1.505  0.725 0.470556    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.238   Deviance explained = 26.2%
# UBRE = 0.68619  Scale est. = 1         n = 181
# ................................................
# concurvity does not look problematic, however dropping ice difference because not significant (and AIC went up)

# dropped ice difference, adding NPP
KGI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(mixed_layer,sp=0.1,k=4) + s(productivity_0,sp=0.1,k=4),
               family = binomial, weights=weights,data = KGI_binned)
# AIC: 232.2211
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1) + s(mixed_layer, sp = 0.1, 
#                                         k = 4) + s(productivity_0, sp = 0.1, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.2134     0.6273  -3.529 0.000418 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(ice_conc)       1.125  1.233 10.524 0.00123 ** 
#   s(mixed_layer)    1.617  1.981  0.631 0.72160    
# s(productivity_0) 1.622  1.979 50.772 < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.465   Deviance explained = 44.5%
# UBRE = 0.28299  Scale est. = 1         n = 181
# ...........................................
# concurvity doesn't look problematic, but not keeping mixed layer depth (since it is no longer significant)

# Dropping mixed layer depth
KGI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(productivity_0,sp=0.1,k=4),
               family = binomial, weights=weights,data = KGI_binned)
# AIC: 231.6787
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1) + s(productivity_0, sp = 0.1, 
#                                         k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -2.230      0.578  -3.858 0.000114 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       1.140  1.259  13.18 0.000314 ***
#   s(productivity_0) 1.634  1.998  67.76  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.461   Deviance explained = 43.9%
# UBRE = 0.27999  Scale est. = 1         n = 181

# Adding in oxygen concentration
KGI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(productivity_0,sp=0.1,k=4) + s(o2_0,k=4,sp=0.1),
               family = binomial, weights=weights,data = KGI_binned)
# AIC: 233.5662
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1) + s(productivity_0, sp = 0.1, 
#                                         k = 4) + s(o2_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.3996     0.6359  -3.774 0.000161 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       1.094  1.179 11.970 0.000468 ***
#   s(productivity_0) 1.466  1.767 49.361  < 2e-16 ***
#   s(o2_0)           1.554  1.885  1.999 0.294160    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.461   Deviance explained = 44.1%
# UBRE = 0.29042  Scale est. = 1         n = 181
# ......................................................
# oxygen is not significant, scrapping it for next model

# Replacing oxygen concentration with julian day
KGI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(productivity_0,sp=0.1,k=4) + s(julian_day,k=4,sp=0.1),
               family = binomial, weights=weights,data = KGI_binned)
# AIC: 232.8101
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1) + s(productivity_0, sp = 0.1, 
#                                         k = 4) + s(julian_day, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.2568     0.5834  -3.868  0.00011 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       1.135  1.250 13.062 0.000327 ***
#   s(productivity_0) 1.578  1.913 21.710 1.48e-05 ***
#   s(julian_day)     1.297  1.523  0.026 0.962129    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.459   Deviance explained = 44.2%
# UBRE = 0.28624  Scale est. = 1         n = 181
# ....................................................
# julian day and ice concentration have problematic concurvity
# julian day not significant, removing it

# Replacing julian day with salinity
KGI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(productivity_0,sp=0.1,k=4) + s(salinity_0,k=4,sp=0.1),
               family = binomial, weights=weights,data = KGI_binned)
# AIC: 232.1874
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1) + s(productivity_0, sp = 0.1, 
#                                         k = 4) + s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.1706     0.5907  -3.675 0.000238 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(ice_conc)       1.142  1.261  6.764 0.00962 ** 
#   s(productivity_0) 1.615  1.969 55.141 < 2e-16 ***
#   s(salinity_0)     1.655  1.997  1.641 0.43257    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.464   Deviance explained = 44.6%
# UBRE = 0.2828  Scale est. = 1         n = 181
# ..................................................
# salinity and ice concentration have problematic concurvity
# salinity not significant, removing it

# Replacing salinity with sea surface height
KGI_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(productivity_0,sp=0.1,k=4) + s(SSH,k=4,sp=0.1),
               family = binomial, weights=weights,data = KGI_binned)
# AIC: 231.6068
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1) + s(productivity_0, sp = 0.1, 
#                                         k = 4) + s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.4301     0.6616  -3.673  0.00024 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(ice_conc)       1.104  1.195  10.08 0.00127 ** 
#   s(productivity_0) 1.618  1.978  48.32 < 2e-16 ***
#   s(SSH)            1.379  1.646   1.13 0.33579    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.463   Deviance explained = 44.6%
# UBRE = 0.2796  Scale est. = 1         n = 181
# ...............................................
# SSH not significant, dropping it from final model

# final KGI model: ice concentration and net primary production are significant (43.9% deviance explained)
KGI_final <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(productivity_0,sp=0.1,k=4),
               family = binomial, weights=weights,data = KGI_binned)
# AIC: 231.6787
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4, sp = 0.1) + s(productivity_0, sp = 0.1, 
#                                         k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -2.230      0.578  -3.858 0.000114 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       1.140  1.259  13.18 0.000314 ***
#   s(productivity_0) 1.634  1.998  67.76  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.461   Deviance explained = 43.9%
# UBRE = 0.27999  Scale est. = 1         n = 181

# Scrapped strategies below:
# Trying 4 day bins, similar variable-by-variable initial approach
KGI_gam <- gam(Gm ~ s(ice_conc, k=4), family=binomial, method='REML',data=KGI_4day)
# Bins make the data coarser and weren't doing anything to make the model look more reasonable
#   - so, we are sticking with ACF binned data


# --------------------- Step 5c: Clarence Island GAM ------------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Similar to KGI: weighing to reach roughly 1:1 1s to 0s and setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in CI_pred, include FSLE, SSH, mixed layer, ice concentration, 
#   ice concentration difference, EKE, oxygen concentration, net primary production

# adding weights for CI
# weighing presence (1s) at 10 because there are roughly 10x more 0s than 1s in the CI dataframe
CI_binned$weights <- ifelse(CI_binned$Gm == 1, 10, 1)

# Starting with FSLE and FSLE orientation interaction
CI_gam <- gam(Gm ~ s(FSLE,fsle_orient,k=4),
              family=binomial,weights=weights,data=CI_binned, method = 'REML')
# summary:
# Formula:
#   Gm ~ s(FSLE, fsle_orient, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.06481    0.12253  -0.529    0.597
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(FSLE,fsle_orient) 2.364  2.596  9.685  0.0155 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0131   Deviance explained = 2.75%
# -REML = 193.36  Scale est. = 1         n = 155

# Sea surface height
CI_gam <- gam(Gm ~ s(SSH,k=4),
              method='REML',family=binomial,weights=weights,data=CI_binned)
# summary:
# Formula:
#   Gm ~ s(SSH, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1631     0.1308  -1.247    0.213
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(SSH)   2  2.397  20.32 8.43e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0529   Deviance explained = 6.22%
# -REML = 186.21  Scale est. = 1         n = 155

# mixed layer depth
CI_gam <- gam(Gm ~ s(mixed_layer,k=4),
              method='REML',family=binomial,weights=weights,data=CI_binned)
# summary:
# Formula:
#   Gm ~ s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.3862     0.1439  -2.684  0.00727 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(mixed_layer) 2.928  2.996  47.13  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.149   Deviance explained = 15.9%
# -REML = 170.78  Scale est. = 1         n = 155

# sea ice concentration
CI_gam <- gam(Gm ~ s(ice_conc,k=4),
              method='REML',family=binomial,weights=weights,data=CI_binned)
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1810     0.1442  -1.255    0.209
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc) 2.444  2.785  11.28 0.00665 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0299   Deviance explained = 5.15%
# -REML = 188.88  Scale est. = 1         n = 155

# ice regime
CI_gam <- gam(Gm ~ ice_regime,
              method='REML', family=binomial,weights=weights,data=CI_binned)
# summary:
# Formula:
#   Gm ~ ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            0.7340     0.2483   2.956  0.00312 ** 
#   ice_regimeincreasing  -0.8737     0.3939  -2.218  0.02654 *  
#   ice_regimenone        -1.8491     0.3579  -5.167 2.38e-07 ***
#   ice_regimestable      -0.3185     0.3346  -0.952  0.34125    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# R-sq.(adj) =  0.0731   Deviance explained = 8.93%
# -REML = 179.13  Scale est. = 1         n = 155

# eddy kinetic energy
CI_gam <- gam(Gm ~ s(EKE_0,k=4),
              method='REML',family=binomial,weights=weights,data=CI_binned)
# summary:
# Formula:
#   Gm ~ s(EKE_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.09444    0.12851  -0.735    0.462
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(EKE_0) 2.509  2.798  6.985  0.0566 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0125   Deviance explained = 2.73%
# -REML =  194.4  Scale est. = 1         n = 155

# eddy kinetic energy variation
CI_gam <- gam(Gm ~ s(EKE_0_mad,k=4),
              method='REML',family=binomial,weights=weights,data=CI_binned)
# summary:
# Formula:
#   Gm ~ s(EKE_0_mad, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.0374     0.1208   -0.31    0.757
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0_mad) 1.908  2.298  5.443   0.113
# 
# R-sq.(adj) =  0.00613   Deviance explained = 1.68%
# -REML = 195.12  Scale est. = 1         n = 155

# oxygen concentration
CI_gam <- gam(Gm ~ s(o2_0,k=4),
              method='REML',family=binomial,weights=weights,data=CI_binned)
# summary:
# Formula:
#   Gm ~ s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -4.776      1.102  -4.333 1.47e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(o2_0) 2.905  2.991  34.02  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.204   Deviance explained = 23.1%
# -REML = 156.86  Scale est. = 1         n = 155

# net primary production
CI_gam <- gam(Gm ~ s(productivity_0,k=4),
              method='REML',family=binomial,weights=weights,data=CI_binned)
# Formula:
#   Gm ~ s(productivity_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -3.229      1.225  -2.636  0.00839 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(productivity_0) 2.801  2.962  20.52 0.000122 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.133   Deviance explained =   17%
# -REML = 167.57  Scale est. = 1         n = 155

# significant variables: FSLE (with interaction), SSH, mixed layer, ice concentration, 
#    ice regime, oxygen, primary production

# MAKING MULTIVARIABLE GAMs
# Starting model with sea ice concentration and FSLE
CI_gam <- gam(Gm ~ s(ice_conc,k=4) + s(FSLE,fsle_orient,k=4), 
              method='REML',family=binomial, weights=weights, data=CI_binned)
# AIC: 369.6846
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(FSLE, fsle_orient, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.2431     0.1454  -1.672   0.0946 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc)         2.528  2.843  17.33 0.00238 **
#   s(FSLE,fsle_orient) 2.538  2.787   9.16 0.01413 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0476   Deviance explained =  8.5%
# -REML = 185.27  Scale est. = 1         n = 155
# ......................................................
# concurvity, plots, and gam.check all look good, moving on

# adding SSH
CI_gam <- gam(Gm ~ s(ice_conc,k=4) + s(FSLE,fsle_orient,k=4) + s(SSH,k=4), 
              method='REML',family=binomial, weights=weights, data=CI_binned)
# AIC: 310.3127
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(FSLE, fsle_orient, k = 4) + s(SSH, 
#                                                             k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.8014     0.2046  -3.918 8.94e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)         2.610  2.889  16.57  0.00106 ** 
#   s(FSLE,fsle_orient) 2.000  2.000  30.77 1.54e-06 ***
#   s(SSH)              2.722  2.940  38.69  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.227   Deviance explained = 24.9%
# -REML = 156.14  Scale est. = 1         n = 155
# ....................................................
# all aspects of model look good, moving on

# adding mixed layer depth
CI_gam <- gam(Gm ~ s(ice_conc,k=4) + s(FSLE,fsle_orient,k=4) + s(SSH,k=4) + s(mixed_layer,k=4), 
              method='REML',family=binomial, weights=weights, data=CI_binned)
# AIC: 296.5749
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(FSLE, fsle_orient, k = 4) + s(SSH, 
#                                                             k = 4) + s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.9183     0.2109  -4.354 1.33e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)         2.653  2.907   9.31 0.014537 *  
#   s(FSLE,fsle_orient) 2.000  2.000  22.77  1.1e-05 ***
#   s(SSH)              2.471  2.800  14.51 0.002152 ** 
#   s(mixed_layer)      2.873  2.987  17.46 0.000997 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.299   Deviance explained = 29.9%
# -REML = 150.41  Scale est. = 1         n = 155
# ...................................................
# concurvity and plots look good

# adding oxygen concentration
CI_gam <- gam(Gm ~ s(ice_conc,k=4) + s(FSLE,fsle_orient,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) + 
                s(o2_0,k=4), 
              method='REML',family=binomial, weights=weights, data=CI_binned)
# AIC: 227.2797
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(FSLE, fsle_orient, k = 4) + s(SSH, 
#                                                             k = 4) + s(mixed_layer, k = 4) + s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -3.846      1.171  -3.284  0.00103 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)         2.815  2.971  11.30  0.00948 ** 
#   s(FSLE,fsle_orient) 2.922  2.989  23.25 4.84e-05 ***
#   s(SSH)              2.165  2.584  22.00 7.68e-05 ***
#   s(mixed_layer)      2.813  2.970  11.35  0.00618 ** 
#   s(o2_0)             2.701  2.928  12.86  0.00402 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.489   Deviance explained = 49.6%
# -REML = 115.48  Scale est. = 1         n = 155
# ..........................................................
# concurvity potential problematic later on, plots looking good, moving on

# adding primary production
CI_gam <- gam(Gm ~ s(ice_conc,k=4) + s(FSLE,fsle_orient,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) + 
                s(o2_0,k=4) + s(productivity_0,k=4), 
              method='REML',family=binomial, weights=weights, data=CI_binned)
# AIC: 210.9987
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(FSLE, fsle_orient, k = 4) + s(SSH, 
#                                                             k = 4) + s(mixed_layer, k = 4) + s(o2_0, k = 4) + s(productivity_0, 
#                                                                                                                 k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -3.840      1.062  -3.616 0.000299 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)         2.846  2.980 17.816 0.000397 ***
#   s(FSLE,fsle_orient) 2.938  2.993 22.951 2.89e-05 ***
#   s(SSH)              2.597  2.891 25.644 2.00e-05 ***
#   s(mixed_layer)      2.883  2.985 17.920 0.000356 ***
#   s(o2_0)             2.110  2.470  7.141 0.065884 .  
# s(productivity_0)   1.000  1.001 12.002 0.000534 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.569   Deviance explained = 54.2%
# -REML = 105.95  Scale est. = 1         n = 155
# ...................................................
# problematic concurvity for o2, no longer significant

# dropping oxygen concentration
CI_gam <- gam(Gm ~ s(ice_conc,k=4) + s(FSLE,fsle_orient,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) + 
                s(productivity_0,k=4), 
              method='REML',family=binomial, weights=weights, data=CI_binned)
# AIC: 222.5486
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(FSLE, fsle_orient, k = 4) + s(SSH, 
#                                                             k = 4) + s(mixed_layer, k = 4) + s(productivity_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.3856     0.4185  -5.701 1.19e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)         2.794  2.968  12.17  0.00469 ** 
#   s(FSLE,fsle_orient) 2.925  2.992  31.89  < 2e-16 ***
#   s(SSH)              2.801  2.970  35.18  < 2e-16 ***
#   s(mixed_layer)      2.882  2.988  21.09 8.22e-05 ***
#   s(productivity_0)   1.000  1.000  33.40  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.535   Deviance explained =   50%
# -REML = 113.88  Scale est. = 1         n = 155
# ....................................................
# concurvity resolved, plots looking better

# adding ice regime
CI_gam <- gam(Gm ~ s(ice_conc,k=4) + s(FSLE,fsle_orient,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) + 
                s(productivity_0,k=4) + ice_regime, 
              method='REML',family=binomial, weights=weights, data=CI_binned)
# AIC: 215.2103
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(FSLE, fsle_orient, k = 4) + s(SSH, 
#                                                             k = 4) + s(mixed_layer, k = 4) + s(productivity_0, k = 4) + 
#   ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)          -3.56635    1.57252  -2.268  0.02333 * 
#   ice_regimeincreasing -2.98784    1.03944  -2.874  0.00405 **
#   ice_regimenone       -2.16660    1.55629  -1.392  0.16387   
# ice_regimestable      0.03773    0.63093   0.060  0.95232   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)         2.795  2.968  11.00 0.007729 ** 
#   s(FSLE,fsle_orient) 2.925  2.991  29.02 3.03e-06 ***
#   s(SSH)              2.822  2.976  29.93 2.98e-06 ***
#   s(mixed_layer)      2.846  2.982  20.88 0.000187 ***
#   s(productivity_0)   2.625  2.864  11.52 0.010439 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.575   Deviance explained = 54.4%
# -REML = 105.94  Scale est. = 1         n = 155
# ....................................................
# concurvity went up, plot error bars drastically increased
# since ice regime is not very significant, dropping it and reverting to model at line 1932

# final variables: ice concentration, FSLE, sea surface height, mixed layer depth, primary production
CI_final <- gam(Gm ~ s(ice_conc,k=4) + s(FSLE,fsle_orient,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) + 
                s(productivity_0,k=4), 
              method='REML',family=binomial, weights=weights, data=CI_binned)
# running gam.check
# no p-values significant in gam.check not messing with smooothing


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
    plot_annotation(title = paste0("Long-finned Pilot Whale at ",
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
KGI_pred <- c('ice_conc','productivity_0')
KGI_plots <- visualizeGAM(KGI_final, KGI_pred, 'KGI')

EI_pred <- c('mixed_layer')
EI_plots <- visualizeGAM(EI_final, EI_pred, 'EI')


CI_pred <- c('ice_conc','SSH','mixed_layer','productivity_0')
CI_plots <- visualizeGAM(CI_final, CI_pred, 'CI')
# adding plot for FSLE interaction
vis.gam(CI_final, view = c('FSLE','fsle_orient'), plot.type = 'contour', too.far = 0.1,
        color='gray', type = 'response',xlab="FSLE magnitude",ylab='FSLE Orientation (°)',
        main = 'Probability of Long-Finned Pilot Whale Presence (color)')