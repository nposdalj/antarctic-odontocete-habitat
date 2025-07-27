library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
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
allData <- read.csv("/Users/trisha/R/antarctic-odontocete-habitat/Data/allData_40km.csv")
allData <- allData %>% subset(select=-X)
allData$date <- as.Date(allData$date, "%Y-%m-%d")
# Filter by species relevant data
if(species =='BW29') {
  depths <- c(0, 768)
  sp_specific <- allData %>% subset(select=-c(BW37,BW58,Oo,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_768,salinity_768,EKE_768,
                    chla_0,o2_0,productivity_0,chla_768,o2_768,productivity_768))
} else if(species =='BW37') {
  depths <- c(0, 67, 920) 
  sp_specific <- allData %>% subset(select=-c(BW29,BW58,Oo,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_67,salinity_67,EKE_67,
                    temperature_920,salinity_920,EKE_920, chla_0,o2_0,productivity_0,chla_67,
                    o2_67,productivity_67, chla_920,o2_920,productivity_920))
} else if(species =='Oo') {
  depths <- c(0, 11, 455) 
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_11,salinity_11,EKE_11,
                    temperature_455,salinity_455,EKE_455, chla_0,o2_0,productivity_0,chla_11,
                    o2_11,productivity_11, chla_455,o2_455,productivity_455))
} else if(species =='Pm') {
  depths <- c(0, 375, 1665)
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Oo,Gm)) %>%
    subset(select=c(date,julian_day,Site,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_375,salinity_375,EKE_375,
                    temperature_1665,salinity_1665,EKE_1665, chla_0,o2_0,productivity_0,chla_375,
                    o2_375,productivity_375, chla_1665,o2_1665,productivity_1665))
} else if(species =='Gm') {
  depths <- c(0, 16, 635) 
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Oo,Pm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_16,salinity_16,EKE_16,
                    temperature_635,salinity_635,EKE_635, chla_0,o2_0,productivity_0,chla_16,
                    o2_16,productivity_16, chla_635,o2_635,productivity_635))
} else
  print('Species code not valid. Check inputs.')


# ------------- Step 2: Average by ACF ------------
acf_table <- read.csv("/Users/trisha/R/antarctic-odontocete-habitat/Autocorrelation/acf_table.csv")
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
                         o2_0 = mean_col('o2_0'), productivity_0 = mean_col('productivity_0'))
  
  # Adding shallow dive depth for applicable variables
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
intl_predictors <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
                     "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(intl_predictors, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula), family = binomial, data = EI_binned)

# KING GEORGE ISLAND
# adding julian day into initial formula
intl_predictors <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
                     "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0', 'julian_day')
mod_formula <- paste(species, "~", paste(intl_predictors, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), family = binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0 
# 1.501086       2.069468       7.185479       3.138806       1.568101       4.337320      35.906240 
# EKE_0           o2_0         chla_0 productivity_0     julian_day 
# 1.074483      13.137926       7.994408      14.278968       3.180173

# removing temp_0
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
                    "EKE_0",'o2_0','chla_0','productivity_0', 'julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), family = binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0 
# 1.454792       2.028160       4.876820       2.908102       1.357751       4.092606       1.072813 
# o2_0         chla_0 productivity_0     julian_day 
# 2.755725       7.688033       6.222557       3.105939 

# removing chlorophyll
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
              "EKE_0",'o2_0','productivity_0', 'julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula), family = binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0 
# 1.434252       1.877656       4.388842       2.821834       1.358126       3.954242       1.064858 
# o2_0 productivity_0     julian_day 
# 2.577644       2.466533       2.559630 

# Final predictors for KGI: FSLE, SSH, mixed layer depth, sea ice concentration, difference
# in ice concentration, salinity, EKE, oxygen, net primary production, julian day

# CLARENCE ISLAND
intl_predictors <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
                     "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(intl_predictors, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula), family = binomial, data = CI_binned)


# -------------- Step 5: Build GAMs ------------------------
# Function to visualize GAMs on a probability scale with the proper confidence interval
# Run this for each iteration of the model to plot smooth terms
# Currently running into issues with plots: they look more reasonable without including model intercept
#      -> remove shift argument and seWithMean argument?
plotGam <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],seWithMean=TRUE))
}
# Run this if all plots should be in one figure
plotGam1 <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],seWithMean=TRUE,pages=1))
}

# ELEPHANT ISLAND

# KING GEORGE ISLAND
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


# CLARENCE ISLAND