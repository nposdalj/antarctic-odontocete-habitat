# Currently working with 40 km radius for all variables EXCEPT FSLE
# Change the loaded file when FSLE is corrected
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
allData <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/allData_40km.csv")
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
# Run this for each iteration of the model
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

# Variables that were significant (p < 0.05): SSH, salinity, oxygen, net primary production, julian day
# Variables that were insignifacnt: FSLE, mixed layer depth, EKE, ice concentration, ice concentration difference
# Still include ice concentration (p=0.15) and/or ice concentration difference (p=0.417) because ecologically important

# Next steps: one by one build final model for binned by ACF (two days)
# Repeat process for other binning lengths (3 and 5 days) to reduce amount of 0s
# Pick best model


# CLARENCE ISLAND