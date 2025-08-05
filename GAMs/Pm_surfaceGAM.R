library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)
# ------------- Step 0: Choose Species ----------------
# Modeling Pm for all sites it is present (EI & CI), 40 km radius environmental data
species <- c('Pm') # options: BW29, BW37, Oo, Pm, Gm
# BW29 = Southern bottlenose whale, BW37 = Gray's and strap-toothed whales
# Oo = Killer whale, Pm = Sperm Whale, Gm = Long-finned pilot whale
# Note: not enough data to model BW58 at any site
sites <- c('EI','CI')

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
binned_plot <- grid.arrange(binnedTimeseries(EI_binned,'EI',EI_acf), 
                            binnedTimeseries(CI_binned,'CI',CI_acf), nrow=2, 
                            top = paste('ACF Binned Species Presence for ', name(species), sep=''))

# -------------- Step 4: VIF for Correlation -------------------
# Not including any variables at depth for initial predictors.
# Also, not modeling AAO (varies on yearly timescales) and ice thickness
# Not including julian day since EI and CI cover less than a year

# ELEPHANT ISLAND
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0          EKE_0 
# 2.360611       3.206931      10.759774       4.394141       2.224918       3.001744      42.743276       2.760488 
# o2_0         chla_0 productivity_0 
# 36.578570      20.222276      13.903403 

# removing surface temperature
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0           o2_0 
# 2.415926       3.026538      10.946946       3.022225       2.247418       3.036003       2.393427      20.531634 
# chla_0 productivity_0 
# 15.464146      14.221125 

# removing oxygen
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0",'chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0         chla_0 
# 2.135270       2.556825       7.248250       2.138651       2.073489       2.879637       1.437701      15.658254 
# productivity_0 
# 8.190306

# removing chlorophyll
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0",'productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0          EKE_0 productivity_0 
# 2.255186       1.899690       5.259062       2.041016       1.828531       2.375809       1.381850       3.743075 

# removing mixed layer depth
EI_pred <- c("FSLE", "SSH", "ice_conc", 'ice_diff', "salinity_0", 
             "EKE_0",'productivity_0')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE            SSH       ice_conc       ice_diff     salinity_0          EKE_0 productivity_0 
# 1.224019       1.925341       2.126818       1.373726       1.570107       1.326052       2.474364 

# final predictors for EI: FSLE, SSH, ice concentration, difference in ice concentration, salinity,
#   eddy kinetic energy, net primary production

# CLARENCE ISLAND
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0          EKE_0 
# 2.536083       5.642813       9.765472       6.340287       2.049896      30.013032      11.797766       4.073984 
# o2_0         chla_0 productivity_0 
# 105.749775      60.508001      23.655332 

# removing oxygen
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'chla_0','productivity_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc       ice_diff     salinity_0  temperature_0          EKE_0 
# 6.124522       4.719274       5.171847       2.376885       1.512895       7.474652       6.562428       1.539592 
# chla_0 productivity_0 
# 7.584033      11.922138  

# removing primary production
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', "salinity_0", 
             "temperature_0", "EKE_0",'chla_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE           SSH   mixed_layer      ice_conc      ice_diff    salinity_0 temperature_0         EKE_0 
# 5.534956      4.665268      4.453090      2.386922      1.478937      6.773729      4.229930      1.515038 
# chla_0 
# 1.655864 

# removing salinity
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff', 
             "temperature_0", "EKE_0",'chla_0')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE           SSH   mixed_layer      ice_conc      ice_diff temperature_0         EKE_0        chla_0 
# 4.601692      2.635957      2.717384      2.366324      1.361689      4.756979      1.331402      1.587046 

# final predictors for CI: FSLE, SSH, mixed layer depth, ice concentration, ice concentration difference,
#   temperature, EKE, chlorophyll

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
# NOT weighing for sperm whales since majority of bins have presence
# List of predictors are stored in EI_pred include FSLE, SSH, ice concentration, 
#   difference in ice concentration, salinity, eddy kinetic energy, net primary production

# starting with FSLE
EI_gam <- gam(Pm ~ s(FSLE,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Pm ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   0.5662     0.3193   1.773   0.0761 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(FSLE) 1.422  1.716  14.47 0.000367 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.398   Deviance explained =   35%
# UBRE = -0.04596  Scale est. = 1         n = 69

# sea surface height
EI_gam <- gam(Pm ~ s(SSH,k=4,sp=0.1), family=binomial, data=EI_binned)
# Formula:
#   Pm ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3267     0.2456    1.33    0.183
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(SSH) 1.74  2.124  0.962    0.67
# 
# R-sq.(adj) =  -0.00924   Deviance explained = 1.29%
# UBRE = 0.42263  Scale est. = 1         n = 69

# sea ice concentration
EI_gam <- gam(Pm ~ s(ice_conc,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Pm ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.4050     0.2711   1.494    0.135
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(ice_conc) 1.264  1.477  3.394  0.0758 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0842   Deviance explained = 7.67%
# UBRE = 0.32206  Scale est. = 1         n = 69

# ice concentration difference
EI_gam <- gam(Pm ~ s(ice_diff,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Pm ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3276     0.2452   1.336    0.182
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_diff) 1.433  1.734  0.744   0.712
# 
# R-sq.(adj) =  -0.006   Deviance explained =  1.3%
# UBRE = 0.41358  Scale est. = 1         n = 69

# salinity
EI_gam <- gam(Pm ~ s(salinity_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Pm ~ s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.4087     0.3158   1.294    0.196
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(salinity_0) 1.195  1.361  15.54 0.00023 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.325   Deviance explained = 29.3%
# UBRE = 0.026351  Scale est. = 1         n = 69

# eddy kinetic energy
EI_gam <- gam(Pm ~ s(EKE_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Pm ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3530     0.2688   1.313    0.189
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(EKE_0) 1.532  1.862  9.051 0.00519 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.218   Deviance explained = 18.3%
# UBRE = 0.18573  Scale est. = 1         n = 69

# net primary production
EI_gam <- gam(Pm ~ s(productivity_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# summary:
# Formula:
#   Pm ~ s(productivity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1040     0.4063   0.256    0.798
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(productivity_0) 1.2  1.369  17.97 6.35e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.553   Deviance explained = 45.7%
# UBRE = -0.19692  Scale est. = 1         n = 69

# significant variables: FSLE, salinity, EKE, net primary production
# adding these one by one to make final model

# starting with FSLE and salinity
EI_gam <- gam(Pm ~ s(FSLE,k=4,sp=0.1) + s(salinity_0,k=4,sp=0.1), 
              family=binomial, data=EI_binned)
# AIC: 56.98019
# summary:
# Formula:
#   Pm ~ s(FSLE, k = 4, sp = 0.1) + s(salinity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.5833     0.3677   1.586    0.113
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(FSLE)       1.338  1.593  7.529 0.00817 **
#   s(salinity_0) 1.241  1.437  8.026 0.00771 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.518   Deviance explained = 46.9%
# UBRE = -0.1742  Scale est. = 1         n = 69
# ...................................................
# concurvity and plot look good

# adding EKE
EI_gam <- gam(Pm ~ s(FSLE,k=4,sp=0.1) + s(salinity_0,k=4,sp=0.1) + s(EKE_0,k=4,sp=0.1), 
              family=binomial, data=EI_binned)
# AIC: 56.35371
# summary:
# Formula:
#   Pm ~ s(FSLE, k = 4, sp = 0.1) + s(salinity_0, k = 4, sp = 0.1) + 
#   s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.6054     0.3795   1.595    0.111
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(FSLE)       1.330  1.581  6.108  0.0163 *
#   s(salinity_0) 1.228  1.415  5.396  0.0304 *
#   s(EKE_0)      1.461  1.760  1.386  0.3172  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.535   Deviance explained = 50.7%
# UBRE = -0.18328  Scale est. = 1         n = 69
# .................................................
# EKE not significant, AIC largely unchanged, EKE plot not significant
# dropping EKE

# dropping EKE, adding in net primary production
EI_gam <- gam(Pm ~ s(FSLE,k=4,sp=0.1) + s(salinity_0,k=4,sp=0.1) + s(productivity_0,k=4,sp=0.1), 
              family=binomial, data=EI_binned)
# AIC: 47.28446
# summary: 
# Formula:
#   Pm ~ s(FSLE, k = 4, sp = 0.1) + s(salinity_0, k = 4, sp = 0.1) + 
#   s(productivity_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1940     0.4827   0.402    0.688
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(FSLE)           1.235  1.428  1.916 0.14598   
# s(salinity_0)     1.121  1.229  2.386 0.12908   
# s(productivity_0) 1.138  1.259  7.643 0.00716 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.668   Deviance explained = 59.2%
# UBRE = -0.31472  Scale est. = 1         n = 69
# ................................................
# going to run AIC on model with just productivity since this one suggests dropping FSLE and salinity


EI_gam <- gam(Pm ~ s(productivity_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# AIC: 55.41229
# for summary output, see productivity model at line 386. Deviance explained was 45.7%

# AIC is lower for GAM with just productivity than it is for GAM with FSLE and salinity
# Deviance explained is comparable (46.9% for FSLE/salinity model)
# Since AIC is lower, using just productivity in final model (will now check different basis functions and smoothing)

# Increasing k
EI_gam <- gam(Pm ~ s(productivity_0,k=8,sp=0.1), family=binomial, data=EI_binned)
# AIC: 54.4717
# While AIC is slightly lower with 8 basis functions, the plot is not complex enough to justify this change
# Sticking with k=4

# changing smoothing parameter
EI_gam <- gam(Pm ~ s(productivity_0,k=4,sp=0.01), family=binomial, data=EI_binned)
# AIC: 56.65483
EI_gam <- gam(Pm ~ s(productivity_0,k=4,sp=1), family=binomial, data=EI_binned)
# AIC: 55.70425
# Error bars too high if smoothing parameter is decreased
# No significant enough reason to increase smoothing, sticking with sp-0.1

# final model for EI: uses just primary production
EI_final <- gam(Pm ~ s(productivity_0,k=4,sp=0.1), family=binomial, data=EI_binned)
# -------------------- Step 5b: CLarence Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing at 0.1
# NOT weighing for sperm whales since majority of bins have presence
# List of predictors are stored in EI_pred include FSLE, SSH, mixed layer depth, ice concentration, 
#   difference in ice concentration, temperature, eddy kinetic energy, chlorophyll

# starting with FSLE
CI_gam <- gam(Pm ~ s(FSLE,k=4,sp=0.1), family=binomial, data=CI_binned)
# summary:
# Formula:
#   Pm ~ s(FSLE, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.8644     0.3299   5.651  1.6e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(FSLE) 1.473   1.78  4.979   0.126
# 
# R-sq.(adj) =  0.0345   Deviance explained = 8.02%
# UBRE = -0.13827  Scale est. = 1         n = 105

# sea surface height
CI_gam <- gam(Pm ~ s(SSH,k=4,sp=0.1), family=binomial, data=CI_binned)
# summary:
# Formula:
#   Pm ~ s(SSH, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.7978     0.3012   5.968  2.4e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(SSH) 1.494  1.816  3.304   0.106
# 
# R-sq.(adj) =  0.0367   Deviance explained = 7.86%
# UBRE = -0.13647  Scale est. = 1         n = 105

# mixed layer depth
CI_gam <- gam(Pm ~ s(mixed_layer,k=4,sp=0.1), family=binomial, data=CI_binned)
# summary:
# Formula:
#   Pm ~ s(mixed_layer, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.6641     0.2683   6.202 5.56e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(mixed_layer) 1.435  1.733  1.128   0.585
# 
# R-sq.(adj) =  0.00295   Deviance explained = 2.21%
# UBRE = -0.087563  Scale est. = 1         n = 105

# sea ice concentration
CI_gam <- gam(Pm ~ s(ice_conc,k=4,sp=0.1), family=binomial, data=CI_binned)
# Formula:
#   Pm ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.9366     0.3818   5.072 3.93e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.223  1.409  3.167   0.114
# 
# R-sq.(adj) =   0.03   Deviance explained = 6.74%
# UBRE = -0.13168  Scale est. = 1         n = 105

# difference in sea ice concentration
CI_gam <- gam(Pm ~ s(ice_diff,k=4,sp=0.1), family=binomial, data=CI_binned)
# Formula:
#   Pm ~ s(ice_diff, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.6472     0.2653   6.209 5.33e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_diff) 1.341  1.603   0.25   0.861
# 
# R-sq.(adj) =  -0.00787   Deviance explained = 0.648%
# UBRE = -0.075525  Scale est. = 1         n = 105

# temperature
CI_gam <- gam(Pm ~ s(temperature_0,k=4,sp=0.1), family=binomial, data=CI_binned)
# Formula:
#   Pm ~ s(temperature_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.7062     0.2744   6.217 5.06e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(temperature_0) 1.631  1.982  4.319     0.1
# 
# R-sq.(adj) =  0.051   Deviance explained =  7.3%
# UBRE = -0.12887  Scale est. = 1         n = 105

# eddy kinetic energy
CI_gam <- gam(Pm ~ s(EKE_0,k=4,sp=0.1), family=binomial, data=CI_binned)
# Formula:
#   Pm ~ s(EKE_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    1.665      0.270   6.166 7.02e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0) 1.505  1.826  0.665   0.633
# 
# R-sq.(adj) =  -0.00408   Deviance explained = 1.09%
# UBRE = -0.076327  Scale est. = 1         n = 105

# chlorophyll
CI_gam <- gam(Pm ~ s(chla_0,k=4,sp=0.1), family=binomial, data=CI_binned)
# Formula:
#   Pm ~ s(chla_0, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.7876     0.2963   6.033 1.61e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(chla_0) 1.739  2.102  6.631  0.0387 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0729   Deviance explained = 8.84%
# UBRE = -0.14045  Scale est. = 1         n = 105
# .............................
# AIC: 90.252261

# only significant variable was chlorophyll-a concentration
# trying different k and smoothing for chlorophyll to get model with better fit

# increasing k
CI_gam <- gam(Pm ~ s(chla_0,k=8,sp=0.1), family=binomial, data=CI_binned)
# AIC: 88.37996
# Change in plot does not justify higher number of basis functions, leaving k at 4

# increasing sp
CI_gam <- gam(Pm ~ s(chla_0,k=4,sp=1), family=binomial, data=CI_binned)
# AIC: 90.51375
# decreasing sp
CI_gam <- gam(Pm ~ s(chla_0,k=4,sp=0.01), family=binomial, data=CI_binned)
# AIC: 88.6712
# more sansitive trendline is not justified because of the larger error bars, keeping sp at 0.1

# final model for Clarence Island
CI_final <- gam(Pm ~ s(chla_0,k=4,sp=0.1), family=binomial, data=CI_binned)

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
EI_pred <- c('productivity_0')
EI_plots <- visualizeGAM(EI_final, EI_pred, 'EI')


CI_pred <- c('chla_0')
CI_plots <- visualizeGAM(CI_final, CI_pred, 'CI')
