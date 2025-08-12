library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)

# ------------- Step 0: Choose Species ----------------
species <- c('BW29','Gm','BW37') # species with enough data to model all sites
# BW29 = Southern bottlenose whale, BW37 = Gray's and strap-toothed whales
# Oo = Killer whale, Pm = Sperm Whale, Gm = Long-finned pilot whale
# Note: not enough data to model BW58 at any site
sites <- c('EI','KGI','CI') # Modeling across all sites

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
allData <- read.csv("/Users/trisha/scripps/antarctic-odontocete-habitat/Data/allData_40km.csv")
allData <- allData %>% subset(select=-X)
allData$date <- as.Date(allData$date, "%Y-%m-%d")
# Filter by species relevant data
# Only adding standard deviations of surface variables, feel free to change that if needed
depths <- c(0, 768)
surface <- allData %>%
  subset(select=c(date,Site,julian_day,Gm,Pm,Oo,BW29,BW37,BW58,
                  AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                  temperature_0,salinity_0,EKE_0,
                  chla_0,o2_0,productivity_0,
                  ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                  chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime,bathymetry))


# ------------- Step 2: Average by ACF ------------
# Binning each site by its ACF value, then joining the three sites into one dataframe
acf_table <- read.csv("/Users/trisha/scripps/antarctic-odontocete-habitat/Autocorrelation/acf_table.csv")
acfVal <- function(site,sp) {
  row_idx <- which(acf_table$site == site) # Row index for the site
  acf_val <- acf_table[row_idx,sp][[1]]
  return(acf_val)
}

# binning an individual site by its acf value
binByACF <- function(site, bin,sp) {
  # Filtering by site and adding bin_start date
  bin <- ifelse(bin==0, 1, bin) # if the acf value is 0 (indicating no data for that site), bin by 1 instead
  filtered <- surface %>% filter(Site==site) %>% # Filtering by site and setting floor date based on binning
    mutate(bin_start = floor_date(date, unit = paste(bin, 'days')))
  
  # Function for taking mean with correct syntax for a column
  mean_col <- function(col) expr(mean(!!sym(col), na.rm = TRUE))
  # Expression to summarize species presence
  sp_expr <- set_names(list(expr(as.integer(any(!!sym(sp) == 1)))), sp)
  
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
                         salinity_0_sd = mean_col('salinity_sd_0'), bathymetry = mean_col('bathymetry'))
  
  
  # Joining together dataframes and grouping to make final binned data
  all_summaries <- c(sp_expr, summarize_cols)
  binned <- filtered %>%
    group_by(bin_start, Site) %>%
    summarise(!!!all_summaries, .groups = "drop")
  
  return(binned)
}

combineSites <- function(sp) {
  # dataframes for each site for given species
  EI_acf <- acfVal('EI', sp)
  EI_binned <- binByACF('EI',EI_acf, sp)
  KGI_acf <- acfVal('KGI', sp)
  KGI_binned <- binByACF('KGI',KGI_acf, sp)
  CI_acf <- acfVal('CI', sp)
  CI_binned <- binByACF('CI',CI_acf, sp)
  
  # Creating all site, acf binned dataframe
  sp_binned <- rbind(EI_binned, KGI_binned, CI_binned)
  
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

# Running through all given species to generate all site, binned dataframe
for(s in species) {
  assign(x = paste0(s, "_binned"), value = combineSites(s), envir = .GlobalEnv)
}



# ------------- Step 3: Plot Presence Timeseries --------------
binnedTimeseries <- function(data,sp) { # Function to create a timeseries plot
  # Making timeseries 
  ggplot(data = data, mapping = aes(x = bin_start, y = get(sp))) + geom_col(width = 1, color = "slateblue") +
    scale_x_date(date_labels = "%b %Y")+
    labs(subtitle = name(sp), y = NULL, x = NULL) + 
    annotate("rect", xmin=as.Date('2014-7-20'), xmax=as.Date('2015-02-01'), 
             ymin=-Inf, ymax=Inf, alpha=0.6, fill="dimgray") +
    theme(plot.subtitle = element_text(size = 9, face = "bold"), 
          plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), units = "line"))
  
}

all_presence <- list()
for(s in species) {
  all_presence[[length(all_presence)+1]] <- binnedTimeseries(get(paste0(s,'_binned')),s)
}

final_presence <- wrap_plots(all_presence, nrow = length(all_presence), guides = "collect") 
# -------------- Step 4: VIF for Correlation -------------------
# Including AAO and julian day because have multiple years of data now
# skipping ice difference and ice thickness in models
# not doing ice regime in vif because it's categorical, but will include in GAMs
# including median absolute deviation of EKE to quantify variability in EKE

# LONG-FINNED PILOT WHALE
Gm_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc','temperature_0','salinity_0',
          'EKE_0','chla_0','o2_0','productivity_0','EKE_0_mad', 'bathymetry')
mod_formula <- paste('Gm', "~", paste(Gm_pred, collapse = " + "))
Gm_vif <- glm(as.formula(mod_formula),family=binomial, data = Gm_binned)
vif(Gm_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc  temperature_0 
# 5.728525       1.172656       8.339980       1.457207       4.354805       2.547655      16.491266 
# salinity_0          EKE_0         chla_0           o2_0 productivity_0      EKE_0_mad     bathymetry 
# 16.465565       1.100879       5.096060       7.807557       7.815921       1.868155       7.232413 

# dropping surface temperature
Gm_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc','salinity_0',
             'EKE_0','chla_0','o2_0','productivity_0','EKE_0_mad','bathymetry')
mod_formula <- paste('Gm', "~", paste(Gm_pred, collapse = " + "))
Gm_vif <- glm(as.formula(mod_formula),family=binomial, data = Gm_binned)
vif(Gm_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc     salinity_0 
# 5.456144       1.126727       8.108424       1.448854       4.349912       1.883413      15.779471 
# EKE_0         chla_0           o2_0 productivity_0      EKE_0_mad     bathymetry 
# 1.081411       5.040646       4.731775       4.960660       1.841614       7.083851 

# dropping salinity
Gm_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc',
             'EKE_0','chla_0','o2_0','productivity_0','EKE_0_mad','bathymetry')
mod_formula <- paste('Gm', "~", paste(Gm_pred, collapse = " + "))
Gm_vif <- glm(as.formula(mod_formula),family=binomial, data = Gm_binned)
vif(Gm_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc          EKE_0 
# 3.805776       1.121054       4.948258       1.429686       2.109671       1.813396       1.081518 
# chla_0           o2_0 productivity_0      EKE_0_mad     bathymetry 
# 5.076688       3.430069       4.793530       1.737491       5.271448

# dropping bathymetry
Gm_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc',
             'EKE_0','chla_0','o2_0','productivity_0','EKE_0_mad')
mod_formula <- paste('Gm', "~", paste(Gm_pred, collapse = " + "))
Gm_vif <- glm(as.formula(mod_formula),family=binomial, data = Gm_binned)
vif(Gm_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc          EKE_0 
# 3.735738       1.132143       1.751263       1.364929       2.116814       1.713405       1.080046 
# chla_0           o2_0 productivity_0      EKE_0_mad 
# 4.926068       3.316443       4.946642       1.689770 

# all VIF under five, initial predictors for Gm GAM are julian day, AAO, SSH, FSLE,
#  mixed layer depth, ice concentration, EKE, chlorophyll, oxygen, primary production, EKE variability





# SOUTHERN BOTTLENOSE WHALE
BW29_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc','temperature_0','salinity_0',
             'EKE_0','chla_0','o2_0','productivity_0','EKE_0_mad', 'bathymetry')
mod_formula <- paste('BW29', "~", paste(BW29_pred, collapse = " + "))
BW29_vif <- glm(as.formula(mod_formula),family=binomial, data = BW29_binned)
vif(BW29_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc  temperature_0 
# 7.150689       1.273256       4.971735       1.764209       2.701471       2.250175      10.899265 
# salinity_0          EKE_0         chla_0           o2_0 productivity_0      EKE_0_mad     bathymetry 
# 11.142979       1.114449       5.755051      10.944039       8.699705       1.379283       5.651130 

# dropping surface salinity
BW29_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc','temperature_0',
               'EKE_0','chla_0','o2_0','productivity_0','EKE_0_mad', 'bathymetry')
mod_formula <- paste('BW29', "~", paste(BW29_pred, collapse = " + "))
BW29_vif <- glm(as.formula(mod_formula),family=binomial, data = BW29_binned)
vif(BW29_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc  temperature_0 
# 5.498033       1.248847       3.190564       1.721469       2.055399       2.260342       9.467231 
# EKE_0         chla_0           o2_0 productivity_0      EKE_0_mad     bathymetry 
# 1.098468       5.702439       8.174772       8.036270       1.367164       3.766588 

# dropping surface temperature
BW29_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc',
               'EKE_0','chla_0','o2_0','productivity_0','EKE_0_mad', 'bathymetry')
mod_formula <- paste('BW29', "~", paste(BW29_pred, collapse = " + "))
BW29_vif <- glm(as.formula(mod_formula),family=binomial, data = BW29_binned)
vif(BW29_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc          EKE_0 
# 5.403330       1.205052       2.924650       1.683020       2.085344       2.031606       1.080340 
# chla_0           o2_0 productivity_0      EKE_0_mad     bathymetry 
# 5.384540       5.420739       5.646201       1.375201       3.987482 

# dropping primary production
BW29_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc',
               'EKE_0','chla_0','o2_0','EKE_0_mad', 'bathymetry')
mod_formula <- paste('BW29', "~", paste(BW29_pred, collapse = " + "))
BW29_vif <- glm(as.formula(mod_formula),family=binomial, data = BW29_binned)
vif(BW29_vif)
# julian_day         AAO         SSH        FSLE mixed_layer    ice_conc       EKE_0      chla_0        o2_0 
# 4.479190    1.213812    2.571181    1.593786    2.034043    2.132971    1.076805    1.940039    4.475326 
# EKE_0_mad  bathymetry 
# 1.345338    3.627142  

# all VIF under five, initial predictors for BW29 GAM are julian day, AAO, SSH, FSLE,
#  mixed layer depth, ice concentration, EKE, chlorophyll, oxygen, EKE variability, bathymetry




# GRAY'S AND STRAP-TOOTHED WHALE
BW37_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc','temperature_0','salinity_0',
               'EKE_0','chla_0','o2_0','productivity_0','EKE_0_mad','bathymetry')
mod_formula <- paste('BW37', "~", paste(BW37_pred, collapse = " + "))
BW37_vif <- glm(as.formula(mod_formula),family=binomial, data = BW37_binned)
vif(BW37_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc  temperature_0 
# 17.725603       1.446230       3.862211       1.760413       2.898430       1.557771       7.534608 
# salinity_0          EKE_0         chla_0           o2_0 productivity_0      EKE_0_mad     bathymetry 
# 15.327656       1.296805       7.731137      17.371399      11.502201       1.351688       8.532204 

# removing oxygen concentration
BW37_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc','temperature_0','salinity_0',
               'EKE_0','chla_0','productivity_0','EKE_0_mad','bathymetry')
mod_formula <- paste('BW37', "~", paste(BW37_pred, collapse = " + "))
BW37_vif <- glm(as.formula(mod_formula),family=binomial, data = BW37_binned)
vif(BW37_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc  temperature_0 
# 6.519766       1.425190       3.874463       1.634441       2.397643       1.513046       6.278406 
# salinity_0          EKE_0         chla_0 productivity_0      EKE_0_mad     bathymetry 
# 12.656489       1.287685       7.397869       9.101380       1.314622       7.774989 

# removing salinity
BW37_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc','temperature_0',
               'EKE_0','chla_0','productivity_0','EKE_0_mad','bathymetry')
mod_formula <- paste('BW37', "~", paste(BW37_pred, collapse = " + "))
BW37_vif <- glm(as.formula(mod_formula),family=binomial, data = BW37_binned)
vif(BW37_vif)
# julian_day            AAO            SSH           FSLE    mixed_layer       ice_conc  temperature_0 
# 4.811460       1.425247       3.352221       1.635862       2.172295       1.509936       5.455867 
# EKE_0         chla_0 productivity_0      EKE_0_mad     bathymetry 
# 1.207020       7.159117       9.017689       1.340326       4.114878 

# removing primary production
BW37_pred <- c('julian_day','AAO','SSH','FSLE','mixed_layer','ice_conc','temperature_0',
               'EKE_0','chla_0','EKE_0_mad','bathymetry')
mod_formula <- paste('BW37', "~", paste(BW37_pred, collapse = " + "))
BW37_vif <- glm(as.formula(mod_formula),family=binomial, data = BW37_binned)
vif(BW37_vif)
# julian_day           AAO           SSH          FSLE   mixed_layer      ice_conc temperature_0 
# 3.827289      1.459493      2.615022      1.688835      2.090521      1.620431      3.655930 
# EKE_0        chla_0     EKE_0_mad    bathymetry 
# 1.177380      2.637888      1.341291      3.853325 

# all VIF under five, initial predictors for BW37 GAM are julian day, AAO, SSH, FSLE,
#  mixed layer depth, ice concentration, temperature, EKE, chlorophyll, EKE variability, and bathymetry

# -------------- Step 5: Build GAMs ------------------------
# Function to visualize GAMs on a probability scale with the proper confidence interval
# Run this for each iteration of the model to quickly visualize smooth terms
plotGam <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],seWithMean=TRUE))
}
# Run this if all plots should be in one figure
plotGam1 <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],seWithMean=TRUE,pages=1))
}

# -------------------- Step 5a: Southern Bottlenose Whale GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing method at REML
# Using list of predictors in BW29_pred, adding ice regime as well
# Weighing presence values at 2 (ratio of 0s to 1s is about 2.5)
# running summary() for all individual models,
#   summary(), plotGam1(), AIC(), and concurvity() for all combined models
#   adding on gam.check() for final model

# weighing 1s by 2 to make the ratio of 1s to 0s near 1:1
BW29_binned$weights <- ifelse(BW29_binned$BW29 == 1, 2, 1)

# INDIVIDUAL MODELS

# Julian Day
BW29_gam <- gam(BW29 ~ s(julian_day,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(julian_day, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.4377     0.1174  -3.728 0.000193 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(julian_day) 2.915  2.994  53.95  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.147   Deviance explained = 12.4%
# -REML = 237.96  Scale est. = 1         n = 301

# Antarctic Oscillation Index
BW29_gam <- gam(BW29 ~ s(AAO,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(AAO, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.2567     0.1028  -2.496   0.0125 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(AAO)   1  1.001   0.32   0.572
# 
# R-sq.(adj) =  -0.00255   Deviance explained = 0.0608%
# -REML = 266.28  Scale est. = 1         n = 301

# sea surface height
BW29_gam <- gam(BW29 ~ s(SSH,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(SSH, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.2582     0.1029   -2.51   0.0121 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(SSH) 1.613  1.971  1.025   0.595
# 
# R-sq.(adj) =  -0.00102   Deviance explained = 0.326%
# -REML = 266.27  Scale est. = 1         n = 301

# FSLE
BW29_gam <- gam(BW29 ~ s(FSLE,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(FSLE, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.2760     0.1042  -2.648  0.00809 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(FSLE) 1.734  2.106  7.551  0.0249 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0157   Deviance explained = 1.64%
# -REML = 262.95  Scale est. = 1         n = 301

# mixed layer depth
BW29_gam <- gam(BW29 ~ s(mixed_layer,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.4539     0.1342  -3.382  0.00072 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(mixed_layer) 2.534  2.839  22.78 3.73e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.084   Deviance explained = 7.96%
# -REML = 247.52  Scale est. = 1         n = 301

# sea ice concentration
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -3.512      2.801  -1.254     0.21
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(ice_conc) 1.8  1.976  5.063  0.0607 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.165   Deviance explained = 17.3%
# -REML =  219.9  Scale est. = 1         n = 301
# .............................................
# rerunning ice concentration, manually setting sp to 0.1 to see if it's significant
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1), weights=weights,
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.1166     0.3219  -3.469 0.000522 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc) 1.255  1.453  12.34 0.000689 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.16   Deviance explained = 16.5%
# UBRE = 0.47881  Scale est. = 1         n = 301

# eddy kinetic energy
BW29_gam <- gam(BW29 ~ s(EKE_0,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(EKE_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.2561     0.1028  -2.492   0.0127 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_0)   1  1.001  0.078   0.781
# 
# R-sq.(adj) =  -0.00315   Deviance explained = 0.0149%
# -REML = 266.39  Scale est. = 1         n = 301

# oxygen concentration
BW29_gam <- gam(BW29 ~ s(o2_0,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.2939     0.1055  -2.787  0.00533 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(o2_0) 2.79  2.967  15.51 0.00323 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0305   Deviance explained = 3.21%
# -REML = 261.37  Scale est. = 1         n = 301

# chlorophyll-a concentration
BW29_gam <- gam(BW29 ~ s(chla_0,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(chla_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.3081     0.1078  -2.858  0.00427 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(chla_0) 2.427  2.768  13.32 0.00846 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0306   Deviance explained = 3.11%
# -REML =  260.2  Scale est. = 1         n = 301

# variability in eddy kinetic energy
BW29_gam <- gam(BW29 ~ s(EKE_0_mad,k=4), weights=weights, method='REML',
                family=binomial, data=BW29_binned)
# Formula:
#   BW29 ~ s(EKE_0_mad, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.2714     0.1048  -2.589  0.00962 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(EKE_0_mad) 2.742  2.948  8.803  0.0532 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0195   Deviance explained = 2.33%
# -REML = 263.25  Scale est. = 1         n = 301

# ice regime
BW29_gam <- gam(BW29 ~ ice_regime, weights=weights, method='REML',
                family=binomial, data=BW29_binned) # not a smooth, running as categorical variable
# Formula:
#   BW29 ~ ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           -1.0498     0.3105  -3.381 0.000723 ***
#   ice_regimeincreasing   0.5299     0.4110   1.289 0.197242    
# ice_regimenone         1.1639     0.3424   3.400 0.000675 ***
#   ice_regimestable       0.5592     0.3875   1.443 0.149012    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# R-sq.(adj) =  0.0296   Deviance explained = 3.01%
# -REML = 257.99  Scale est. = 1         n = 301

# bathymetry
BW29_gam <- gam(BW29 ~ bathymetry, weights=weights, method='REML',
                family=binomial, data=BW29_binned) # not including as a smooth because only two depths
# Formula:
#   BW29 ~ bathymetry
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept) -1.8296749  0.8041846  -2.275   0.0229 *
#   bathymetry   0.0019432  0.0009845   1.974   0.0484 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# R-sq.(adj) =  0.00665   Deviance explained = 0.742%
# -REML = 269.13  Scale est. = 1         n = 301




# BUILDING FINAL MODEL
# significant variables: julian day, FSLE, mixed layer depth, ice concentration, oxygen concentration
#   bathymetry, ice regime, chlorophyll-a concentration

# starting with ice concentration and julian day
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(julian_day,k=4), 
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# AIC: 402.9761
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.6301     0.4485  -3.635 0.000278 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.082  1.157  11.91 0.000452 ***
#   s(julian_day) 2.646  2.899  36.71  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.267   Deviance explained = 25.4%
# -REML = 201.55  Scale est. = 1         n = 301
# ...................................................
# concurvity and plots look good

# adding FSLE
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(julian_day,k=4) + s(FSLE,k=4), 
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# AIC: 403.1978
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(FSLE, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.6807     0.4608  -3.647 0.000265 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(ice_conc)   1.073  1.140  11.74 0.00047 ***
#   s(julian_day) 2.661  2.908  39.26 < 2e-16 ***
#   s(FSLE)       1.000  1.000   1.67 0.19635    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.267   Deviance explained = 25.8%
# -REML = 201.83  Scale est. = 1         n = 301
# .................................................
# dropping FSLE because not significant

# dropping FSLE and adding mixed layer depth
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(julian_day,k=4) + s(mixed_layer,k=4), 
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# AIC: 405.0934
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.5254     0.4372  -3.489 0.000485 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(ice_conc)    1.102  1.194  9.577 0.00164 ** 
#   s(julian_day)  2.606  2.875 32.973 2.3e-07 ***
#   s(mixed_layer) 1.786  2.170  2.058 0.41636    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.27   Deviance explained = 25.9%
# -REML = 201.59  Scale est. = 1         n = 301
# .................................................
# removing mixed layer depth because AIC went up and it is not significant

# removing mixed layer, adding oxygen concentration
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(julian_day,k=4) + s(o2_0,k=4), 
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# AIC: 390.7851
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.5442     0.3862  -3.999 6.37e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.183  1.335  20.52 1.22e-05 ***
#   s(julian_day) 2.786  2.954  34.35 3.35e-07 ***
#   s(o2_0)       2.320  2.694  12.13  0.00428 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.296   Deviance explained = 28.9%
# -REML = 195.34  Scale est. = 1         n = 301
# .................................................
# plots look good, keep an eye on concurvity (starting to get potentially problematic)

# adding bathymetry
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(julian_day,k=4) + s(o2_0,k=4) + bathymetry,
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# AIC: 379.6375
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(o2_0, k = 4) + bathymetry
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -5.881744   1.307978  -4.497  6.9e-06 ***
#   bathymetry   0.005125   0.001417   3.616    3e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.107  1.202  17.38 3.48e-05 ***
#   s(julian_day) 2.723  2.928  34.43 1.09e-06 ***
#   s(o2_0)       1.777  2.181  10.34  0.00695 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.317   Deviance explained = 31.2%
# -REML = 193.41  Scale est. = 1         n = 301
# ....................................................
# AIC went down, bathymetry is significant
# however, not keeping it in the model because the effect size of all the smooths shrunk
#    by approx. a factor of 10

# dropping bathymetry, adding chlorophyll-a
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(julian_day,k=4) + s(o2_0,k=4) + 
                  s(chla_0,k=4),
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# AIC: 355.408
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(o2_0, k = 4) + s(chla_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.9089     0.4086  -4.672 2.98e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.174  1.319  22.78 3.95e-06 ***
#   s(julian_day) 2.883  2.987  25.07 1.83e-05 ***
#   s(o2_0)       2.867  2.983  25.18 1.89e-05 ***
#   s(chla_0)     2.822  2.973  19.46 0.000477 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.37   Deviance explained = 36.9%
# -REML = 180.72  Scale est. = 1         n = 301
# ....................................................
# chlorophyll-a has high concurvity with julian day and oxygen concentration
#    not too bad as of right now, but may have to drop chlorophyll later
# AIC went down

# adding ice regime
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(julian_day,k=4) + s(o2_0,k=4) + 
                  s(chla_0,k=4) + ice_regime,
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# AIC: 352.3295
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(o2_0, k = 4) + s(chla_0, k = 4) + ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)           -2.2146     0.7345  -3.015  0.00257 **
#   ice_regimeincreasing   1.0575     0.6691   1.580  0.11402   
# ice_regimenone        -0.1174     0.7767  -0.151  0.87984   
# ice_regimestable      -0.6450     0.5975  -1.079  0.28038   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.095  1.182  20.12 9.03e-06 ***
#   s(julian_day) 2.870  2.985  23.31 3.29e-05 ***
#   s(o2_0)       2.859  2.982  21.88 8.33e-05 ***
#   s(chla_0)     2.818  2.972  18.15 0.000758 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.387   Deviance explained = 38.6%
# -REML = 174.78  Scale est. = 1         n = 301
# ...................................................
# ice regime not significant, not keeping it in model
# going to try removing chlorophyll-a because concurvity continues to be significant

# removing ice regime and chlorophyll-a
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(julian_day,k=4) + s(o2_0,k=4),
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# AIC: 390.7851
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.5442     0.3862  -3.999 6.37e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.183  1.335  20.52 1.22e-05 ***
#   s(julian_day) 2.786  2.954  34.35 3.35e-07 ***
#   s(o2_0)       2.320  2.694  12.13  0.00428 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.296   Deviance explained = 28.9%
# -REML = 195.34  Scale est. = 1         n = 301
# ...................................................
# AIC went up, deviance explained went down
# still, taking away chlorophyll because concurvity is largely resolved + plots have smaller intervals
# will compare to keeping chlorophyll and removing julian day to see if that is better

# keeping chlorophyll, removing julian day
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(chla_0,k=4) + s(o2_0,k=4),
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# AIC: 394.4104
# summary:
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 0.1) + s(chla_0, k = 4) + s(o2_0, 
#                                                              k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.4743     0.3614  -4.079 4.52e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc) 1.199  1.361  19.38 2.22e-05 ***
#   s(chla_0)   2.753  2.950  34.77 7.96e-07 ***
#   s(o2_0)     2.876  2.986  23.25 5.68e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.295   Deviance explained = 28.4%
# -REML = 199.31  Scale est. = 1         n = 301
# ....................................................
# AIC is a little higher, plots are less interpretable
# moving forward with removing chlorophyll, keepign julian day




# REFINING FINAL MODEL
# final variables: sea ice concentration, oxygen concentration, julian day
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.1) + s(julian_day,k=4) + s(o2_0,k=4),
                weights=weights, method='REML', family=binomial, data=BW29_binned)
# gam.check to see if k needs to be adjusted
# Method: REML   Optimizer: outer newton
# full convergence after 5 iterations.
# Gradient range [-9.50591e-09,4.966914e-08]
# (score 195.3367 & scale 1).
# Hessian positive definite, eigenvalue range [0.4256077,0.7819978].
# Model rank =  10 / 10 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                 k'  edf k-index p-value    
# s(ice_conc)   3.00 1.18    0.81  <2e-16 ***
#   s(julian_day) 3.00 2.79    1.06    0.92    
# s(o2_0)       3.00 2.32    0.98    0.50    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# trying higher k for ice concentration because significant in gam.check
BW29_gam <- gam(BW29 ~ s(ice_conc,k=8,sp=0.1) + s(julian_day,k=4) + s(o2_0,k=4),
                weights=weights, method='REML', family=binomial, data=BW29_binned) # AIC: 391.3156
# AIC went up (391), keeping k at 1
# adjusting smoothing parameter for ice concentration (keeping others the same because using REML method)
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=0.01) + s(julian_day,k=4) + s(o2_0,k=4),
                weights=weights, method='REML', family=binomial, data=BW29_binned) # AIC: 391.083
BW29_gam <- gam(BW29 ~ s(ice_conc,k=4,sp=1) + s(julian_day,k=4) + s(o2_0,k=4),
                weights=weights, method='REML', family=binomial, data=BW29_binned) # AIC: 390.471
# AIC decreased and confidence interval narrowed, making sp 1

# final model for BW29 at all sites:
BW29_final <- gam(BW29 ~ s(ice_conc,k=4,sp=1) + s(julian_day,k=4) + s(o2_0,k=4),
                weights=weights, method='REML', family=binomial, data=BW29_binned) # AIC: 390.471
# BW29 final predictors:
BW29_pred <- c('ice_conc','julian_day','o2_0')
# Formula:
#   BW29 ~ s(ice_conc, k = 4, sp = 1) + s(julian_day, k = 4) + s(o2_0, 
#                                                                k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.4866     0.3311   -4.49 7.13e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.025  1.049  28.99  < 2e-16 ***
#   s(julian_day) 2.788  2.954  34.77 9.16e-07 ***
#   s(o2_0)       2.311  2.686  12.17  0.00417 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.296   Deviance explained = 28.9%
# -REML = 195.36  Scale est. = 1         n = 301




# -------------------- Step 5b: Gray's and Strap-Toothed Whale GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing method at REML
# Using list of predictors in BW37_pred, adding ice regime as well
# Weighing presence values at 2 (ratio of 0s to 1s is about 2.5)
# running summary() for all individual models,
#   summary(), plotGam1(), AIC(), and concurvity() for all combined models
#   adding on gam.check() for final model

# weighing 1s by 2 to make the ratio of 1s to 0s near 1:1
BW37_binned$weights <- ifelse(BW37_binned$BW37 == 1, 2, 1)

# INDIVIDUAL MODELS

# julian day
BW37_gam <- gam(BW37 ~ s(julian_day,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(julian_day, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.4619     0.1701  -2.715  0.00662 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(julian_day) 2.907  2.993  28.68 3.34e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.152   Deviance explained = 14.3%
# -REML = 121.03  Scale est. = 1         n = 152

# antarctic oscillation index
BW37_gam <- gam(BW37 ~ s(AAO,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(AAO, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.3233     0.1560  -2.072   0.0383 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(AAO) 2.496  2.828  7.217  0.0924 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0274   Deviance explained = 4.11%
# -REML =  132.4  Scale est. = 1         n = 152

# FSLE
BW37_gam <- gam(BW37 ~ s(FSLE,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(FSLE, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.2686     0.1477  -1.818    0.069 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(FSLE) 2.501  2.819  7.417  0.0337 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0447   Deviance explained = 4.38%
# -REML = 132.32  Scale est. = 1         n = 152

# mixed layer
BW37_gam <- gam(BW37 ~ s(mixed_layer,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.3943     0.1627  -2.424   0.0154 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(mixed_layer)   1      1   18.2 1.99e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.101   Deviance explained = 8.83%
# -REML = 123.63  Scale est. = 1         n = 152

# sea ice concentration
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -26.05      25.88  -1.007    0.314
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc) 1.762  1.944  2.239   0.257
# 
# R-sq.(adj) =  0.185   Deviance explained = 19.6%
# -REML = 108.21  Scale est. = 1         n = 152
# ............................................
# manually setting sp to fix large confidence interval
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1), 
                family=binomial,weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.4173     0.5033  -2.816  0.00487 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc) 1.023  1.045  9.239 0.00169 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.172   Deviance explained = 17.7%
# UBRE = 0.47611  Scale est. = 1         n = 152

# sea surface temperature
BW37_gam <- gam(BW37 ~ s(temperature_0,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(temperature_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.5655     0.1849  -3.058  0.00223 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(temperature_0) 2.798  2.969   37.6  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.222   Deviance explained = 19.9%
# -REML = 112.36  Scale est. = 1         n = 152

# eddy kinetic energy
BW37_gam <- gam(BW37 ~ s(EKE_0,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(EKE_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.2540     0.1461  -1.738   0.0822 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(EKE_0) 1.001  1.002  4.455  0.0348 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0169   Deviance explained = 1.77%
# -REML = 133.45  Scale est. = 1         n = 152

# chlorophyll-a concentration
BW37_gam <- gam(BW37 ~ s(chla_0,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(chla_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.2579     0.1462  -1.764   0.0778 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(chla_0) 2.071  2.446  4.706   0.173
# 
# R-sq.(adj) =  0.0129   Deviance explained = 2.02%
# -REML = 134.45  Scale est. = 1         n = 152

# EKE variability
BW37_gam <- gam(BW37 ~ s(EKE_0_mad,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ s(EKE_0_mad, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.3639     0.1584  -2.296   0.0217 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(EKE_0_mad) 2.657  2.925  26.89 1.86e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.126   Deviance explained = 11.5%
# -REML = 122.89  Scale est. = 1         n = 152

# bathymetry
BW37_gam <- gam(BW37 ~ bathymetry, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ bathymetry
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  7.438383   1.136190   6.547 5.88e-11 ***
#   bathymetry  -0.008838   0.001322  -6.686 2.30e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# R-sq.(adj) =  0.255   Deviance explained = 20.6%
# -REML = 112.77  Scale est. = 1         n = 152

# ice regime
BW37_gam <- gam(BW37 ~ ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# Formula:
#   BW37 ~ ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)           -1.6094     0.5477  -2.938  0.00330 **
#   ice_regimeincreasing   1.2040     0.6625   1.817  0.06916 . 
# ice_regimenone         1.7838     0.5864   3.042  0.00235 **
#   ice_regimestable       1.2379     0.6201   1.996  0.04590 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# R-sq.(adj) =  0.0392   Deviance explained = 4.67%
# -REML = 128.27  Scale est. = 1         n = 152



# BUILDING FINAL MODEL
# significant variables: julian day, FSLE, mixed layer depth, ice concentration, temperature
#   EKE, EKE variability, bathymetry, ice regime

# starting with julian day and ice concentration
BW37_gam <- gam(BW37 ~ s(julian_day,k=4) + s(ice_conc,k=4,sp=0.1), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 214.4342
# summary:
# Formula:
#   BW37 ~ s(julian_day, k = 4) + s(ice_conc, k = 4, sp = 0.1)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.3590     0.4601  -2.954  0.00314 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(julian_day) 2.819  2.973 11.720 0.00545 **
#   s(ice_conc)   1.036  1.071  7.261 0.00596 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.242   Deviance explained = 23.6%
# -REML = 107.09  Scale est. = 1         n = 152
# ................................................
# concurvity and plot looks good

# adding FSLE
BW37_gam <- gam(BW37 ~ s(julian_day,k=4) + s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 208.8272
# summary:
# Formula:
#   BW37 ~ s(julian_day, k = 4) + s(ice_conc, k = 4, sp = 0.1) + 
#   s(FSLE, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.4785     0.4826  -3.064  0.00218 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(julian_day) 2.831  2.976 12.833 0.00318 **
#   s(ice_conc)   1.032  1.063  7.605 0.00486 **
#   s(FSLE)       1.000  1.000  6.845 0.00889 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.275   Deviance explained = 26.5%
# -REML = 104.14  Scale est. = 1         n = 152
# ......................................................
# FSLE trend looks weak in plot, consider dropping later on
# concurvity looks good

# adding mixed layer depth
BW37_gam <- gam(BW37 ~ s(julian_day,k=4) + s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4) +
                  s(mixed_layer,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 210.5139
# summary:
# Formula:
#   BW37 ~ s(julian_day, k = 4) + s(ice_conc, k = 4, sp = 0.1) + 
#   s(FSLE, k = 4) + s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.4851     0.4875  -3.046  0.00232 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(julian_day)  2.798  2.963 10.400 0.00960 **
#   s(ice_conc)    1.030  1.058  7.836 0.00429 **
#   s(FSLE)        1.000  1.000  5.782 0.01620 * 
#   s(mixed_layer) 1.699  2.048  1.646 0.46923   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.276   Deviance explained = 27.4%
# -REML =  103.9  Scale est. = 1         n = 152
# ..................................................
# removing mixed layer, AIC went up and not significant

# removing mixed layer depth, adding temperature
BW37_gam <- gam(BW37 ~ s(julian_day,k=4) + s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4) +
                  s(temperature_0,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 204.4069
# summary:
# Formula:
#   BW37 ~ s(julian_day, k = 4) + s(ice_conc, k = 4, sp = 0.1) + 
#   s(FSLE, k = 4) + s(temperature_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.2435     0.4655  -2.671  0.00756 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(julian_day)    2.556  2.851  6.493  0.0460 *
#   s(ice_conc)      1.043  1.085  3.367  0.0688 .
# s(FSLE)          1.001  1.001  4.672  0.0307 *
#   s(temperature_0) 2.581  2.874  6.764  0.0842 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.306   Deviance explained = 30.2%
# -REML = 100.78  Scale est. = 1         n = 152
# ......................................................
# temperature has high convurvity with temperature
# will try dropping julian day and see how it compares to model without temperature

# dropping julian day
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4) +
                  s(temperature_0,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 210.0211
# summary:
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1) + s(FSLE, k = 4) + s(temperature_0, 
#                                                            k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -1.049      0.391  -2.682  0.00731 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc)      1.077  1.148  2.921 0.10247   
# s(FSLE)          1.000  1.000  8.522 0.00351 **
#   s(temperature_0) 2.734  2.947 13.596 0.00389 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.262   Deviance explained = 26.1%
# -REML = 103.95  Scale est. = 1         n = 152
# ...................................................
# keeping julian day and dropping temperature
#   higher AIC if we drop julian day, ice concentration loses its significance

# dropping temperature, keeping julian day, adding EKE
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4) +
                  s(julian_day,k=4) + s(EKE_0,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 209.6743
# summary:
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1) + s(FSLE, k = 4) + s(julian_day, 
#                                                            k = 4) + s(EKE_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.6045     0.5354  -2.997  0.00273 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc)   1.022  1.043  7.389 0.00546 **
#   s(FSLE)       1.000  1.000  5.393 0.02023 * 
#   s(julian_day) 2.802  2.967 11.078 0.00706 **
#   s(EKE_0)      2.107  2.491  3.813 0.34638   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.285   Deviance explained =   28%
# -REML = 104.16  Scale est. = 1         n = 152
# .......................................................
# dropping EKE, weak trend and insignificant

# dropping EKE strength, adding EKE variability
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4) +
                  s(julian_day,k=4) + s(EKE_0_mad,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 201.7507
# summary:
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1) + s(FSLE, k = 4) + s(julian_day, 
#                                                            k = 4) + s(EKE_0_mad, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.5326     0.5416   -2.83  0.00466 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(ice_conc)   1.023  1.045  5.872 0.01332 * 
#   s(FSLE)       1.000  1.000  3.264 0.07082 . 
# s(julian_day) 2.711  2.932  5.995 0.07359 . 
# s(EKE_0_mad)  2.452  2.811 12.235 0.00789 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.331   Deviance explained = 31.2%
# -REML = 99.651  Scale est. = 1         n = 152
# ....................................................
# FSLE and julian day dropped in significance
# removing FSLE since it has a weak trend

# removing FSLE
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + 
                  s(julian_day,k=4) + s(EKE_0_mad,k=4), 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 205.8705
# summary:
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(EKE_0_mad, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.6451     0.5715  -2.879  0.00399 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.016  1.031  8.062 0.003492 ** 
#   s(julian_day) 1.000  1.001  0.569 0.450876    
# s(EKE_0_mad)  2.586  2.892 20.841 0.000235 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.285   Deviance explained = 27.5%
# -REML = 100.95  Scale est. = 1         n = 152
# .............................................
# AIC went up slightly, plot shows julian day trend disappearing
# since there is large confidence intervals for EKE variability, removing that variable
#    this reverts back to model at line 1187
# will also add FSLE back in (model from line 1213) since it had better AIC


# removing EKE variability, adding FSLE back, adding bathymetry
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + 
                  s(julian_day,k=4) + s(FSLE,k=4) + bathymetry, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 186.3372
# summary:
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(FSLE, k = 4) + bathymetry
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  7.343483   1.768585   4.152 3.29e-05 ***
#   bathymetry  -0.009900   0.002009  -4.928 8.31e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(ice_conc)   1.030  1.058  4.715  0.0269 *
#   s(julian_day) 2.198  2.553  1.299  0.4536  
# s(FSLE)       1.162  1.303  0.948  0.3672  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.366   Deviance explained = 35.5%
# -REML =  95.16  Scale est. = 1         n = 152
# .................................................
# plots look crazy (flat lines, weird confidence intervals)
# removing bathymetry (there are only two bathymetry levels anyways)

# removing bathymetry, adding ice regime
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + 
                  s(julian_day,k=4) + s(FSLE,k=4) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 193.9988
# summary:
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   s(FSLE, k = 4) + ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)          -1.99037    1.06124  -1.876   0.0607 .
# ice_regimeincreasing  1.95667    0.92015   2.126   0.0335 *
#   ice_regimenone       -2.10162    1.02451  -2.051   0.0402 *
#   ice_regimestable     -0.09972    0.83230  -0.120   0.9046  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.006  1.013  7.699 0.005024 ** 
#   s(julian_day) 2.892  2.990 18.890 0.000245 ***
#   s(FSLE)       1.000  1.000 10.090 0.001492 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.343   Deviance explained = 34.2%
# -REML = 92.066  Scale est. = 1         n = 152
# ....................................................
# FSLE trend weak -> will try removing
# keeping ice regime because trends are significant, plots look good, AIC went down

# removing FSLE
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + 
                  s(julian_day,k=4) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC: 203.715
# summary:
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1) + s(julian_day, k = 4) + 
#   ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)           -1.7103     0.9850  -1.736   0.0825 .
# ice_regimeincreasing   1.7692     0.9176   1.928   0.0538 .
# ice_regimenone        -1.8222     0.9996  -1.823   0.0683 .
# ice_regimestable      -0.2554     0.8297  -0.308   0.7583  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.009  1.019  7.423 0.005831 ** 
#   s(julian_day) 2.880  2.988 16.726 0.000636 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.295   Deviance explained = 29.9%
# -REML =  97.14  Scale est. = 1         n = 152
# ......................................................
# removing FSLE made AIC go up and got rid of significance of ice regime
# keeping FSLE in the model



# REFINING FINAL MODEL
# final variables: sea ice concentration, julian day, ice regime, FSLE
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4) + # AIC: 193.9988
                  s(julian_day,k=4) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# first, running gam.check():
# Method: REML   Optimizer: outer newton
# full convergence after 9 iterations.
# Gradient range [-2.874799e-05,7.225422e-07]
# (score 92.06591 & scale 1).
# Hessian positive definite, eigenvalue range [2.874665e-05,0.7954052].
# Model rank =  13 / 13 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                 k'  edf k-index p-value
# s(ice_conc)   3.00 1.01    0.97    0.33
# s(FSLE)       3.00 1.00    1.05    0.75
# s(julian_day) 3.00 2.89    1.01    0.57
# ..................................................
# gam.check output looks good, not messing with k values

# trying different smoothing for ice concentration
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.01) + s(FSLE,k=4) + # AIC: 194.1404
                  s(julian_day,k=4) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=1) + s(FSLE,k=4) + # AIC: 193.9806
                  s(julian_day,k=4) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC stays same or increases, keeping current sp for ice concentration

# trying different smoothing for FSLE (confidence interval is large)
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4,sp=0.1) + # AIC: 195.9873
                  s(julian_day,k=4) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4,sp=1) + # AIC: 194.6002
                  s(julian_day,k=4) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC increases, keeping REML sp for FSLE

# trying different smoothing for julian day (confidence interval is large at certain points)
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4) + # AIC: 208.2606
                  s(julian_day,k=4,sp=0.1) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
BW37_gam <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4) + # AIC: 213.0995
                  s(julian_day,k=4,sp=1) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# AIC increases, keeping REML sp for julian day

# final model for BW37 click type
BW37_final <- gam(BW37 ~ s(ice_conc,k=4,sp=0.1) + s(FSLE,k=4) + # AIC: 193.9988
                  s(julian_day,k=4) + ice_regime, 
                family=binomial, method='REML', weights=weights, data=BW37_binned)
# final predictors:
BW37_pred <- c('ice_conc','FSLE','julian_day') 
# not adding ice regime to list because can't be visualized as a smooth term
# summary:
# Formula:
#   BW37 ~ s(ice_conc, k = 4, sp = 0.1) + s(FSLE, k = 4) + s(julian_day, 
#                                                            k = 4) + ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)          -1.99037    1.06124  -1.876   0.0607 .
# ice_regimeincreasing  1.95667    0.92015   2.126   0.0335 *
#   ice_regimenone       -2.10162    1.02451  -2.051   0.0402 *
#   ice_regimestable     -0.09972    0.83230  -0.120   0.9046  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)   1.006  1.013  7.699 0.005024 ** 
#   s(FSLE)       1.000  1.000 10.090 0.001492 ** 
#   s(julian_day) 2.892  2.990 18.890 0.000245 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.343   Deviance explained = 34.2%
# -REML = 92.066  Scale est. = 1         n = 152


# -------------------- Step 5c: Long-Finned Pilot Whale GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Setting initial knots at 4, smoothing method at REML
# Using list of predictors in Gm_pred, adding ice regime as well
# Weighing presence values at 5 (ratio of 0s to 1s is about 5.3)
# running summary() for all individual models,
#   summary(), plotGam1(), AIC(), and concurvity() for all combined models
#   adding on gam.check() for final model

# weighing 1s by 5 to make the ratio of 1s to 0s near 1:1
Gm_binned$weights <- ifelse(Gm_binned$Gm == 1, 5, 1)

# INDIVIDUAL MODELS
# julian day
Gm_gam <- gam(Gm ~ s(julian_day,k=4), 
                family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(julian_day, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.5943     0.1266  -4.693 2.69e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(julian_day) 2.856  2.984     91  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.204   Deviance explained = 17.8%
# -REML = 347.16  Scale est. = 1         n = 366

# antarctic oscillation index
Gm_gam <- gam(Gm ~ s(AAO,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(AAO, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.1152     0.0851  -1.353    0.176
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(AAO) 2.836   2.98  15.05 0.00341 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.017   Deviance explained = 2.32%
# -REML = 411.33  Scale est. = 1         n = 366

# sea surface height
Gm_gam <- gam(Gm ~ s(SSH,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(SSH, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept) -0.16161    0.08623  -1.874   0.0609 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(SSH) 2.488  2.815  42.04  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0617   Deviance explained = 5.78%
# -REML = 395.77  Scale est. = 1         n = 366

# FSLE
Gm_gam <- gam(Gm ~ s(FSLE,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(FSLE, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.06297    0.08195  -0.768    0.442
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(FSLE)   1  1.001  1.207   0.272
# 
# R-sq.(adj) =  -0.00101   Deviance explained = 0.146%
# -REML = 416.79  Scale est. = 1         n = 366


# mixed layer depth
Gm_gam <- gam(Gm ~ s(mixed_layer,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept) -0.15582    0.08927  -1.746   0.0809 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(mixed_layer) 2.776  2.962  19.12 0.000417 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0274   Deviance explained = 3.27%
# -REML = 407.04  Scale est. = 1         n = 366

# sea ice concentration
Gm_gam <- gam(Gm ~ s(ice_conc,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.5586     0.1744  -3.202  0.00136 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc) 2.496  2.803  31.78 2.33e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.107   Deviance explained =   11%
# -REML = 373.58  Scale est. = 1         n = 366

# EKE
Gm_gam <- gam(Gm ~ s(EKE_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(EKE_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.08698    0.08308  -1.047    0.295
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(EKE_0) 2.393  2.759  11.48  0.0144 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0111   Deviance explained = 1.48%
# -REML = 413.15  Scale est. = 1         n = 366

# chlorophyll
Gm_gam <- gam(Gm ~ s(chla_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(chla_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.4810     0.1105  -4.353 1.34e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(chla_0) 1.002  1.004  72.86  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.146   Deviance explained = 13.5%
# -REML = 360.88  Scale est. = 1         n = 366

# oxygen
Gm_gam <- gam(Gm ~ s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.9031     0.1764  -5.119 3.07e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(o2_0) 2.884  2.989  85.61  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.219   Deviance explained =   20%
# -REML = 337.95  Scale est. = 1         n = 366

# net primary production
Gm_gam <- gam(Gm ~ s(productivity_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(productivity_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.4635     0.1158  -4.004 6.23e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(productivity_0) 2.355  2.688  86.64  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.158   Deviance explained = 14.7%
# -REML = 358.01  Scale est. = 1         n = 366

# variation in EKE
Gm_gam <- gam(Gm ~ s(EKE_0_mad,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ s(EKE_0_mad, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.13536    0.08499  -1.593    0.111
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(EKE_0_mad) 2.602  2.885  30.96 7.09e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0457   Deviance explained = 4.36%
# -REML = 401.96  Scale est. = 1         n = 366

# ice regime
Gm_gam <- gam(Gm ~ ice_regime, 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# Formula:
#   Gm ~ ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           -0.1636     0.2027  -0.807 0.419539    
# ice_regimeincreasing  -0.2700     0.2970  -0.909 0.363212    
# ice_regimenone        -0.2523     0.2420  -1.043 0.297125    
# ice_regimestable       0.8655     0.2595   3.336 0.000851 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# R-sq.(adj) =  0.041   Deviance explained = 4.14%
# -REML = 400.39  Scale est. = 1         n = 366

# BUILDING FINAL MODEL
# significant variables: julian day, AAO, sea surface height, mixed layer, sea ice concentration 
#   EKE, chlorophyll concentration, oxygen concentration, net primary production, EKE variability,
#   ice regime
# adding one at a time to build final model

# starting with ice concentration and julian day
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 576.2755
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.1218     0.1684  -6.661 2.73e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(julian_day) 2.921  2.995 118.07  <2e-16 ***
#   s(ice_conc)   1.001  1.001  59.64  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.355   Deviance explained = 31.6%
# -REML = 290.77  Scale est. = 1         n = 366
# ..................................................
# plots and concurvity look great

# adding AAO
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4) + s(AAO,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 571.1718
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4) + s(AAO, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -1.188      0.176  -6.754 1.43e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq p-value    
# s(julian_day) 2.924  2.995 114.680  <2e-16 ***
#   s(ice_conc)   1.001  1.002  59.499  <2e-16 ***
#   s(AAO)        2.761  2.959   8.466  0.0382 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.368   Deviance explained =   33%
# -REML = 289.81  Scale est. = 1         n = 366
# ....................................................
# plots look great, concurvity as well

# adding sea surface height
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4) + s(AAO,k=4) +
                s(SSH,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 559.6934
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4) + s(AAO, k = 4) + 
#   s(SSH, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.1425     0.1692  -6.753 1.45e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value    
# s(julian_day) 2.937  2.997 115.460  < 2e-16 ***
#   s(ice_conc)   1.000  1.001  59.558  < 2e-16 ***
#   s(AAO)        2.703  2.939   7.442 0.052713 .  
# s(SSH)        1.000  1.000  13.291 0.000267 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.376   Deviance explained = 34.6%
# -REML = 284.22  Scale est. = 1         n = 366
# ......................................................
# plots and concurvity looks great
# even though AAO is currently not significant, leaving it be for now
#    because it is barely over p=0.05

# adding mixed layer depth
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4) + s(AAO,k=4) +
                s(SSH,k=4) + s(mixed_layer,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 552.1255
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4) + s(AAO, k = 4) + 
#   s(SSH, k = 4) + s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.1535     0.1665   -6.93 4.22e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(julian_day)  2.936  2.995 91.906  < 2e-16 ***
#   s(ice_conc)    1.000  1.000 53.201  < 2e-16 ***
#   s(AAO)         2.600  2.891  5.908   0.1165    
# s(SSH)         1.000  1.001 16.045 6.26e-05 ***
#   s(mixed_layer) 2.519  2.843  8.941   0.0202 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.387   Deviance explained = 36.2%
# -REML = 280.62  Scale est. = 1         n = 366
# .....................................................
# removing AAO because no longer significant
# plots and concurvity look good

# removing AAO
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4) +
                s(SSH,k=4) + s(mixed_layer,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 551.9513
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4) + s(SSH, k = 4) + 
#   s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.1313     0.1625   -6.96  3.4e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(julian_day)  2.923  2.992  93.84  < 2e-16 ***
#   s(ice_conc)    1.000  1.001  51.49  < 2e-16 ***
#   s(SSH)         2.567  2.862  22.92 0.000605 ***
#   s(mixed_layer) 2.588  2.879  10.13 0.011292 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.381   Deviance explained =   36%
# -REML = 280.47  Scale est. = 1         n = 366
# .......................................................
# AIC went down, ready to add next variable

# adding EKE
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4) +
                s(SSH,k=4) + s(mixed_layer,k=4) + s(EKE_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 542.9903
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4) + s(SSH, k = 4) + 
#   s(mixed_layer, k = 4) + s(EKE_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.1834     0.1672  -7.079 1.46e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(julian_day)  2.929  2.993  96.62 < 2e-16 ***
#   s(ice_conc)    1.000  1.001  51.92 < 2e-16 ***
#   s(SSH)         2.519  2.832  21.58 0.00118 ** 
#   s(mixed_layer) 2.602  2.885  12.01 0.00481 ** 
#   s(EKE_0)       2.197  2.502  12.96 0.00461 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.397   Deviance explained = 37.6%
# -REML = 276.45  Scale est. = 1         n = 366
# .................................................................
# EKE trend weak (high error at end), bu keeping it for now
# concurvity good

# adding chlorophyll
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4) + s(chla_0,k=4) +
                s(SSH,k=4) + s(mixed_layer,k=4) + s(EKE_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 521.6007
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4) + s(chla_0, k = 4) + 
#   s(SSH, k = 4) + s(mixed_layer, k = 4) + s(EKE_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.2435     0.1733  -7.175 7.21e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(julian_day)  2.814  2.972 17.480 0.000585 ***
#   s(ice_conc)    1.000  1.000 65.078  < 2e-16 ***
#   s(chla_0)      1.000  1.000 19.916 8.55e-06 ***
#   s(SSH)         1.001  1.002  6.849 0.008877 ** 
#   s(mixed_layer) 2.144  2.548 10.061 0.007502 ** 
#   s(EKE_0)       2.122  2.429  8.686 0.026536 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.415   Deviance explained = 39.9%
# -REML = 262.42  Scale est. = 1         n = 366
# ...................................................
# potentially concerning chlorophyll and julian day concurvity, not an issue yet
# plots look good

# adding oxygen
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4) + s(chla_0,k=4) +
                s(SSH,k=4) + s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 491.9063
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4) + s(chla_0, k = 4) + 
#   s(SSH, k = 4) + s(mixed_layer, k = 4) + s(EKE_0, k = 4) + 
#   s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.5362     0.2081  -7.381 1.57e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(julian_day)  2.742  2.947 17.216 0.000661 ***
#   s(ice_conc)    1.000  1.000 52.279  < 2e-16 ***
#   s(chla_0)      2.290  2.680 15.112 0.002275 ** 
#   s(SSH)         1.001  1.002  4.717 0.029861 *  
#   s(mixed_layer) 2.641  2.900 11.819 0.004196 ** 
#   s(EKE_0)       2.123  2.428  8.248 0.026794 *  
#   s(o2_0)        2.612  2.870 30.694 7.23e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.463   Deviance explained = 44.7%
# -REML = 248.69  Scale est. = 1         n = 366
# ...............................................
# oxygen & chlorophyll and julian day & chlorophyll have concurvity
# removing chlorophyll, its plot also has large confidence interval

# removing chlorophyll
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4) + 
                s(SSH,k=4) + s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 499.2806
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4) + s(SSH, k = 4) + 
#   s(mixed_layer, k = 4) + s(EKE_0, k = 4) + s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.5712     0.2392  -6.568 5.09e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(julian_day)  2.912  2.992  46.39 < 2e-16 ***
#   s(ice_conc)    2.519  2.822  49.42 < 2e-16 ***
#   s(SSH)         1.000  1.000  10.61 0.00112 ** 
#   s(mixed_layer) 2.734  2.944  11.40 0.00685 ** 
#   s(EKE_0)       2.243  2.553  14.16 0.00253 ** 
#   s(o2_0)        2.726  2.943  32.91 < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.456   Deviance explained = 43.7%
# -REML = 255.21  Scale est. = 1         n = 366
# ....................................................
# AIC went slightly up, but okay
# plots look good

# adding net primary production
Gm_gam <- gam(Gm ~ s(julian_day,k=4) + s(ice_conc,k=4) + s(productivity_0,k=4) +
                s(SSH,k=4) + s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 490.5631
# summary:
# Formula:
#   Gm ~ s(julian_day, k = 4) + s(ice_conc, k = 4) + s(productivity_0, 
#                                                      k = 4) + s(SSH, k = 4) + s(mixed_layer, k = 4) + s(EKE_0, 
#                                                                                                         k = 4) + s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7436     0.2632  -6.625 3.47e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(julian_day)     2.673  2.935 19.939 0.000123 ***
#   s(ice_conc)       1.000  1.001 52.438  < 2e-16 ***
#   s(productivity_0) 2.694  2.922 14.544 0.002541 ** 
#   s(SSH)            1.000  1.000  5.713 0.016847 *  
#   s(mixed_layer)    2.674  2.918  9.806 0.011363 *  
#   s(EKE_0)          2.101  2.400  7.823 0.031000 *  
#   s(o2_0)           2.616  2.887 33.203 2.09e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.458   Deviance explained = 44.9%
# -REML = 248.53  Scale est. = 1         n = 366
# .......................................................
# plots look great
# problematic julian day and productivity concurvity
# going to try removing julian day

# removing julian day
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) +
                s(SSH,k=4) + s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 503.1032
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(productivity_0, k = 4) + s(SSH, k = 4) + 
#   s(mixed_layer, k = 4) + s(EKE_0, k = 4) + s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7316     0.2728  -6.348 2.18e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       2.641  2.899 47.231  < 2e-16 ***
#   s(productivity_0) 2.723  2.936 44.122  < 2e-16 ***
#   s(SSH)            2.501  2.825  8.185 0.098752 .  
# s(mixed_layer)    2.736  2.945  8.068 0.055175 .  
# s(EKE_0)          2.142  2.436  9.973 0.013314 *  
#   s(o2_0)           2.573  2.876 22.034 0.000927 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.451   Deviance explained = 43.6%
# -REML = 256.24  Scale est. = 1         n = 366
# .................................................
# concurvity all looks great
# plots are good
# removing SSH since no longer significant and plot trend isn't great

# removing SSH
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 505.6918
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(productivity_0, k = 4) + s(mixed_layer, 
#                                                          k = 4) + s(EKE_0, k = 4) + s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7096     0.2727  -6.269 3.64e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       2.607  2.879 52.358  < 2e-16 ***
#   s(productivity_0) 2.698  2.924 40.740  < 2e-16 ***
#   s(mixed_layer)    2.779  2.961  8.827 0.024929 *  
#   s(EKE_0)          2.111  2.403  8.939 0.020043 *  
#   s(o2_0)           2.607  2.891 25.348 0.000209 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.44   Deviance explained = 42.6%
# -REML = 256.91  Scale est. = 1         n = 366
# .....................................................
# AIC went up slightly but plot trends are generally clear (Except for EKE)
# keeping SSH in the model for now and adding more variables
#   will evaluate later if it needs to be removed

# adding variability in EKE
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 471.6393
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(productivity_0, k = 4) + s(EKE_0_mad, 
#                                                          k = 4) + s(mixed_layer, k = 4) + s(EKE_0, k = 4) + s(o2_0, 
#                                                                                                               k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.8790     0.2694  -6.974 3.08e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       1.055  1.107 46.978  < 2e-16 ***
#   s(productivity_0) 2.696  2.925 38.395  < 2e-16 ***
#   s(EKE_0_mad)      2.719  2.940 34.705 1.34e-06 ***
#   s(mixed_layer)    2.856  2.982 14.182 0.002091 ** 
#   s(EKE_0)          2.073  2.366  7.023 0.046831 *  
#   s(o2_0)           2.578  2.869 23.989 0.000388 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.472   Deviance explained =   47%
# -REML = 240.48  Scale est. = 1         n = 366
# ....................................................
# AIC significantly went down
# concurvity looks good, plots are good

# adding ice regime
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4) + ice_regime, 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 470.7947
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(productivity_0, k = 4) + s(EKE_0_mad, 
#                                                          k = 4) + s(mixed_layer, k = 4) + s(EKE_0, k = 4) + s(o2_0, 
#                                                                                                               k = 4) + ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           -1.8750     0.4813  -3.895  9.8e-05 ***
#   ice_regimeincreasing  -0.7416     0.4642  -1.597    0.110    
# ice_regimenone         0.1080     0.6437   0.168    0.867    
# ice_regimestable       0.4054     0.4369   0.928    0.353    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       1.000  1.000  36.18  < 2e-16 ***
#   s(productivity_0) 2.692  2.923  26.29 9.69e-06 ***
#   s(EKE_0_mad)      2.728  2.944  34.19 8.23e-07 ***
#   s(mixed_layer)    2.860  2.982  14.09  0.00215 ** 
#   s(EKE_0)          2.106  2.412   7.33  0.04154 *  
#   s(o2_0)           2.576  2.865  24.73  0.00036 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.479   Deviance explained = 47.8%
# -REML = 236.85  Scale est. = 1         n = 366
# .................................................
# Ice regime not significant, not keeping it in the model

# try dropping EKE because it has large error, almost not significant
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 474.8803
# summary:
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(productivity_0, k = 4) + s(EKE_0_mad, 
#                                                          k = 4) + s(mixed_layer, k = 4) + s(o2_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.8143     0.2652  -6.841 7.89e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       1.000  1.000  50.19  < 2e-16 ***
#   s(productivity_0) 2.732  2.941  40.07  < 2e-16 ***
#   s(EKE_0_mad)      2.724  2.942  35.17 1.28e-06 ***
#   s(mixed_layer)    2.844  2.979  12.90 0.003757 ** 
#   s(o2_0)           2.567  2.864  23.43 0.000489 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.463   Deviance explained =   46%
# -REML =    242  Scale est. = 1         n = 366

# REFINING FINAL MODEL
# final significant variables: sea ice concentration, primary production,
#   variability in EKE, mixed layer depth, EKE, oxygen concentration
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# AIC: 470.7947
# running gam.check()
# gam.check outputs look good, not messing with k values

# trying different smoothing for each variable

# sea ice concentration
Gm_gam <- gam(Gm ~ s(ice_conc,k=4,sp=0.1) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned) # AIC: 471.3038
Gm_gam <- gam(Gm ~ s(ice_conc,k=4,sp=1) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned) # AIC: 471.7232
# not changing sea ice smoothing because AIC increased

# primary production
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4,sp=0.1) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), # AIC: 476.5538
              family=binomial, method='REML', weights=weights, data=Gm_binned)
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4,sp=1) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), # AIC: 478.5523
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# not changing smoothing because AIC increased

# variability in EKE
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4,sp=0.1) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), # AIC: 475.2021
              family=binomial, method='REML', weights=weights, data=Gm_binned)
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4,sp=1) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), # AIC: 493.6372
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# not changing because AIC increased

# mixed layer depth
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4,sp=0.1) + s(EKE_0,k=4) + s(o2_0,k=4), # AIC: 481.0849
              family=binomial, method='REML', weights=weights, data=Gm_binned)
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4,sp=1) + s(EKE_0,k=4) + s(o2_0,k=4), # AIC: 486.955
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# not changing because AIC increased

# eddy kinetic energy
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4,sp=0.1) + s(o2_0,k=4), # AIC: 472.4249
              family=binomial, method='REML', weights=weights, data=Gm_binned)
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4,sp=1) + s(o2_0,k=4), # AIC: 475.9879
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# not changing because AIC increased

# oxygen concentration
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4,sp=0.1), # AIC: 475.2331
              family=binomial, method='REML', weights=weights, data=Gm_binned)
Gm_gam <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4,sp=1), # AIC: 486.2469
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# not changing because AIC increased


# final model: k=4, smoothing via REML for all variables
Gm_final <- gam(Gm ~ s(ice_conc,k=4) + s(productivity_0,k=4) + s(EKE_0_mad,k=4) +
                s(mixed_layer,k=4) + s(EKE_0,k=4) + s(o2_0,k=4), 
              family=binomial, method='REML', weights=weights, data=Gm_binned)
# final predictors:
Gm_pred <- c('ice_conc','productivity_0','mixed_layer','EKE_0_mad','o2_0')
# Formula:
#   Gm ~ s(ice_conc, k = 4) + s(productivity_0, k = 4) + s(EKE_0_mad, 
#                                                          k = 4) + s(mixed_layer, k = 4) + s(EKE_0, k = 4) + s(o2_0, 
#                                                                                                               k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.8790     0.2694  -6.974 3.08e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(ice_conc)       1.055  1.107 46.978  < 2e-16 ***
#   s(productivity_0) 2.696  2.925 38.395  < 2e-16 ***
#   s(EKE_0_mad)      2.719  2.940 34.705 1.34e-06 ***
#   s(mixed_layer)    2.856  2.982 14.182 0.002091 ** 
#   s(EKE_0)          2.073  2.366  7.023 0.046831 *  
#   s(o2_0)           2.578  2.869 23.989 0.000388 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.472   Deviance explained =   47%
# -REML = 240.48  Scale est. = 1         n = 366
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
    plot_annotation(title = paste0(sp," Presence (",deviance,"% Deviance Explained)"))
  
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
  } else if(paste(var) == 'bathymetry') {
    return('Bathymetry (m)')
  } else if(paste(var) == 'ice_regime') {
    return('Sea Ice Regime (increasing/decreasing/stable/not present)')
  } else if(paste(var) == 'EKE_0_mad') {
    return('Variability in EKE (cm^2/s^2)')
  } else if(paste(var) == 'AAO') {
    return("Antarctic Oscillation Index")
  }
}

# Generating visualizations for each site's final model
Gm_plots <- visualizeGAM(Gm_final, Gm_pred, 'Long-finned Pilot Whale')
BW37_plots <- visualizeGAM(BW37_final, BW37_pred, "Gray's and Strap-toothed Whale")
# note: effect of ice regime not visualized for BW37 because it's not a smooth term
BW29_plots <- visualizeGAM(BW29_final, BW29_pred, 'Southern Bottlenose Whale')

