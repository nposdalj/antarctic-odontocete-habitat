library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)

# ------------- Step 0: Choose Species ----------------
# Modeling Pm for all sites it is present (EI & CI), 40 km radius environmental data
species <- c('BW37')
# Note: not enough data to model BW58 at any site
sites <- c('EI','KGI','CI') # sites where BW37 can be modeled

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
# only adding lags for chlorophyll and primary production, change that later if needed
#   other variables with lags available: temperature, salinity, EKE
depths <- c(0, 67, 920) # depths for BW37
sp_specific <- allData %>% subset(select=-c(BW29,BW58,Oo,Pm,Gm)) %>%
  subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                  temperature_0,salinity_0,EKE_0,temperature_67,salinity_67,EKE_67,
                  temperature_920,salinity_920,EKE_920, chla_0,o2_0,productivity_0,chla_67,
                  o2_67,productivity_67, chla_920,o2_920,productivity_920,
                  ssh_sd, mixed_layer_sd, fsle_sd, temp_sd_0, salinity_sd_0, EKE_mad_0, 
                  chla_sd_0,o2_sd_0,productivity_sd_0,ice_regime,fsle_orient, chla_1mon,
                  chla_2mon, chla_3mon, chla_4mon, chla_5mon, chla_6mon, productivity_1mon,
                  productivity_2mon, productivity_3mon, productivity_4mon, 
                  productivity_5mon, productivity_6mon))


# ------------- Step 2: Average by ACF ------------
acf_table <- read.csv("/Users/trisha/scripps/antarctic-odontocete-habitat/Autocorrelation/acf_table.csv")
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
  
  # Taking the average of environnmental variables at surface and at 920 m depth
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
                         o2_0_sd = mean_col('o2_sd_0'),
                         salinity_0_sd = mean_col('salinity_sd_0'),temperature_920 = mean_col('temperature_920'),
                         salinity_920 = mean_col('salinity_920'), EKE_920 = mean_col('EKE_920'),
                         o2_920 = mean_col('o2_920'), fsle_orient= mean_col('fsle_orient'), 
                         temperature_67 = mean_col('temperature_67'),
                         salinity_67 = mean_col('salinity_67'), EKE_67 = mean_col('EKE_67'),
                         o2_67 = mean_col('o2_67'))
  
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


# ------------- Step 4: Find Relevant Lags ---------------------
# first: ACF bin lagged data
binLagged <- function(site, bin) {
  # Filtering by site and adding bin_start date
  sp_filtered <- sp_specific %>% filter(Site==site) %>%
    mutate(bin_start = floor_date(date, unit = paste(bin, 'days')))
  
  # Function for taking mean with correct syntax for a column
  mean_col <- function(col) expr(mean(!!sym(col), na.rm = TRUE))
  
  # Expression to summarize species presence
  species_expr <- set_names(list(expr(as.integer(any(!!sym(species) == 1)))), species)
  
  # Taking the average of lagged variables at surface
  summarize_cols <- list(chla = mean_col('chla_0'), productivity = mean_col('productivity_0'),
                         chla_1mon = mean_col('chla_1mon'), chla_2mon  = mean_col('chla_2mon'), 
                         chla_3mon = mean_col('chla_3mon'), chla_4mon = mean_col('chla_4mon'), 
                         chla_5mon = mean_col('chla_5mon'), chla_6mon = mean_col('chla_6mon'), 
                         productivity_1mon = mean_col('productivity_1mon'),
                         productivity_2mon = mean_col('productivity_2mon'), 
                         productivity_3mon = mean_col('productivity_3mon'), 
                         productivity_4mon = mean_col('productivity_4mon'), 
                         productivity_5mon = mean_col('productivity_5mon'), 
                         productivity_6mon = mean_col('productivity_6mon'))
  
  # Joining together dataframes and grouping to make final binned data
  all_summaries <- c(species_expr, summarize_cols)
  sp_binned <- sp_filtered %>%
    group_by(bin_start, Site) %>%
    summarise(!!!all_summaries, .groups = "drop")
  
  return(sp_binned)
}
EI_lagged <- binLagged('EI',EI_acf)
KGI_lagged <- binLagged('KGI',KGI_acf)
CI_lagged <- binLagged('CI',CI_acf)

# ELEPHANT ISLAND
# ...PRIMARY PRODUCTION...
# 1 month
EI_lag <- gam(BW37 ~ s(productivity_1mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value: 0.247   deviance explained: 7.96%  AIC: 66.54915

# 2 month
EI_lag <- gam(BW37 ~ s(productivity_2mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value: 0.748   deviance explained: 0.16%  AIC: 68.51961

# 3 month
EI_lag <- gam(BW37 ~ s(productivity_3mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value: 0.466   deviance explained: 3.46%  AIC: 68.71169

# 4 month
EI_lag <- gam(BW37 ~ s(productivity_4mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value: 0.657   deviance explained: 0.307%  AIC: 68.42453

# 5 month
EI_lag <- gam(BW37 ~ s(productivity_5mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value: 0.0361   deviance explained: 7.51%  AIC: 63.76936

# 6 month
EI_lag <- gam(BW37 ~ s(productivity_6mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value: 0.207   deviance explained: 2.56%  AIC: 66.9691

# ...CHLOROPHYLL...
# 1 month
EI_lag <- gam(BW37 ~ s(chla_1mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value:  0.509  deviance explained: 4.66%  AIC: 68.27719

# 2 month
EI_lag <- gam(BW37 ~ s(chla_2mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value:  0.953  deviance explained: 0.0056%  AIC: 68.61948

# 3 month
EI_lag <- gam(BW37 ~ s(chla_3mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value:  0.119  deviance explained: 10.2%  AIC: 65.07419

# 4 month
EI_lag <- gam(BW37 ~ s(chla_4mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value:  0.64  deviance explained: 2.28%  AIC: 68.81738

# 5 month
EI_lag <- gam(BW37 ~ s(chla_5mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value:  0.111  deviance explained: 12.9%  AIC: 64.16685

# 6 month
EI_lag <- gam(BW37 ~ s(chla_6mon,k=4), 
              family = binomial, data = EI_lagged, method = 'REML')
# p-value:  0.259  deviance explained: 2.05%  AIC: 67.301






# KING GEORGE ISLAND
# ...PRIMARY PRODUCTION...
# no lag
KGI_lag <- gam(BW37 ~ s(productivity,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: .0988    deviance explained: 10.6%   AIC: 32.39145

# 1 month
KGI_lag <- gam(BW37 ~ s(productivity_1mon,k=4), 
              family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.0463    deviance explained: 60.3%  AIC: 20.45698

# 2 month
KGI_lag <- gam(BW37 ~ s(productivity_2mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.127    deviance explained: 66.9%  AIC: 17.9588

# 3 month
KGI_lag <- gam(BW37 ~ s(productivity_3mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.0313    deviance explained: 23.5%  AIC: 28.29084

# 4 month
KGI_lag <- gam(BW37 ~ s(productivity_4mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.429    deviance explained: 58.9%  AIC: 19.25084

# 5 month
KGI_lag <- gam(BW37 ~ s(productivity_5mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.358    deviance explained: 0.358%  AIC: 34.854

# 6 month
KGI_lag <- gam(BW37 ~ s(productivity_6mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.353    deviance explained: 2.79%  AIC: 34.86869

# ...CHLOROPHYLL...
# no lag
KGI_lag <- gam(BW37 ~ s(chla,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.129   deviance explained: 15.8%  AIC: 31.75674

# 1 month
KGI_lag <- gam(BW37 ~ s(chla_1mon,k=4), 
              family = binomial, data = KGI_lagged, method = 'REML')
# did not converge

# 2 month
KGI_lag <- gam(BW37 ~ s(chla_2mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.031   deviance explained: 46.8%  AIC: 24.63693

# 3 month
KGI_lag <- gam(BW37 ~ s(chla_3mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.147   deviance explained: 34.3%  AIC: 28.76052

# 4 month
KGI_lag <- gam(BW37 ~ s(chla_4mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.597   deviance explained: 76.8%  AIC: 13.65604

# 5 month
KGI_lag <- gam(BW37 ~ s(chla_5mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.246   deviance explained: 38.3%  AIC: 26.27662

# 6 month
KGI_lag <- gam(BW37 ~ s(chla_6mon,k=4), 
               family = binomial, data = KGI_lagged, method = 'REML')
# p-value: 0.102   deviance explained: 40.4%  AIC: 26.54112




# CLARENCE ISLAND
# ...PRIMARY PRODUCTION...
# no lag
CI_lag <- gam(BW37 ~ s(productivity,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.113   deviance explained: 16.6%   AIC: 46.93949

# 1 month
CI_lag <- gam(BW37 ~ s(productivity_1mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.132   deviance explained: 19.9%   AIC: 45.43269

# 2 month
CI_lag <- gam(BW37 ~ s(productivity_2mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.115   deviance explained: 12.6%   AIC: 47.99688

# 3 month
CI_lag <- gam(BW37 ~ s(productivity_3mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.155   deviance explained: 4.21%   AIC: 49.64837

# 4 month
CI_lag <- gam(BW37 ~ s(productivity_4mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.385   deviance explained: 1.56%   AIC: 50.91207

# 5 month
CI_lag <- gam(BW37 ~ s(productivity_5mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.567   deviance explained: 0.68%   AIC: 51.33258

# 6 month
CI_lag <- gam(BW37 ~ s(productivity_6mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.114   deviance explained: 7%   AIC: 48.32256

# ...CHLOROPHYLL...
# no lag
CI_lag <- gam(BW37 ~ s(chla,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value:  0.272  deviance explained: 27%   AIC: 41.91113

# 1 month
CI_lag <- gam(BW37 ~ s(chla_1mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.132    deviance explained: 25.5%   AIC: 43.37097

# 2 month
CI_lag <- gam(BW37 ~ s(chla_2mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.153    deviance explained: 12.3%   AIC: 47.66694

# 3 month
CI_lag <- gam(BW37 ~ s(chla_3mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.0141    deviance explained: 33.6%   AIC: 39.54774

# 4 month
CI_lag <- gam(BW37 ~ s(chla_4mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.171    deviance explained: 3.74%   AIC: 49.87255

# 5 month
CI_lag <- gam(BW37 ~ s(chla_5mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.354    deviance explained: 4.74%   AIC: 50.68253

# 6 month
CI_lag <- gam(BW37 ~ s(chla_6mon,k=4), 
              family = binomial, data = CI_lagged, method = 'REML')
# p-value: 0.32    deviance explained: 2.4%   AIC: 50.51354





# -------------- Step 5: VIF for Correlation -------------------
# Including variables at depth
# Also, not modeling AAO (varies on yearly timescales) and ice thickness
# Including julian day only for KGI because it covers almost a whole year
# Ice regime not included in VIF but will be present in modeling step
# including just FSLE magnitude in vif analysis, but it will be interaction term with fsle orientation in GAM

# ELEPHANT ISLAND
EI_pred <- c()
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)

# KING GEORGE ISLAND
KGI_pred <- c()
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)

# CLARENCE ISLAND
CI_pred <- c()
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)



# -------------- Step 6: Build GAMs ------------------------
# Function to visualize GAMs on a probability scale with the proper confidence interval
# Run this for each iteration of the model to plot smooth terms
plotGam <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],scheme=2,seWithMean=TRUE,pages=1))
}

# -------------------- Step 6a: Elephant Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Similar process: weighing to reach roughly 1:1 1s to 0s and setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in EI_pred include



# SINGLE VARIABLE GAMS

# MULTIPLE VARIABLES

# REFINING MODEL


# -------------------- Step 6a: King George Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Similar process: weighing to reach roughly 1:1 1s to 0s and setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in KGI_pred include



# SINGLE VARIABLE GAMS

# MULTIPLE VARIABLES

# REFINING MODEL

# -------------------- Step 6a: Clarence Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Similar process: weighing to reach roughly 1:1 1s to 0s and setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in CI_pred include



# SINGLE VARIABLE GAMS

# MULTIPLE VARIABLES

# REFINING MODEL