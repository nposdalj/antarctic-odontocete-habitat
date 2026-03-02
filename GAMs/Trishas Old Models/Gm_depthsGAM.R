library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
library(gratia)
library(patchwork)

# ------------- Step 0: Choose Species ----------------
# Modeling Pm for all sites it is present (EI & CI), 40 km radius environmental data
species <- c('Gm')
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
depths <- c(0, 16, 635) # Gm depths
sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Oo,Pm)) %>%
  subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                  temperature_0,salinity_0,EKE_0,temperature_16,salinity_16,EKE_16,
                  temperature_635,salinity_635,EKE_635, chla_0,o2_0,productivity_0,chla_16,
                  o2_16,productivity_16, chla_635,o2_635,productivity_635,
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

# making lagged GAMs to pick best lags for chlorophyll and primary production

# ELEPHANT ISLAND
# ...PRIMARY PRODUCTION...
# no lag
EI_lag <- gam(Gm ~ s(productivity,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.00291    deviance explained: 41.9%   AIC: 27.8412

# 1 month
EI_lag <- gam(Gm ~ s(productivity_1mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.00336    deviance explained: 33.4%   AIC: 31.35764

# 2 month
EI_lag <- gam(Gm ~ s(productivity_2mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.0112    deviance explained: 28%   AIC: 33.55119

# 3 month
EI_lag <- gam(Gm ~ s(productivity_3mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.00974    deviance explained: 45.9%   AIC: 26.2287

# 4 month
EI_lag <- gam(Gm ~ s(productivity_4mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.0875    deviance explained: 20.9%   AIC: 39.1987

# 5 month
EI_lag <- gam(Gm ~ s(productivity_5mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.0615    deviance explained: 41.9%   AIC: 31.68455

# 6 month
EI_lag <- gam(Gm ~ s(productivity_6mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.0116    deviance explained: 26.6%   AIC: 34.14929

# ...CHLOROPHYLL...
# no lag
EI_lag <- gam(Gm ~ s(chla,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.00234   deviance explained: 35.8%   AIC: 30.36378

# 1 month
EI_lag <- gam(Gm ~ s(chla_1mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.137   deviance explained: 37.7%   AIC: 33.44175

# 2 month
EI_lag <- gam(Gm ~ s(chla_2mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.386   deviance explained: 8.69%   AIC: 44.23984

# 3 month
EI_lag <- gam(Gm ~ s(chla_3mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: .0694   deviance explained: 20.6%   AIC: 39.1904

# 4 month
EI_lag <- gam(Gm ~ s(chla_4mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.151   deviance explained: 39%   AIC: 31.46814

# 5 month
EI_lag <- gam(Gm ~ s(chla_5mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.00325   deviance explained: 38%   AIC: 29.47171

# 6 month
EI_lag <- gam(Gm ~ s(chla_6mon,k=4), family = binomial, method = 'REML', data = EI_lagged)
# p-value: 0.251   deviance explained: 10.5%   AIC: 43.46904




# KING GEORGE ISLAND
# ...PRIMARY PRODUCTION...
# no lag
KGI_lag <- gam(Gm ~ s(productivity,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.00253   deviance explained: 9.12%   AIC: 142.5871

# 1 month
KGI_lag <- gam(Gm ~ s(productivity_1mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.147   deviance explained: 5.19%   AIC: 151.4417

# 2 month
KGI_lag <- gam(Gm ~ s(productivity_2mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.000179   deviance explained: 17.3%   AIC: 134.0152

# 3 month
KGI_lag <- gam(Gm ~ s(productivity_3mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.000355   deviance explained: 19.7%   AIC: 129.7531

# 4 month
KGI_lag <- gam(Gm ~ s(productivity_4mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.000169   deviance explained: 21.8%   AIC: 126.7335

# 5 month
KGI_lag <- gam(Gm ~ s(productivity_5mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 2.16e-05   deviance explained: 19.6%   AIC: 128.2386

# 5 month
KGI_lag <- gam(Gm ~ s(productivity_5mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 2.16e-05   deviance explained: 19.6%   AIC: 128.2386

# 6 month
KGI_lag <- gam(Gm ~ s(productivity_6mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.000282  deviance explained: 14.4%   AIC: 137.0094

# ...CHLOROPHYLL...
# no lag
KGI_lag <- gam(Gm ~ s(chla,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.126   deviance explained: 6.67%    AIC: 148.9245

# 1 month
KGI_lag <- gam(Gm ~ s(chla_1mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.118   deviance explained: 14.3%    AIC: 138.6544

# 2 month
KGI_lag <- gam(Gm ~ s(chla_2mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.000305   deviance explained: 25.5%    AIC: 121.619

# 3 month
KGI_lag <- gam(Gm ~ s(chla_3mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.00355   deviance explained: 18.1%    AIC: 132.4659

# 4 month
KGI_lag <- gam(Gm ~ s(chla_4mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.00095   deviance explained: 16.6%    AIC: 134.7946

# 5 month
KGI_lag <- gam(Gm ~ s(chla_5mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 0.0118   deviance explained: 8.86%    AIC: 145.5996

# 6 month
KGI_lag <- gam(Gm ~ s(chla_6mon,k=4), family = binomial, method = 'REML', data = KGI_lagged)
# p-value: 1.62e-05   deviance explained: 21.7%    AIC: 127.3692




# CLARENCE ISLAND
# ...PRIMARY PRODUCTION...
# no lag
CI_lag <- gam(Gm ~ s(productivity,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.0302   deviance explained: 8.66%   AIC: 89.87685

# 1 month
CI_lag <- gam(Gm ~ s(productivity_1mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.048   deviance explained: 23.3%   AIC: 80.02578

# 2 month
CI_lag <- gam(Gm ~ s(productivity_2mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.0581   deviance explained: 12.2%   AIC: 90.18451

# 3 month
CI_lag <- gam(Gm ~ s(productivity_3mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.259   deviance explained: 7.57%   AIC: 93.56982

# 4 month
CI_lag <- gam(Gm ~ s(productivity_4mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.284   deviance explained: 5.17%   AIC: 96.19852

# 5 month
CI_lag <- gam(Gm ~ s(productivity_5mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.053   deviance explained: 12.8%   AIC: 89.06576

# 6 month
CI_lag <- gam(Gm ~ s(productivity_6mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.0333   deviance explained: 10.6%   AIC: 91.5196

# ...CHLOROPHYLL...
# no lag
CI_lag <- gam(Gm ~ s(chla,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.00954   deviance explained: 12.9%   AIC: 85.91753

# 1 month
CI_lag <- gam(Gm ~ s(chla_1mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.00756   deviance explained: 13%   AIC: 85.79459

# 2 month
CI_lag <- gam(Gm ~ s(chla_2mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.0727   deviance explained: 10.7%   AIC: 91.81761

# 3 month
CI_lag <- gam(Gm ~ s(chla_3mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.591   deviance explained: 2.09%   AIC: 98.372

# 4 month
CI_lag <- gam(Gm ~ s(chla_4mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.0889   deviance explained: 9.86%   AIC: 92.62715

# 5 month
CI_lag <- gam(Gm ~ s(chla_5mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.186   deviance explained: 22.5%   AIC: 79.06849

# 6 month
CI_lag <- gam(Gm ~ s(chla_6mon,k=4), family = binomial, method = 'REML', data = CI_lagged)
# p-value: 0.395   deviance explained: 3.97%   AIC: 97.30111


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



# -------------------- Step 6b: King George Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Similar process: weighing to reach roughly 1:1 1s to 0s and setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in KGI_pred include



# SINGLE VARIABLE GAMS

# MULTIPLE VARIABLES

# REFINING MODEL

# -------------------- Step 6b: Clarence Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Similar process: weighing to reach roughly 1:1 1s to 0s and setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in CI_pred include



# SINGLE VARIABLE GAMS

# MULTIPLE VARIABLES

# REFINING MODEL