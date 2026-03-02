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
allData <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/allData.csv")
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
acf_table <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Autocorrelation/acf_table.csv")
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
  # replace this:
  # species_expr <- set_names(list(expr(as.integer(any(!!sym(species) > 0)))), species)
  
  # with this: mean of the raw species values per bin
  species_expr <- set_names(
    list(expr(mean(!!sym(species), na.rm = TRUE))),
    species
  )
  
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
                     "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- reformulate(KGI_pred, response = species)
KGI_vif <- lm(mod_formula, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0  temperature_0          EKE_0 
# 1.732701       5.564749       4.053836      11.374234       1.245090      11.944795      17.466446       1.106210 
# o2_0         chla_0 productivity_0 
# 7.529085      12.121116      13.377388 

# removing temp_0
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
                    "EKE_0",'o2_0','chla_0','productivity_0')
mod_formula <- reformulate(KGI_pred, response = species)
KGI_vif <- lm(mod_formula, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0          EKE_0           o2_0 
# 1.731489       5.431432       3.840026       9.240507       1.209360       9.878726       1.105461       2.652091 
# chla_0 productivity_0 
# 10.893597       7.337584

# removing chlorophyll
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
              "EKE_0",'o2_0','productivity_0')
mod_formula <- reformulate(KGI_pred, response = species)
KGI_vif <- lm(mod_formula, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0          EKE_0           o2_0 
# 1.589022       4.536424       3.614101       7.471798       1.148549       9.677623       1.097909       2.461717 
# productivity_0 
# 3.339818 

# removing salinity
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', 
              "EKE_0",'o2_0','productivity_0')
mod_formula <- reformulate(KGI_pred, response = species)
KGI_vif <- lm(mod_formula, data = KGI_binned)
vif(KGI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient          EKE_0           o2_0 productivity_0 
# 1.585085       3.518265       2.235293       4.730680       1.148037       1.090434       2.019236       3.154052 

# Final predictors for KGI: FLSE, SSH, mixed_layer, ince_conc, fsle_orient, EKE_0, 02_0, productivity_0

# CLARENCE ISLAND
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
                     "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0','EKE_0_mad')
mod_formula <- reformulate(CI_pred, response = species)
CI_vif <- lm(mod_formula, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0  temperature_0          EKE_0 
# 3.592319       4.811339       3.492701       2.532728       1.487290       6.785843      12.015794       1.129707 
# o2_0         chla_0 productivity_0      EKE_0_mad 
# 5.600827       5.351261      13.132636       1.401548 

# removing SST
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', "salinity_0", 
                     "EKE_0",'o2_0','chla_0','productivity_0','EKE_0_mad')
mod_formula <- reformulate(CI_pred, response = species)
CI_vif <- lm(mod_formula, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient     salinity_0          EKE_0           o2_0 
# 3.286721       4.685039       3.487291       2.062552       1.396326       6.654851       1.126196       3.656945 
# chla_0 productivity_0      EKE_0_mad 
# 4.740739       5.089223       1.379581 

# removing salinity
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 'fsle_orient', 
             "EKE_0",'o2_0','chla_0','productivity_0','EKE_0_mad')
mod_formula <- reformulate(CI_pred, response = species)
CI_vif <- lm(mod_formula, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc    fsle_orient          EKE_0           o2_0         chla_0 
# 3.039547       2.011872       2.087581       2.023663       1.209645       1.122853       3.521689       4.737174 
# productivity_0      EKE_0_mad 
# 4.794004       1.379580 

# Final Clarence Island predictors: FSLE, SSH, mixed layer thickness, sea ice concentration,
#   FSLE orientation, EKE, oxygen concentration, Chla, net primary production, EKE_0_mad

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
#Full model
EI_gam <- gam(Gm ~ s(salinity_0,k=4) + s(mixed_layer,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) +
                s(EKE_0,k=4,sp=0.1) + s(ice_conc,k=4) + s(FSLE,k=4,sp=0.1)
              , family = tw(link = "log", a = 1.1, b = 1.9), data = EI_binned, method = "REML")

# Remove FSLE
EI_gam <- gam(Gm ~ s(salinity_0,k=4) + s(mixed_layer,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) +
                s(EKE_0,k=4,sp=0.1) + s(ice_conc,k=4)
              , family = tw(link = "log", a = 1.1, b = 1.9), data = EI_binned, method = "REML")

# Remove EKE
EI_gam <- gam(Gm ~ s(salinity_0,k=4) + s(mixed_layer,k=4,sp=0.1) + s(SSH,k=4,sp=0.1) +
                s(ice_conc,k=4)
              , family = tw(link = "log", a = 1.1, b = 1.9), data = EI_binned, method = "REML")

# Remove SSH
EI_gam <- gam(Gm ~ s(salinity_0,k=4) + s(mixed_layer,k=4,sp=0.1)
              , family = tw(link = "log", a = 1.1, b = 1.9), data = EI_binned, method = "REML")


# final model: mixed layer thickness and salinity
#   could also have to do with weighing - EI did not need to have detections weighed at all to reach 1:1 ratio
EI_final <- gam(Gm ~ s(salinity_0,k=4) + s(mixed_layer,k=4,sp=0.1)
                , family = tw(link = "log", a = 1.1, b = 1.9), data = EI_binned, method = "REML")


# -------------------- Step 5b: King George Island GAM -------------------------
# initial gam based on list of KGI predictors: FLSE, SSH, mixed_layer, ince_conc, fsle_orient, EKE_0, 02_0, productivity_0
KGI_gam <- gam(Gm ~ s(FSLE,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) + s(ice_conc,k=4) + 
                  s(EKE_0,k=4) + s(o2_0,k=4) + s(productivity_0,k=4), 
               family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

# Remove 02
KGI_gam <- gam(Gm ~ s(FSLE,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) + s(ice_conc,k=4) + 
                 s(EKE_0,k=4) + s(productivity_0,k=4),
               family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

# Remove Ice concentration
KGI_gam <- gam(Gm ~ s(FSLE,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) +
                 s(EKE_0,k=4) + s(productivity_0,k=4), 
               family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

#Remove EKE
KGI_gam <- gam(Gm ~ s(FSLE,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) + s(productivity_0,k=4), 
               family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

#Remove FSLE
KGI_gam <- gam(Gm ~ s(SSH,k=4) + s(mixed_layer,k=4) + s(productivity_0,k=4), 
               family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

#Remove ML
KGI_final <- gam(Gm ~ s(SSH,k=4) + s(productivity_0,k=4), 
               family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")


# --------------------- Step 5c: Clarence Island GAM ------------------------------------
# List of predictors are stored in CI_pred, include FSLE, SSH, mixed layer thickness, sea ice concentration,
#   FSLE orientation, EKE, oxygen concentration, Chla, net primary production, EKE_0_mad

CI_gam <- gam(Gm ~ s(FSLE,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) +
                s(ice_conc,k=4) + s(EKE_0_mad,k=4) + s(o2_0,k=4) +
                s(productivity_0,k=4),
              family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

# Remove 02
CI_gam <- gam(Gm ~ s(FSLE,k=4) + s(SSH,k=4) + s(mixed_layer,k=4) +
                s(ice_conc,k=4) + s(EKE_0_mad,k=4) + s(productivity_0,k=4),
              family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

#Remove FSLE
CI_gam <- gam(Gm ~ s(SSH,k=4) + s(mixed_layer,k=4) +
                s(ice_conc,k=4) + s(EKE_0_mad,k=4) + s(productivity_0,k=4),
              family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

#Remove Mixed Layer
CI_gam <- gam(Gm ~ s(SSH,k=4) + s(ice_conc,k=4) + s(EKE_0_mad,k=4) + s(productivity_0,k=4),
              family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

#Remove EKE
CI_gam <- gam(Gm ~ s(SSH,k=4) + s(ice_conc,k=4) + s(productivity_0,k=4),
              family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

#Remove Ice
CI_final <- gam(Gm ~ s(SSH,k=6) + s(productivity_0,k=6),
              family = tw(link = "log", a = 1.1, b = 1.9), data = KGI_binned, method = "REML")

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
KGI_pred <- c('SSH','productivity_0')
KGI_plots <- visualizeGAM(KGI_final, KGI_pred, 'KGI')

EI_pred <- c('mixed_layer','salinity_0')
EI_plots <- visualizeGAM(EI_final, EI_pred, 'EI')


CI_pred <- c('SSH','productivity_0')
CI_plots <- visualizeGAM(CI_final, CI_pred, 'CI')
