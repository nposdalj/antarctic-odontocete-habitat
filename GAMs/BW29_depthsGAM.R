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
allData <- read.csv("/Users/trisha/scripps/antarctic-odontocete-habitat/Data/allData_40km.csv")
allData <- allData %>% subset(select=-X)
allData$date <- as.Date(allData$date, "%Y-%m-%d")
# Filter by species relevant data
# Only adding standard deviations of surface variables, feel free to change that if needed
if(species =='BW29') {
  depths <- c(0, 768)
  sp_specific <- allData %>% subset(select=-c(BW37,BW58,Oo,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_768,salinity_768,EKE_768,EKE_mad_768,
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
                         salinity_0_sd = mean_col('salinity_sd_0'),temperature_768 = mean_col('temperature_768'),
                         salinity_768 = mean_col('salinity_768'), EKE_mad_768 = mean_col('EKE_mad_768'),
                         o2_768 = mean_col('o2_768'), fsle_orient= mean_col('fsle_orient'))
  
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


# -------------- Step 4: VIF for Correlation -------------------
# Including variables at depth
# Also, not modeling AAO (varies on yearly timescales) and ice thickness
# Including julian day only for KGI because it covers almost a whole year
# Ice regime not included in VIF but will be present in modeling step
# including just FSLE magnitude in vif analysis, but it will be interaction term with fsle orientation in GAM

# ELEPHANT ISLAND
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
             "temperature_0", 'EKE_mad_0' ,'o2_0','chla_0','productivity_0',
             'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0   temperature_0       EKE_mad_0            o2_0 
# 4.138371        8.999121        8.664663        8.465678        6.303774       32.991087        2.017693       29.941611 
# chla_0  productivity_0    salinity_768 temperature_768     EKE_mad_768          o2_768 
# 12.574486       21.451406       13.180761       18.481152        1.857527        2.567427 

# removing surface temp
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
              'EKE_mad_0' ,'o2_0','chla_0','productivity_0',
             'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0       EKE_mad_0            o2_0          chla_0 
# 3.920737        8.876270        8.155184        7.947639        6.273319        2.024743       24.940506        9.909413 
# productivity_0    salinity_768 temperature_768     EKE_mad_768          o2_768 
# 20.002890       13.243761       17.559802        1.735246        2.194862 

# removing surface primary production
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
             'EKE_mad_0' ,'o2_0','chla_0',
             'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0       EKE_mad_0            o2_0          chla_0 
# 3.813020        8.853095        6.122185        6.566861        6.224487        1.856676       11.467889        9.290520 
# salinity_768 temperature_768     EKE_mad_768          o2_768 
# 13.129956       15.754712        1.704804        1.752906 

# removing 768m temperature
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
             'EKE_mad_0' ,'o2_0','chla_0',
             'salinity_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE          SSH  mixed_layer     ice_conc   salinity_0    EKE_mad_0         o2_0       chla_0 salinity_768  EKE_mad_768       o2_768 
# 3.813415     2.255628     5.201535     4.814262     5.440616     1.726575     6.528614     9.057792     9.427520     1.642745     1.718622 

# removing 768m salinity
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
             'EKE_mad_0' ,'o2_0','chla_0',
             'EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE         SSH mixed_layer    ice_conc  salinity_0   EKE_mad_0        o2_0      chla_0 EKE_mad_768      o2_768 
# 1.891657    2.108126    4.285818    3.390141    3.888680    1.460162    6.387231    9.797013    1.338864    1.695447 

# removing chlorophyll
EI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
             'EKE_mad_0' ,'o2_0',
             'EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(EI_pred, collapse = " + "))
EI_vif <- glm(as.formula(mod_formula),family=binomial, data = EI_binned)
vif(EI_vif)
# FSLE         SSH mixed_layer    ice_conc  salinity_0   EKE_mad_0        o2_0 EKE_mad_768      o2_768 
# 1.684850    1.954120    4.260477    2.511020    2.447548    1.441115    4.807040    1.342076    1.695196 

# GAM predictors: FSLE, SSH, mixed layer depth, sea ice concentration, surface salinity, surface EKE, surface oxygen
#   768 m EKE, 768 m oxygen




# KING GEORGE ISLAND
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
             "temperature_0", 'EKE_mad_0' ,'o2_0','chla_0','productivity_0',
             'salinity_768','temperature_768','EKE_mad_768','o2_768','julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0   temperature_0       EKE_mad_0            o2_0 
# 4.228240        2.308024       12.208235        3.160563        8.083785       78.979122        2.595131       19.237562 
# chla_0  productivity_0    salinity_768 temperature_768     EKE_mad_768          o2_768      julian_day 
# 10.095445       23.427621        4.972203        8.380225        1.759170        3.382552       21.322022 

# dropping surface temperature
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
              'EKE_mad_0' ,'o2_0','chla_0','productivity_0',
              'salinity_768','temperature_768','EKE_mad_768','o2_768','julian_day')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0       EKE_mad_0            o2_0          chla_0 
# 3.406110        2.179548        6.319869        3.257645        7.788995        2.482338       13.415993       11.494913 
# productivity_0    salinity_768 temperature_768     EKE_mad_768          o2_768      julian_day 
# 19.998231        4.289120        8.681848        1.777048        3.367259       26.176067 

# dropping julian day
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
              'EKE_mad_0' ,'o2_0','chla_0','productivity_0',
              'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0       EKE_mad_0            o2_0          chla_0 
# 3.469999        2.104196        6.081926        3.161453        7.536540        2.388979        6.292098        9.803129 
# productivity_0    salinity_768 temperature_768     EKE_mad_768          o2_768 
# 17.700172        3.743951        7.754695        1.767769        3.094647 

# dropping primary production
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
              'EKE_mad_0' ,'o2_0','chla_0',
              'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0       EKE_mad_0            o2_0          chla_0 
# 3.579525        2.130970        5.623761        2.739497        7.090765        2.321739        4.420736        4.264790 
# salinity_768 temperature_768     EKE_mad_768          o2_768 
# 2.885451        4.804884        1.744649        2.992620 

# dropping surface salinity
KGI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 
              'EKE_mad_0' ,'o2_0','chla_0',
              'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(KGI_pred, collapse = " + "))
KGI_vif <- glm(as.formula(mod_formula),family=binomial, data = KGI_binned)
vif(KGI_vif)
# FSLE             SSH     mixed_layer        ice_conc       EKE_mad_0            o2_0          chla_0    salinity_768 
# 2.163992        2.470691        3.101607        2.813127        1.775634        2.882169        2.817173        2.829212 
# temperature_768     EKE_mad_768          o2_768 
# 4.064711        1.835048        1.940621 

# GAM predictors KGI: FSLE, SSH, mixed layer depth, ice concentration, surface EKE, surface oxygen, chlorophyll, 768m salinity, 
#    768m temperature, 768m salinity, 768m EKE, 768m oxygen




# CLARENCE ISLAND
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
              "temperature_0", 'EKE_mad_0' ,'o2_0','chla_0','productivity_0',
              'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0   temperature_0       EKE_mad_0            o2_0 
# 79.15717       226.42688       321.96831        14.78656      2850.62111      1430.54186        23.59628       641.53011 
# chla_0  productivity_0    salinity_768 temperature_768     EKE_mad_768          o2_768 
# 375.31740       681.58620       274.23184       356.43470        37.51588        26.74891 

# removing surface temperature
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
             'EKE_mad_0' ,'o2_0','chla_0','productivity_0',
             'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0       EKE_mad_0            o2_0          chla_0 
# 57.44217        36.38167        99.70812        10.09728       228.76166        21.35414       667.55090        50.46854 
# productivity_0    salinity_768 temperature_768     EKE_mad_768          o2_768 
# 120.12915       189.77344       237.47447        34.71987        16.86531

# removing surface oxygen
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", "salinity_0", 
             'EKE_mad_0','chla_0','productivity_0',
             'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE             SSH     mixed_layer        ice_conc      salinity_0       EKE_mad_0          chla_0  productivity_0 
# 62.23012        33.50472        59.54789        16.10809       259.59428        23.69887        44.95023        86.38316 
# salinity_768 temperature_768     EKE_mad_768          o2_768 
# 50.57186        90.56652        17.32817        10.81857 

# removing surface salinity
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 
             'EKE_mad_0','chla_0','productivity_0',
             'salinity_768','temperature_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE             SSH     mixed_layer        ice_conc       EKE_mad_0          chla_0  productivity_0    salinity_768 
# 11.733537       11.476112       10.688905        3.028607       10.605493       10.406661       28.982548       12.656938 
# temperature_768     EKE_mad_768          o2_768 
# 44.065643        7.173388        8.293336 

# removing 768m temperature
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 
             'EKE_mad_0','chla_0','productivity_0',
             'salinity_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE            SSH    mixed_layer       ice_conc      EKE_mad_0         chla_0 productivity_0   salinity_768    EKE_mad_768 
# 5.699566       5.417158       5.285425       3.284533       4.714192       7.293697      11.452611       9.969015       2.351136 
# o2_768 
# 3.908949 

# removing primary production
CI_pred <- c("FSLE", "SSH", "mixed_layer", "ice_conc", 
             'EKE_mad_0','chla_0',
             'salinity_768','EKE_mad_768','o2_768')
mod_formula <- paste(species, "~", paste(CI_pred, collapse = " + "))
CI_vif <- glm(as.formula(mod_formula),family=binomial, data = CI_binned)
vif(CI_vif)
# FSLE          SSH  mixed_layer     ice_conc    EKE_mad_0       chla_0 salinity_768  EKE_mad_768       o2_768 
# 3.146149     5.721636     4.542309     3.221460     3.893474     2.477894     3.832778     1.481576     3.447434 

# GAM predictors for CI: FSLE, SSH, mixed layer depth, sea ice concentration, surface EKE, chlorophyll,
#   768m salinity, 768m EKE, 768m oxygen


# -------------- Step 5: Build GAMs ------------------------
# Function to visualize GAMs on a probability scale with the proper confidence interval
# Run this for each iteration of the model to plot smooth terms
plotGam <- function(gam) {
  return(plot(gam,trans=plogis,shift=coef(gam)[1],scheme=2,seWithMean=TRUE,pages=1))
}


# -------------------- Step 5a: Elephant Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Similar process: weighing to reach roughly 1:1 1s to 0s and setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in EI_pred include FSLE, SSH, mixed layer, ice concentration, 
#   surface salinity, surface EKE, surface oxygen, 768m EKE, 768m salinity, ice regime



# SINGLE VARIABLE GAMS

# not adding weights for EI because ratio of 1s already exceed 0s (41:28))
# starting with FSLE magnitude and orientation
EI_gam <- gam(BW29 ~ s(FSLE,fsle_orient,k=4), family=binomial, method='REML', data=EI_binned)
# Formula:
#   BW29 ~ s(FSLE, fsle_orient, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   0.4735     0.2726   1.737   0.0824 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(FSLE,fsle_orient) 2.315   2.53  8.336  0.0349 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0998   Deviance explained = 11.4%
# -REML = 42.571  Scale est. = 1         n = 69

# sea surface height
EI_gam <- gam(BW29 ~ s(SSH,k=4), family=binomial, method='REML', data=EI_binned)
# Formula:
#   BW29 ~ s(SSH, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3906     0.2496   1.565    0.118
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(SSH) 1.457  1.767  3.188   0.257
# 
# R-sq.(adj) =  0.026   Deviance explained = 3.43%
# -REML = 46.441  Scale est. = 1         n = 69

# mixed layer depth
EI_gam <- gam(BW29 ~ s(mixed_layer,k=4), family=binomial, method='REML', data=EI_binned)
# Formula:
#   BW29 ~ s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   0.5069     0.2787   1.819   0.0689 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(mixed_layer)   1      1  8.604 0.00335 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.134   Deviance explained = 13.1%
# -REML = 41.001  Scale est. = 1         n = 69

# sea ice concentration
EI_gam <- gam(BW29 ~ s(ice_conc,k=4), family=binomial, method='REML', data=EI_binned)
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   BW29 ~ s(ice_conc, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.4063     0.2520   1.612    0.107
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(ice_conc)   1      1  1.726   0.189
# 
# R-sq.(adj) =  0.0135   Deviance explained = 2.29%
# -REML = 46.234  Scale est. = 1         n = 69

# surface salinity
EI_gam <- gam(BW29 ~ s(salinity_0,k=4), family=binomial, method='REML', data=EI_binned)
# Formula:
#   BW29 ~ s(salinity_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3836     0.2465   1.556     0.12
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(salinity_0) 1.369  1.646  1.189   0.584
# 
# R-sq.(adj) =  0.00078   Deviance explained = 1.52%
# -REML = 47.234  Scale est. = 1         n = 69

# surface EKE deviation
EI_gam <- gam(BW29 ~ s(EKE_mad_0,k=4), family=binomial, method='REML', data=EI_binned)
# Formula:
#   BW29 ~ s(EKE_mad_0, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   0.4893     0.2732   1.791   0.0733 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(EKE_mad_0)   1      1  3.028  0.0819 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0483   Deviance explained =  5.3%
# -REML = 44.311  Scale est. = 1         n = 69

# 768m EKE deviation
EI_gam <- gam(BW29 ~ s(EKE_mad_768,k=4), family=binomial, method='REML', data=EI_binned)
# Formula:
#   BW29 ~ s(EKE_mad_768, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3821     0.2455   1.557     0.12
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(EKE_mad_768)   1      1   0.17    0.68
# 
# R-sq.(adj) =  -0.0123   Deviance explained = 0.182%
# -REML = 47.493  Scale est. = 1         n = 69

# 768m salinity
EI_gam <- gam(BW29 ~ s(salinity_768,k=4), family=binomial, method='REML', data=EI_binned)
# Formula:
#   BW29 ~ s(salinity_768, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.3828     0.2457   1.558    0.119
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(salinity_768) 1.318  1.565  0.575   0.758
# 
# R-sq.(adj) =  -0.00617   Deviance explained = 0.993%
# -REML = 47.451  Scale est. = 1         n = 69

# sea ice regime
EI_gam <- gam(BW29 ~ ice_regime, family=binomial, method='REML', data=EI_binned)
# Formula:
#   BW29 ~ ice_regime
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)             1.609      1.095   1.469    0.142
# ice_regimeincreasing    0.470      1.525   0.308    0.758
# ice_regimenone         -1.721      1.145  -1.503    0.133
# ice_regimestable       -1.157      1.197  -0.967    0.334
# 
# 
# R-sq.(adj) =  0.0586   Deviance explained = 8.21%
# -REML = 40.767  Scale est. = 1         n = 69

# individually significant predictors: FSLE and mixed layer depth

# MULTIVARIATE GAM
# adding FSLE and mixed layer depth
EI_gam <- gam(BW29 ~ s(mixed_layer,k=4) + s(FSLE,fsle_orient,k=4), 
              family=binomial, method='REML', data=EI_binned)
# AIC: 84.97974
# summary:
# BW29 ~ s(mixed_layer, k = 4) + s(FSLE, fsle_orient, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   0.6280     0.3191   1.968   0.0491 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# s(mixed_layer)        1      1  5.138  0.0234 *
#   s(FSLE,fsle_orient)   2      2  3.498  0.1739  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.151   Deviance explained = 17.4%
# -REML = 39.205  Scale est. = 1         n = 69
# .................................................
# concurvity and plots not problematic, but removing FSLE because not significant

EI_gam <- gam(BW29 ~ s(mixed_layer,k=4),
              family=binomial, method='REML', data=EI_binned)
# running gam.check
# model residuals look good, p-values good (not messing with k)

# final model: just keeping mixed layer depth
EI_final <- gam(BW29 ~ s(mixed_layer,k=4),
              family=binomial, method='REML', data=EI_binned)
# AIC: 85.02098
# Formula:
#   BW29 ~ s(mixed_layer, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   0.5069     0.2787   1.819   0.0689 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(mixed_layer)   1      1  8.604 0.00335 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.134   Deviance explained = 13.1%
# -REML = 41.001  Scale est. = 1         n = 69


# -------------------- Step 5b: King George Island GAM ------------------------------
# starting by building GAMs one predictor at a time to find significant variables
# Similar process: weighing to reach roughly 1:1 1s to 0s and setting initial knots at 4, smoothing at 0.1
# List of predictors are stored in KGI_pred include FSLE, SSH, mixed layer, ice concentration, 
#   surface EKE, surface oxygen, surface chlorophyll, 768m salinity, 768m temp, 
#   768m EKE, 768m o2, ice regime

# weighing data (curren ratio 157 0s to 24 1s)
KGI_binned$weights <- ifelse(KGI_binned$BW29 == 1,6,1)

# SINGLE VARIABLE GAMS
# FSLE, FSLE orientation
KGI_gam <- gam(BW29 ~ s(FSLE,fsle_orient,k=4),
               data = KGI_binned, weights = weights, family = 'binomial', method = 'REML')

# sea surface height
# mixed layer depth
# sea ice concentration
# surface EKE
# surface oxygen
# surface chlorophyll
# 768m salinity
# 768m temperature
# 768m EKE
# 768m oxygen
# ice regime
