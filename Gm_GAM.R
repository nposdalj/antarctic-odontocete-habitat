library(tidyverse)
library(mgcv)
library(car)
library(rlang)
library(gridExtra)
# ------------- Step 0: Choose Species ----------------
# Modeling Gm for all sites
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
allData <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/allData.csv")
allData <- allData %>% subset(select=-X)
allData$date <- as.Date(allData$date, "%Y-%m-%d")
# Filter by species relevant data
if(species =='BW29') {
  depths <- c(0, 763)
  sp_specific <- allData %>% subset(select=-c(BW37,BW58,Oo,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_763,salinity_763,EKE_763))
} else if(species =='BW37') {
  depths <- c(0, 66, 902) 
  sp_specific <- allData %>% subset(select=-c(BW29,BW58,Oo,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_66,salinity_66,EKE_66,
                    temperature_902,salinity_902,EKE_902))
} else if(species =='Oo') {
  depths <- c(0, 11, 454) 
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Pm,Gm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_11,salinity_11,EKE_11,
                    temperature_454,salinity_454,EKE_454))
} else if(species =='Pm') {
  depths <- c(0, 380, 1684)
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Oo,Gm)) %>%
    subset(select=c(date,julian_day,Site,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_380,salinity_380,EKE_380,
                    temperature_1684,salinity_1684,EKE_1684))
} else if(species =='Gm') {
  depths <- c(0, 16, 644) 
  sp_specific <- allData  %>% subset(select=-c(BW29,BW37,BW58,Oo,Pm)) %>%
    subset(select=c(date,Site,julian_day,get(species),AAO,SSH,mixed_layer,ice_conc,ice_thickness,ice_diff,FSLE,
                    temperature_0,salinity_0,EKE_0,temperature_16,salinity_16,EKE_16,
                    temperature_644,salinity_644,EKE_644))
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
    EKE_0 = mean_col("EKE_0"))
  
  # Adding shallow dive depth for applicable variables
  summarize_cols[[paste0("temperature_", depths[2])]] <- mean_col(paste0("temperature_", depths[2]))
  summarize_cols[[paste0("salinity_", depths[2])]] <- mean_col(paste0("salinity_", depths[2]))
  summarize_cols[[paste0("EKE_", depths[2])]] <- mean_col(paste0("EKE_", depths[2]))
  # Adding deep dive depth for species with multiple depths
  if (species != "BW29") {
    summarize_cols[[paste0("temperature_", depths[3])]] <- mean_col(paste0("temperature_", depths[3]))
    summarize_cols[[paste0("salinity_", depths[3])]] <- mean_col(paste0("salinity_", depths[3]))
    summarize_cols[[paste0("EKE_", depths[3])]] <- mean_col(paste0("EKE_", depths[3]))
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
# might need to modify to account for the fact that sp_binned is already binned while weekly is not
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
if(species != 'BW29') {
  intl_predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                  "salinity_0", "temperature_0", "EKE_0", 'julian_day',
                  paste0("salinity_", depths[2]), paste0("temperature_", depths[2]), paste0("EKE_", depths[2]),
                  paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
} else
  intl_predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                  "salinity_0", "temperature_0", "EKE_0", 'julian_day',
                  paste0("salinity_", depths[2]), paste0("temperature_", depths[2]), paste0("EKE_", depths[2]))
mod_formula <- paste(species, "~", paste(intl_predictors, collapse = " + "))

# ELEPHANT ISLAND
all_pred <- lm(as.formula(mod_formula), data = EI_binned)
vif(all_pred)
# all predictors included, VIF output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 4.927098        1.929727        9.773431       70.848542       36.443809       66.631590        9.523491 
# salinity_0   temperature_0           EKE_0      julian_day     salinity_16  temperature_16          EKE_16 
# 6236.209553    47755.793924       57.095117      186.515271     5535.993680    50447.750268       70.099264 
# salinity_644 temperature_644         EKE_644 
# 15.565532       21.584381       21.695226 

# Dropping mid-depth temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "salinity_0", "temperature_0", "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 4.880147        1.928873        8.657852       64.619713       36.363479       65.579715        9.515880 
# salinity_0   temperature_0           EKE_0      julian_day          EKE_16     salinity_16    salinity_644 
# 6166.734113      362.546303       50.354678      177.167519       57.064375     5481.752298       15.532845 
# temperature_644         EKE_644 
# 21.209879       19.916999 

# Dropping surface salinity
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "temperature_0", "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 2.989885        1.921956        8.096509       11.797814       22.689703       31.349292        4.894037 
# temperature_0           EKE_0      julian_day          EKE_16     salinity_16    salinity_644 temperature_644 
# 134.300135       42.124713      149.987357       56.047212       22.497066       15.518558       14.971091 
# EKE_644 
# 17.874919 

# Dropping julian day
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "temperature_0", "EKE_0", paste0("EKE_", depths[2]),
                paste0("salinity_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 2.912674        1.810639        7.556851       11.777619       20.169151       30.098091        4.849339 
# temperature_0           EKE_0          EKE_16     salinity_16    salinity_644 temperature_644         EKE_644 
# 60.156089       34.831710       50.766216       15.625792       11.550008       10.080809       17.300940 

# Dropping surface temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "EKE_0", paste0("EKE_", depths[2]),
                paste0("salinity_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 2.808513        1.756818        7.228430        7.601202       20.168923       23.018660        3.832377 
# EKE_0          EKE_16     salinity_16    salinity_644 temperature_644         EKE_644 
# 34.801781       50.229028        7.673175        8.158809        8.979628       17.032715 

# Dropping 16m EKE
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                paste0("salinity_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 2.562432        1.738579        6.870006        5.980210       19.057371       22.866227        2.908408 
# salinity_16    salinity_644 temperature_644         EKE_644 
# 5.138381        6.817948        7.998965        5.712491 

# Dropping ice thickness
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',
                paste0("salinity_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc        ice_diff     salinity_16 
# 2.559042        1.602485        6.513656        5.978821        4.267003        1.517004        5.064940 
# salinity_644 temperature_644         EKE_644 
# 6.799469        7.998319        5.510122 

# Dropping 644m temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',
                paste0("salinity_", depths[2]),
                paste0("salinity_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO         FSLE          SSH  mixed_layer     ice_conc     ice_diff  salinity_16 salinity_644      EKE_644 
# 2.079559     1.582135     3.151588     3.976107     3.507459     1.387602     3.173060     6.746640     2.283287

# Dropping 644m salinity
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',
                paste0("salinity_", depths[2]),
                paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO        FSLE         SSH mixed_layer    ice_conc    ice_diff salinity_16     EKE_644 
# 1.892168    1.580119    2.383937    2.219766    3.196145    1.194159    2.488474    1.608633 

# Final predictors for EI: AAO, FSLE, SSH, mixed layer depth, sea ice concentration, 
# daily difference in sea ice concentration, 16m salinity, 644m EKE
EI_predictors <- predictors



# KING GEORGE ISLAND
# Resetting model formula with initial predictors
mod_formula <- paste(species, "~", paste(intl_predictors, collapse = " + "))
all_pred <- lm(as.formula(mod_formula), data = KGI_binned)
vif(all_pred)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 1.642640        2.001346        9.082796        5.480220        8.625662       11.630452        1.191885 
# salinity_0   temperature_0           EKE_0      julian_day     salinity_16  temperature_16          EKE_16 
# 4490.336814     3135.273419        3.496817        4.990624     4510.308495     3241.691148        3.830063 
# salinity_644 temperature_644         EKE_644 
# 4.694717        6.604567        1.677522 

# Dropping 16m salinity
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "salinity_0", "temperature_0", "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[2]), paste0("temperature_", depths[3]),
                paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
KGI_vif <- lm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 1.628329        1.945211        9.082784        4.235576        8.236587       10.133419        1.179299 
# salinity_0   temperature_0           EKE_0      julian_day          EKE_16    salinity_644  temperature_16 
# 18.027630     3091.159881        3.496563        4.963000        3.824639        4.624503     3198.939699 
# temperature_644         EKE_644 
# 6.573822        1.664748  

# Dropping 16m temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "salinity_0", "temperature_0", "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]),
                paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
KGI_vif <- lm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 1.606782        1.765666        8.670186        4.198423        8.205081       10.021353        1.178679 
# salinity_0   temperature_0           EKE_0      julian_day          EKE_16    salinity_644 temperature_644 
# 17.224202       11.769338        3.274702        4.431687        3.761038        4.377699        6.562551 
# EKE_644 
# 1.662724 

# Dropping surface salinity
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "temperature_0", "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]),
                paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
KGI_vif <- lm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 1.538808        1.740245        8.609337        3.275024        6.703994        9.676998        1.174046 
# temperature_0           EKE_0      julian_day          EKE_16    salinity_644 temperature_644         EKE_644 
# 9.847539        3.220655        3.617189        3.596415        4.117196        5.245144        1.593458 

# Dropping surface temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]),
                paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
KGI_vif <- lm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 1.431294        1.626187        7.880023        2.918126        5.775639        6.447710        1.155809 
# EKE_0      julian_day          EKE_16    salinity_644 temperature_644         EKE_644 
# 3.206158        2.392978        3.570128        3.914101        3.894990        1.592896 

# Dropping ice thickness
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',
                "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]),
                paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
KGI_vif <- lm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc        ice_diff           EKE_0 
# 1.385979        1.603957        7.871014        2.725034        2.969906        1.149827        3.185822 
# julian_day          EKE_16    salinity_644 temperature_644         EKE_644 
# 2.387199        3.532840        3.832753        3.894548        1.577152 

# Final predictors for KGI: AAO, FSLE, SSH, mixed layer, sea ice concentration, ice concentration difference,
# EKE at all depths, julian day, 644 m salinity, and 644m temperature
KGI_predictors <- predictors



# CLARENCE ISLAND
# Resetting model formula with initial predictors
mod_formula <- paste(species, "~", paste(intl_predictors, collapse = " + "))
all_pred <- lm(as.formula(mod_formula), data = CI_binned)
vif(all_pred)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 2.288971        3.645411        3.738286        5.590789        6.296037        7.511069        1.579476 
# salinity_0   temperature_0           EKE_0      julian_day     salinity_16  temperature_16          EKE_16 
# 6192.248868     2993.056379        8.465314       28.208526     6067.205181     3103.553896        9.172470 
# salinity_644 temperature_644         EKE_644 
# 4.546499       24.016870        2.209480 

# Dropping surface salinity
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
               "temperature_0", "EKE_0", 'julian_day',
                paste0("salinity_", depths[2]), paste0("temperature_", depths[2]), paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 2.070446        3.586361        3.738115        4.511984        6.212274        7.505484        1.579226 
# temperature_0           EKE_0      julian_day     salinity_16  temperature_16          EKE_16    salinity_644 
# 2984.440938        8.455724       27.146048       14.008281     3095.418240        9.160126        4.233609 
# temperature_644         EKE_644 
# 22.710537        2.077174 

# Dropping 16m temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "temperature_0", "EKE_0", 'julian_day',
                paste0("salinity_", depths[2]), paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 2.050438        3.550654        3.735546        4.261158        5.812253        7.175590        1.575740 
# temperature_0           EKE_0      julian_day     salinity_16          EKE_16    salinity_644 temperature_644 
# 20.949234        8.441581       23.545641       13.615765        9.159036        4.203531       21.618133 
# EKE_644 
# 2.042239 

# Dropping julian day
# Not dropping surface temperature because it could be relevant based on timeseries
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "temperature_0", "EKE_0",
                paste0("salinity_", depths[2]), paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output: 
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 1.985927        3.542925        3.734861        4.020154        4.912431        6.915978        1.501892 
# temperature_0           EKE_0     salinity_16          EKE_16    salinity_644 temperature_644         EKE_644 
# 10.387105        8.441365       13.595355        9.124783        4.203372        6.649620        2.042183 

# Dropping 16m salinity
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "temperature_0", "EKE_0",
                paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output: 
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 1.981833        2.993650        3.281177        2.524620        4.749710        6.825823        1.498428 
# temperature_0           EKE_0          EKE_16    salinity_644 temperature_644         EKE_644 
# 4.878233        8.413480        8.881340        3.486637        5.993844        1.924609

# Dropping 16m EKE
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "temperature_0", "EKE_0",
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output: 
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 1.882573        2.983573        3.261736        2.518913        4.704156        6.823117        1.477403 
# temperature_0           EKE_0    salinity_644 temperature_644         EKE_644 
# 4.671390        1.598681        3.431271        5.993541        1.722165 

# Dropping ice thickness
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',
                "temperature_0", "EKE_0",
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output: 
# AAO            FSLE             SSH     mixed_layer        ice_conc        ice_diff   temperature_0 
# 1.868202        2.477801        3.221849        2.481793        1.501734        1.467765        4.592898 
# EKE_0    salinity_644 temperature_644         EKE_644 
# 1.538548        3.369883        5.181950        1.720147 

# Dropping 644m temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',
                "temperature_0", "EKE_0",
                paste0("salinity_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output: 
# AAO          FSLE           SSH   mixed_layer      ice_conc      ice_diff temperature_0         EKE_0 
# 1.575418      2.452815      2.198003      2.408869      1.501375      1.466250      4.432262      1.510546 
# salinity_644       EKE_644 
# 1.823191      1.697593 

# Final predictors for CI: AAO, FSLE, SSH, mixed layer, ice concentration, ice concentration difference,
# surface temperature/EKE, 644m salinity/EKE
CI_predictors <- predictors

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
# initial gam based on list of EI predictors
# going to deal with later because less data
EI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(SSH) + s(mixed_layer) + s(ice_conc) + 
                s(ice_diff) + s(salinity_16) + s(EKE_644),
              data = EI_binned, family=binomial, method='REML')

# KING GEORGE ISLAND
# initial gam based on list of KGI predictors
KGI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(ice_diff) + s(EKE_0) + s(julian_day) + 
                 s(EKE_16) + s(salinity_644) + s(temperature_644) + s(EKE_644), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 104.1551
# ................................................
# summary()
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(ice_diff) + 
#   s(EKE_0) + s(julian_day) + s(EKE_16) + s(salinity_644) + 
#   s(temperature_644) + s(EKE_644)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -6.553      2.429  -2.698  0.00698 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(AAO)             5.332  6.305 12.518 0.05509 . 
# s(FSLE)            1.000  1.000  3.089 0.07882 . 
# s(mixed_layer)     1.000  1.000  2.407 0.12077   
# s(ice_conc)        1.719  2.124  2.969 0.18604   
# s(ice_diff)        1.000  1.000  0.487 0.48509   
# s(EKE_0)           1.000  1.000  7.652 0.00567 **
#   s(julian_day)      1.000  1.000  1.619 0.20320   
# s(EKE_16)          1.000  1.000  7.404 0.00651 **
#   s(salinity_644)    2.706  3.326  2.663 0.41289   
# s(temperature_644) 1.000  1.000  1.829 0.17621   
# s(EKE_644)         1.627  2.045  3.781 0.15174   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.598   Deviance explained = 60.3%
# -REML = 38.712  Scale est. = 1         n = 181
# ..................................................
# gam.check()
# Method: REML   Optimizer: outer newton
# full convergence after 12 iterations.
# Gradient range [-8.971908e-06,1.917744e-06]
# (score 38.71168 & scale 1).
# Hessian positive definite, eigenvalue range [1.053381e-07,0.4808878].
# Model rank =  100 / 100 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                      k'  edf k-index p-value  
# s(AAO)             9.00 5.33    1.01    0.57  
# s(FSLE)            9.00 1.00    1.11    0.92  
# s(mixed_layer)     9.00 1.00    1.04    0.76  
# s(ice_conc)        9.00 1.72    1.10    0.91  
# s(ice_diff)        9.00 1.00    0.96    0.29  
# s(EKE_0)           9.00 1.00    1.04    0.76  
# s(julian_day)      9.00 1.00    0.91    0.14  
# s(EKE_16)          9.00 1.00    1.18    1.00  
# s(salinity_644)    9.00 2.71    0.97    0.38  
# s(temperature_644) 9.00 1.00    1.04    0.67  
# s(EKE_644)         9.00 1.63    0.88    0.04 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Removing ice difference
KGI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(EKE_0) + s(julian_day) +
                 s(EKE_16) + s(temperature_644) + s(EKE_644) + s(salinity_644), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 100.9424
# summary()
# ........................................
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(EKE_0) + 
#   s(julian_day) + s(EKE_16) + s(temperature_644) + s(EKE_644) + 
#   s(salinity_644)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -5.393      1.330  -4.055 5.01e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value   
# s(AAO)             5.368  6.340 13.774 0.04340 * 
#   s(FSLE)            1.000  1.000  3.064 0.08007 . 
# s(mixed_layer)     1.000  1.000  2.993 0.08361 . 
# s(ice_conc)        1.000  1.000  6.206 0.01274 * 
#   s(EKE_0)           1.000  1.000  7.749 0.00537 **
#   s(julian_day)      1.000  1.000  1.591 0.20714   
# s(EKE_16)          1.000  1.000  7.072 0.00783 **
#   s(temperature_644) 1.000  1.000  3.202 0.07355 . 
# s(EKE_644)         1.248  1.455  4.779 0.09613 . 
# s(salinity_644)    2.472  3.071  2.303 0.47957   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.585   Deviance explained = 58.5%
# -REML = 39.035  Scale est. = 1         n = 181
# ..........................................................
# gam.check()
# Method: REML   Optimizer: outer newton
# full convergence after 11 iterations.
# Gradient range [-3.711999e-06,9.265953e-07]
# (score 39.03456 & scale 1).
# Hessian positive definite, eigenvalue range [1.134204e-07,0.5342764].
# Model rank =  91 / 91 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                      k'  edf k-index p-value  
# s(AAO)             9.00 5.37    1.02    0.66  
# s(FSLE)            9.00 1.00    1.10    0.94  
# s(mixed_layer)     9.00 1.00    1.03    0.72  
# s(ice_conc)        9.00 1.00    1.09    0.92  
# s(EKE_0)           9.00 1.00    1.06    0.86  
# s(julian_day)      9.00 1.00    0.90    0.10  
# s(EKE_16)          9.00 1.00    1.17    0.99  
# s(temperature_644) 9.00 1.00    1.05    0.79  
# s(EKE_644)         9.00 1.25    0.87    0.04 *
#   s(salinity_644)    9.00 2.47    0.99    0.41  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Remove 644m salinity
KGI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(EKE_0) + s(julian_day) +
                 s(EKE_16) + s(temperature_644) + s(EKE_644), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 99.52605
# ..............................
# summary()
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(EKE_0) + 
#   s(julian_day) + s(EKE_16) + s(temperature_644) + s(EKE_644)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -4.589      1.011  -4.539 5.65e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(AAO)             5.527  6.522 12.838   0.0564 .  
# s(FSLE)            1.000  1.000  3.908   0.0481 *  
#   s(mixed_layer)     1.000  1.000  1.997   0.1577    
# s(ice_conc)        1.000  1.000  5.724   0.0167 *  
#   s(EKE_0)           1.000  1.000  9.377   0.0022 ** 
#   s(julian_day)      1.000  1.000  0.504   0.4776    
# s(EKE_16)          1.171  1.318  8.496   0.0105 *  
#   s(temperature_644) 1.000  1.000 15.195 9.66e-05 ***
#   s(EKE_644)         1.000  1.000  3.401   0.0652 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.548   Deviance explained = 55.5%
# -REML = 40.634  Scale est. = 1         n = 181
# ...............................
# gam.check()
# Method: REML   Optimizer: outer newton
# full convergence after 13 iterations.
# Gradient range [-9.027089e-06,3.597092e-06]
# (score 40.63419 & scale 1).
# Hessian positive definite, eigenvalue range [4.485224e-07,0.4910171].
# Model rank =  82 / 82 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                      k'  edf k-index p-value  
# s(AAO)             9.00 5.53    1.02    0.62  
# s(FSLE)            9.00 1.00    1.08    0.88  
# s(mixed_layer)     9.00 1.00    1.03    0.70  
# s(ice_conc)        9.00 1.00    1.11    0.95  
# s(EKE_0)           9.00 1.00    1.02    0.59  
# s(julian_day)      9.00 1.00    0.87    0.03 *
#   s(EKE_16)          9.00 1.17    1.18    0.99  
# s(temperature_644) 9.00 1.00    1.04    0.71  
# s(EKE_644)         9.00 1.00    0.87    0.04 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Remove julian day
KGI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(EKE_0) +
                 s(EKE_16) + s(temperature_644) + s(EKE_644), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 98.39254
# ................................
# summary()
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(EKE_0) + 
#   s(EKE_16) + s(temperature_644) + s(EKE_644)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.4737     0.9679  -4.622  3.8e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(AAO)             5.531  6.519 13.792  0.04086 *  
#   s(FSLE)            1.000  1.000  3.921  0.04769 *  
#   s(mixed_layer)     1.000  1.000  1.532  0.21589    
# s(ice_conc)        1.000  1.000  5.773  0.01627 *  
#   s(EKE_0)           1.000  1.000  9.164  0.00247 ** 
#   s(EKE_16)          1.435  1.741  8.861  0.01903 *  
#   s(temperature_644) 1.000  1.000 17.020 3.77e-05 ***
#   s(EKE_644)         1.000  1.000  3.140  0.07639 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.55   Deviance explained = 55.5%
# -REML = 41.188  Scale est. = 1         n = 181
# ....................................
# gam.check()
# Method: REML   Optimizer: outer newton
# full convergence after 11 iterations.
# Gradient range [-6.924208e-06,4.711926e-06]
# (score 41.18787 & scale 1).
# Hessian positive definite, eigenvalue range [2.599424e-07,0.5249838].
# Model rank =  73 / 73 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                      k'  edf k-index p-value  
# s(AAO)             9.00 5.53    1.02    0.62  
# s(FSLE)            9.00 1.00    1.08    0.90  
# s(mixed_layer)     9.00 1.00    1.03    0.68  
# s(ice_conc)        9.00 1.00    1.11    0.94  
# s(EKE_0)           9.00 1.00    1.01    0.61  
# s(EKE_16)          9.00 1.44    1.18    0.98  
# s(temperature_644) 9.00 1.00    1.04    0.73  
# s(EKE_644)         9.00 1.00    0.86    0.03 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Remove mixed layer depth
KGI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(ice_conc) + s(EKE_0) +
                 s(EKE_16) + s(temperature_644) + s(EKE_644), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 100.3554
# AIC increased, add mixed layer depth back in 
# Mess with number of knots (start with edf+2 and round up)
KGI_gam <- gam(Gm ~ s(AAO,k=8) + s(FSLE,k=3) + s(ice_conc,k=3) + s(EKE_0,k=3) + s(mixed_layer,k=3) +
                 s(EKE_16,k=4) + s(temperature_644,k=3) + s(EKE_644,k=3), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 98.53178
# ................................
# summary()
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(AAO, k = 8) + s(FSLE, k = 3) + s(ice_conc, k = 3) + s(EKE_0, 
#                                                                k = 3) + s(mixed_layer, k = 3) + s(EKE_16, k = 4) + s(temperature_644, 
#                                                                                                                      k = 3) + s(EKE_644, k = 3)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.3546     0.9231  -4.717 2.39e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(AAO)             4.892  5.618 14.101 0.02970 *  
#   s(FSLE)            1.000  1.000  3.837 0.05012 .  
# s(ice_conc)        1.000  1.000  5.767 0.01633 *  
#   s(EKE_0)           1.000  1.000  9.466 0.00209 ** 
#   s(mixed_layer)     1.000  1.000  1.486 0.22289    
# s(EKE_16)          1.583  1.925  9.012 0.01768 *  
#   s(temperature_644) 1.000  1.000 16.761 4.3e-05 ***
#   s(EKE_644)         1.000  1.000  3.150 0.07595 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.541   Deviance explained = 54.5%
# -REML = 41.455  Scale est. = 1         n = 181
# ....................................
# gam.check()
# Method: REML   Optimizer: outer newton
# full convergence after 11 iterations.
# Gradient range [-6.65886e-06,1.441089e-06]
# (score 41.45507 & scale 1).
# Hessian positive definite, eigenvalue range [7.908781e-08,0.366068].
# Model rank =  23 / 23 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                      k'  edf k-index p-value  
# s(AAO)             7.00 4.89    1.01    0.59  
# s(FSLE)            2.00 1.00    1.07    0.92  
# s(ice_conc)        2.00 1.00    1.10    0.94  
# s(EKE_0)           2.00 1.00    1.02    0.54  
# s(mixed_layer)     2.00 1.00    1.04    0.71  
# s(EKE_16)          3.00 1.58    1.19    0.99  
# s(temperature_644) 2.00 1.00    1.05    0.78  
# s(EKE_644)         2.00 1.00    0.86    0.04 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Removing mixed layer
KGI_gam <- gam(Gm ~ s(AAO,k=8) + s(FSLE,k=3) + s(ice_conc,k=3) + s(EKE_0,k=3) +
                 s(EKE_16,k=4) + s(temperature_644,k=3) + s(EKE_644,k=3), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 98.49182
# ............................
# summary()
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(AAO, k = 8) + s(FSLE, k = 3) + s(ice_conc, k = 3) + s(EKE_0, 
#                                                                k = 3) + s(EKE_16, k = 4) + s(temperature_644, k = 3) + s(EKE_644, 
#                                                                                                                          k = 3)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.0265     0.8117  -4.961 7.02e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(AAO)             4.770  5.519 13.255  0.03727 *  
#   s(FSLE)            1.000  1.000  2.553  0.11006    
# s(ice_conc)        1.000  1.000  5.649  0.01746 *  
#   s(EKE_0)           1.000  1.000  8.916  0.00283 ** 
#   s(EKE_16)          1.628  1.983  8.125  0.02175 *  
#   s(temperature_644) 1.000  1.000 26.391 1.01e-06 ***
#   s(EKE_644)         1.000  1.000  2.936  0.08663 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.525   Deviance explained = 53.1%
# -REML = 42.887  Scale est. = 1         n = 181
# ......................................
# gam.check()
# Method: REML   Optimizer: outer newton
# full convergence after 10 iterations.
# Gradient range [-7.737902e-06,1.22519e-06]
# (score 42.88697 & scale 1).
# Hessian positive definite, eigenvalue range [8.83043e-08,0.3466696].
# Model rank =  21 / 21 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                      k'  edf k-index p-value  
# s(AAO)             7.00 4.77    1.02   0.580  
# s(FSLE)            2.00 1.00    1.06   0.865  
# s(ice_conc)        2.00 1.00    1.09   0.910  
# s(EKE_0)           2.00 1.00    1.03   0.715  
# s(EKE_16)          3.00 1.63    1.18   0.995  
# s(temperature_644) 2.00 1.00    1.02   0.660  
# s(EKE_644)         2.00 1.00    0.85   0.045 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# CLARENCE ISLAND
# initial gam based on list of CI predictors
CI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(SSH) + s(mixed_layer) + s(ice_conc) + s(temperature_0) +
             s(EKE_0) + s(salinity_644) + s(EKE_644), 
             data = CI_binned, family = binomial,method='REML')
