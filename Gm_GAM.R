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

# Dropping mid-depth temperature and salinity
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "salinity_0", "temperature_0", "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 3.029628        1.924686        8.228406       12.702930       22.309522       32.073832        5.095357 
# salinity_0   temperature_0           EKE_0      julian_day          EKE_16    salinity_644 temperature_644 
# 25.308226      134.952089       40.976986      145.399103       55.595832       15.465427       15.897050 
# EKE_644 
# 18.196025  

# Dropping julian_day and surface temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "salinity_0", "EKE_0", paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 2.822872        1.757347        7.240826        7.329875       20.169989       22.988999        3.845513 
# salinity_0           EKE_0          EKE_16    salinity_644 temperature_644         EKE_644 
# 7.960467       34.457136       50.464982        8.165072        9.220861       17.521791 

# Dropping 16m EKE and ice thickness
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',
                "salinity_0", "EKE_0",
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc        ice_diff      salinity_0 
# 2.710477        1.602465        6.528602        6.773933        4.630897        1.531146        5.253967 
# EKE_0    salinity_644 temperature_644         EKE_644 
# 6.199916        6.957586        8.466284       10.464930 

# Dropping 644m temperature, salinity, and EKE
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',
                "salinity_0", "EKE_0")
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
EI_vif <- lm(as.formula(mod_formula), data = EI_binned)
vif(EI_vif)
# output:
# AAO        FSLE         SSH mixed_layer    ice_conc    ice_diff  salinity_0       EKE_0 
# 1.804394    1.581153    2.245877    2.197384    3.110817    1.187833    2.406984    1.279214

# Final predictors for EI: AAO, FSLE, SSH, mixed layer depth, sea ice concentration, 
# daily difference in sea ice concentration, surface salinity, surface EKE
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

# Dropping 16m salinity and temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "salinity_0", "temperature_0", "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
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

# Dropping surface salinity, surface temp, and ice thickness
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", 'ice_diff',
                "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
KGI_vif <- lm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc        ice_diff           EKE_0 
# 1.385979        1.603957        7.871014        2.725034        2.969906        1.149827        3.185822 
# julian_day          EKE_16    salinity_644 temperature_644         EKE_644 
# 2.387199        3.532840        3.832753        3.894548        1.577152  

# Dropping SSH
predictors <-c("AAO", "FSLE", "mixed_layer", "ice_conc", 'ice_diff',
               "EKE_0", 'julian_day', paste0("EKE_", depths[2]),
               paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
KGI_vif <- lm(as.formula(mod_formula), data = KGI_binned)
vif(KGI_vif)
# output:
# AAO            FSLE     mixed_layer        ice_conc        ice_diff           EKE_0      julian_day 
# 1.281701        1.593272        2.703774        2.832040        1.140709        3.179582        1.552111 
# EKE_16    salinity_644 temperature_644         EKE_644 
# 3.527536        3.142213        1.891018        1.535272

# Final predictors for KGI: AAO, FSLE, mixed  layer, sea ice concentration, ice concentration, difference,
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

# Dropping 16m salinity and temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "salinity_0", "temperature_0", "EKE_0", 'julian_day',
                paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output:
# AAO            FSLE             SSH     mixed_layer        ice_conc   ice_thickness        ice_diff 
# 2.053265        3.533095        3.735524        4.381670        5.803185        7.173994        1.575813 
# salinity_0   temperature_0           EKE_0      julian_day          EKE_16    salinity_644 temperature_644 
# 13.908656       20.969301        8.439893       23.530914        9.153159        4.246715       21.519816 
# EKE_644 
# 2.054408  

# Dropping julian day and 644m temperature
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "salinity_0", "temperature_0", "EKE_0",
                paste0("EKE_", depths[2]),
                paste0("salinity_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output:
# AAO          FSLE           SSH   mixed_layer      ice_conc ice_thickness      ice_diff    salinity_0 
# 1.692651      3.415489      3.162713      3.814860      4.654989      6.242105      1.501851     12.578456 
# temperature_0         EKE_0        EKE_16  salinity_644       EKE_644 
# 10.469495      8.438828      9.103590      3.400284      2.052730 

# Dropping 16m EKE and surface salinity
# Not dropping surface temperature because it could be relevant based on timeseries
predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                "temperature_0", "EKE_0",
                paste0("salinity_", depths[3]), paste0("EKE_", depths[3]))
mod_formula <- paste(species, "~", paste(predictors, collapse = " + "))
CI_vif <- lm(as.formula(mod_formula), data = CI_binned)
vif(CI_vif)
# output: 
# AAO          FSLE           SSH   mixed_layer      ice_conc ice_thickness      ice_diff temperature_0 
# 1.583100      2.970412      2.232839      2.486441      4.293383      5.899192      1.477403      4.445036 
# EKE_0  salinity_644       EKE_644 
# 1.594409      1.874840      1.697775 

# Dropping ice thickness
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


# ELEPHANT ISLAND
# initial gam based on list of EI predictors
# going to deal with later because less data
EI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(SSH) + s(mixed_layer) + s(ice_conc) + s(salinity_0) + s(EKE_0),
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

# Removing nonsigniicant terms ice difference and 644m salinity
KGI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(EKE_0) + s(julian_day) +
                 s(EKE_16) + s(temperature_644) + s(EKE_644), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 99.52605
# summary()
# ........................................
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
# ..........................................................
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
# s(AAO)             9.00 5.53    1.02   0.685  
# s(FSLE)            9.00 1.00    1.08   0.880  
# s(mixed_layer)     9.00 1.00    1.03   0.700  
# s(ice_conc)        9.00 1.00    1.11   0.960  
# s(EKE_0)           9.00 1.00    1.02   0.620  
# s(julian_day)      9.00 1.00    0.87   0.050 *
#   s(EKE_16)          9.00 1.17    1.18   1.000  
# s(temperature_644) 9.00 1.00    1.04   0.805  
# s(EKE_644)         9.00 1.00    0.87   0.055 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Removing julian day because least significant
KGI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(mixed_layer) + s(ice_conc) + s(EKE_0) +
                 s(EKE_16) + s(temperature_644) + s(EKE_644), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 98.39254
# ..........................
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
# .................................
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
# s(AAO)             9.00 5.53    1.02   0.650  
# s(FSLE)            9.00 1.00    1.08   0.900  
# s(mixed_layer)     9.00 1.00    1.03   0.670  
# s(ice_conc)        9.00 1.00    1.11   0.915  
# s(EKE_0)           9.00 1.00    1.01   0.605  
# s(EKE_16)          9.00 1.44    1.18   1.000  
# s(temperature_644) 9.00 1.00    1.04   0.705  
# s(EKE_644)         9.00 1.00    0.86   0.035 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Removing mixed layer depth
KGI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(ice_conc) + s(EKE_0) +
                 s(EKE_16) + s(temperature_644) + s(EKE_644), 
               data = KGI_binned, family=binomial,method='REML')
# AIC()
# 100.3554
# .................................
# summary()
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   Gm ~ s(AAO) + s(FSLE) + s(ice_conc) + s(EKE_0) + s(EKE_16) + 
#   s(temperature_644) + s(EKE_644)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -4.718      1.687  -2.796  0.00517 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(AAO)             5.456  6.463 13.332  0.04668 *  
#   s(FSLE)            1.000  1.000  2.833  0.09233 .  
# s(ice_conc)        1.745  2.197  2.337  0.27728    
# s(EKE_0)           1.000  1.000  8.410  0.00373 ** 
#   s(EKE_16)          1.444  1.754  8.373  0.02338 *  
#   s(temperature_644) 1.320  1.568 20.922 3.25e-05 ***
#   s(EKE_644)         1.000  1.000  2.600  0.10685    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.545   Deviance explained = 55.2%
# -REML =  42.62  Scale est. = 1         n = 181
# .................................
# gam.check()
# Method: REML   Optimizer: outer newton
# full convergence after 24 iterations.
# Gradient range [-5.100595e-06,8.130405e-07]
# (score 42.62003 & scale 1).
# Hessian positive definite, eigenvalue range [2.071747e-06,0.4878189].
# Model rank =  64 / 64 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                      k'  edf k-index p-value  
# s(AAO)             9.00 5.46    1.01    0.62  
# s(FSLE)            9.00 1.00    1.06    0.85  
# s(ice_conc)        9.00 1.75    1.11    0.98  
# s(EKE_0)           9.00 1.00    1.02    0.67  
# s(EKE_16)          9.00 1.44    1.18    0.99  
# s(temperature_644) 9.00 1.32    1.01    0.59  
# s(EKE_644)         9.00 1.00    0.85    0.03 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# CLARENCE ISLAND
# initial gam based on list of CI predictors
CI_gam <- gam(Gm ~ s(AAO) + s(FSLE) + s(SSH) + s(mixed_layer) + s(ice_conc) + s(temperature_0) +
             s(EKE_0) + s(salinity_644) + s(EKE_644), 
             data = CI_binned, family = binomial,method='REML')
