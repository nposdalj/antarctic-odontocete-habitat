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
allData <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/allData_40km_FSLEhfdeg.csv")
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


# ------------- Step 2: Average by ACF EDIT STEP ------------
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
  
  # EDIT HERE
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
if(species != 'BW29') {
  intl_predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                       "salinity_0", "temperature_0", "EKE_0",'o2_0','chla_0','productivity_0', 'julian_day',
                       paste0("salinity_", depths[2]), paste0("temperature_", depths[2]), paste0("EKE_", depths[2]),
                       paste0("salinity_", depths[3]), paste0("temperature_", depths[3]), paste0("EKE_", depths[3]),
                       paste0('chla_',depths[2]), paste0('o2_',depths[2]), paste0('productivity_',depths[2]),
                       paste0('chla_',depths[3]), paste0('o2_',depths[3]), paste0('productivity_',depths[3]))
} else
  intl_predictors <- c("AAO", "FSLE", "SSH", "mixed_layer", "ice_conc", "ice_thickness", 'ice_diff',
                       "salinity_0", "temperature_0", "EKE_0", 'julian_day',
                       paste0("salinity_", depths[2]), paste0("temperature_", depths[2]), paste0("EKE_", depths[2]))
mod_formula <- paste(species, "~", paste(intl_predictors, collapse = " + "))
