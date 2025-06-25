library(tidyverse)
library(ggfortify) # For calculating acf
library(splines) # For creating model to calculate acf
library(photobiology) # For extracting day/night information

# Loading data and correctly formatting date/time 
odontocete <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/Antarc_Odontocetes.csv")
odontocete$Start.time <- as.POSIXct(odontocete$Start.time, 
                                    format ="%Y-%m-%d %H:%M:%S", tz = "GMT")
odontocete$End.time <- as.POSIXct(odontocete$End.time, 
                                  format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
# Filling in cells that were not correctly reformatted
odontocete[1174,4] <- mdy_hms("2/28/2015 0:00:00") 
odontocete[1374,4] <- mdy_hms("1/15/2016 0:00:00")
latlong <- tibble::tibble(lon = -56, lat = -61.5)

# ---------- Step 1: Making the presence dataframe ----------
# Function to create a minute-by-minute species-specific dataframe
HourByHour_df <- function(data, site) {
  # Filtering out the other call types
  filtered <- data %>% 
    filter(Site == site, Call %in% c("Clicks", "regular-click"))
  
  # Setting date range bounds for each site
  if (site == "EI") {
    t1 <- mdy_hms("03-05-2014 00:00:00")
    t2 <- mdy_hms("07-17-2014 12:05:05")
  } else if (site == "KGI") {
    t1 <- mdy_hms("02-10-2015 00:00:00")
    t2 <- mdy_hms("01-29-2016 13:16:40")
  } else if (site == "CI") {
    t1 <- mdy_hms("02-04-2016 0:00:00")
    t2 <- mdy_hms("12-02-2016 5:16:20")
  } 
  
  # Creating a dataframe with the hours and Julian day for the chosen site
  hours <- data.frame(Hour = seq(from = t1, to = t2, by = "60 min"))
  hours$Julian.Day <- yday(hours$Hour)
  
  # Expand each call into several hour-by-hour presences
  expanded <- filtered %>%
    rowwise() %>%
    mutate(Hour = list(seq(floor_date(Start.time, "60 minutes"), 
                             floor_date(End.time, "60 minutes"), by = "60 min"))) %>%
    unnest(cols = c(Hour)) %>%
    select(Species.Code, Hour) %>%
    mutate(Presence = 1)
  
  # Takes the dataframe from the previous step and splits it up by species columns 
  presence <- expanded %>%
    distinct() %>%
    pivot_wider(names_from = Species.Code, values_from = Presence, values_fill = 0)
  
  # Combine with the empty hours dataframe and fill hours without whale presence as 0
  hourbyhour <- hours %>%
    left_join(presence, by = "Hour") %>%
    replace(is.na(.), 0)  # fill missing with 0s (absence)

  # Add day/night column 
  hourbyhour$Diel <- as.factor(ifelse(hourbyhour$Hour > 
                            sunrise_time(date = date(hourbyhour$Hour), 
                                         tz = "GMT", geocode = latlong)
                          & hourbyhour$Hour <
                            sunset_time(date = date(hourbyhour$Hour),
                                        tz = "GMT", geocode = latlong),
                          "day", "night"))

  return(hourbyhour)
}
# Creating 5 minute-by-minute dataframe for each site
EI_presence <- HourByHour_df(odontocete, "EI")
KGI_presence <- HourByHour_df(odontocete, "KGI")
CI_presence <- HourByHour_df(odontocete, "CI")
# Adding columns for species with no detections
CI_presence$BW58 <- 0
KGI_presence$Pm <- 0



# ---------- Step 2: Generating tables -----------
# Creating variables to build three tables:
# 1) Total species detection hours by site 
# 2) ACF values by site and species
# 3) PACF values by site and species
site <- c("EI", "KGI", "CI")
BW29 <- c(0, 0, 0)
BW37 <- c(0, 0, 0)
BW58 <- c(0, 0, 0)
Oo <- c(0, 0, 0)
Gm <- c(0, 0, 0)
Pm <- c(0, 0, 0)
species_list <- list("BW29", "BW37", "BW58", "Oo", "Gm", "Pm")
detection_table <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)
acf_table <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)
pacf_table <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)


# ---------- Step 2a: Detection Table -----------
# Function to build detection table
SiteDetection <- function(site, site_presence, detection_table) {
  row_index <- which(detection_table$site == site) # Row index for the site
  for (x in 1:6) { # Fill in table
    species <- species_list[[x]]
    species_data <- site_presence[[species]]
    detection_table[row_index, species] <- sum(species_data)
  }
  return(detection_table)
}
# Fill in table for each site
detection_table <- SiteDetection("EI", EI_presence, detection_table)
detection_table <- SiteDetection("KGI", KGI_presence, detection_table)
detection_table <- SiteDetection("CI", CI_presence, detection_table)

# ------------ Step 2b: ACF table --------------
# Function to build acf table
FindACF <- function(site, site_presence, acf_table) {
  row_index <- which(acf_table$site == site) # Row index for the site
  # Running loop for each species
  for (x in 1:6) {
    # Creating relevant variables
    species <- species_list[[x]]
    species_data <- site_presence[[species]]
    species_sum <- sum(species_data)
    
    # Base code to generate acf plots/values can be found in Natalie's GitHub:
    # SeasonalityAnalysis\GAM_GEEs\WAT\General\SiteSpecific_GAMGEE_knots4_WAT.R
    
    if (species_sum == 0) { # If no data for this species, return 0 for acf_value
      acf_table[row_index, species] <- 0 }
    else { 
      formula <- as.formula(paste(species, "~ bs(Julian.Day, k=4) + Diel"))
      BlockMod <- glm(formula, data = site_presence, family = binomial)
      ACF = acf(residuals(BlockMod), lag.max = 1500)
      CI = ggfortify:::confint.acf(ACF)
      ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
      ACFval = ACFidx[1]
      acf_table[row_index, species] <- ACFval
    }
  }
  return(acf_table)
}
# Filling in the acf table with data from each site
acf_table <- FindACF("EI", EI_presence, acf_table)
acf_table <- FindACF("KGI", KGI_presence, acf_table)
acf_table <- FindACF("CI", CI_presence, acf_table)

# ------------ Step 2c: PACF table --------------
# Function to build pacf table
FindPACF <- function(site, site_presence, pacf_table) {
  row_index <- which(pacf_table$site == site) # Row index for the site
  # Running loop for each species
  for (x in 1:6) {
    # Creating relevant variables
    species <- species_list[[x]]
    species_data <- site_presence[[species]]
    species_sum <- sum(species_data)
    
    # Base code to generate acf plots/values can be found in Natalie's GitHub:
    # SeasonalityAnalysis\GAM_GEEs\WAT\General\SiteSpecific_GAMGEE_knots4_WAT.R
    
    if (species_sum == 0) { # If no data for this species, return 0 for acf_value
      pacf_table[row_index, species] <- 0 }
    else { 
      formula <- as.formula(paste(species, "~ bs(Julian.Day, k=4) + Diel"))
      BlockMod <- glm(formula, data = site_presence, family = binomial)
      PACF = pacf(residuals(BlockMod), lag.max = 100)
      CI = ggfortify:::confint.acf(PACF)
      PACFidx = which(PACF[["acf"]] < CI, arr.ind=TRUE)
      PACFval = PACFidx[1]
      pacf_table[row_index, species] <- PACFval
    }
  }
  return(pacf_table)
}
# Filling in the pacf table with data from each site
pacf_table <- FindPACF("EI", EI_presence, pacf_table)
pacf_table <- FindPACF("KGI", KGI_presence, pacf_table)
pacf_table <- FindPACF("CI", CI_presence, pacf_table)
