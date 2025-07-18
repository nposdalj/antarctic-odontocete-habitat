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






# ---------- Step 1a: Making hourly presence dataframe ----------
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
  hours$Day <- yday(hours$Hour)
  
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
# Creating hour dataframe for each site
EI_hr <- HourByHour_df(odontocete, "EI")
KGI_hr <- HourByHour_df(odontocete, "KGI")
CI_hr <- HourByHour_df(odontocete, "CI")
# Adding columns for species with no detections
CI_hr$BW58 <- 0
KGI_hr$Pm <- 0


# ---------- Step 1b: Making daily presence dataframe ----------
# Function to create a minute-by-minute species-specific dataframe
DayByDay_df <- function(data, site) {
  # Filtering out the other call types
  filtered <- data %>% 
    filter(Site == site, Call %in% c("Clicks", "regular-click"))
  
  # Setting date range bounds for each site
  if (site == "EI") {
    t1 <- as.Date("03/05/2014",format = "%m/%d/%Y")
    t2 <- as.Date("07/17/2014",format = "%m/%d/%Y")
  } else if (site == "KGI") {
    t1 <- as.Date("02/10/2015",format = "%m/%d/%Y")
    t2 <- as.Date("01/29/2016", format = "%m/%d/%Y")
  } else if (site == "CI") {
    t1 <- as.Date("02/04/2016", format = "%m/%d/%Y")
    t2 <- as.Date("12/02/2016", format = "%m/%d/%Y")
  } 
  
  # Creating a dataframe with the days for the chosen site
  days <- data.frame(Day = seq.Date(from = t1, to = t2, by = 1))

  # Expand each call into several daily presences
  expanded <- filtered %>%
    rowwise() %>%
    mutate(Day = list(seq(floor_date(Start.time, "1 day"), 
                           floor_date(End.time, "1 day"), by = "1 day"))) %>%
    unnest(cols = c(Day)) %>%
    select(Species.Code, Day) %>%
    mutate(Presence = 1)
  
  # Takes the dataframe from the previous step and splits it up by species columns 
  presence <- expanded %>%
    distinct() %>%
    pivot_wider(names_from = Species.Code, values_from = Presence, values_fill = 0)
  
  # Combine with the empty hours dataframe and fill hours without whale presence as 0
  daybyday <- days %>%
    left_join(presence, by = "Day") %>%
    replace(is.na(.), 0)  # fill missing with 0s (absence)
  
  return(daybyday)
}
# Creating hour dataframe for each site
EI_day <- DayByDay_df(odontocete, "EI")
KGI_day <- DayByDay_df(odontocete, "KGI")
CI_day <- DayByDay_df(odontocete, "CI")
# Adding columns for species with no detections
CI_day$BW58 <- 0
KGI_day$Pm <- 0

# ----------------Step 1c: Exporting daily data-----------
# Exporting daily data so a timeseries can be made from it (TimeseriesPlots.R)
CI_day$Site <- "CI"
KGI_day$Site <- "KGI"
EI_day$Site <- "EI"
dailyDetections <- rbind(CI_day, KGI_day, EI_day)
write.csv(dailyDetections, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/dailyDetections.csv",
          row.names = F)



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
hrDetection <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)
dayDetection <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)
hrResidACF <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)
hrModACF <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)
dayResidACF <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)
dayModACF <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)
#pacf_table <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)


# ---------- Step 2a: Detection Tables (hours and days) -----------
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
# Fill in hourly table for each site
hrDetection <- SiteDetection("EI", EI_hr, hrDetection)
hrDetection <- SiteDetection("KGI", KGI_hr, hrDetection)
hrDetection <- SiteDetection("CI", CI_hr, hrDetection)

# Fill in daily table for each site
dayDetection <- SiteDetection("EI", EI_day, dayDetection)
dayDetection <- SiteDetection("KGI", KGI_day, dayDetection)
dayDetection <- SiteDetection("CI", CI_day, dayDetection)

# ------------ Step 2b: ACF table (residuals) --------------
# Function to build acf table on residuals
ResidACF <- function(site, site_presence, acf_table) {
  row_index <- which(acf_table$site == site) # Row index for the site
  # Running loop for each species
  for (x in 1:6) {
    # Creating relevant variables
    species <- species_list[[x]]
    species_data <- site_presence[[species]]
    species_sum <- sum(species_data)
    
    # Base code to generate acf plots/values can be found in Natalie's GitHub:
    # SeasonalityAnalysis/GAM_GEEs/WAT/General/SiteSpecific_GAMGEE_knots4_WAT.R
    
    if (species_sum == 0) { # If no data for this species, return 0 for acf_value
      acf_table[row_index, species] <- 0 }
    else { 
      formula <- as.formula(paste(species, "~ bs(Day, k=4)")) # ADD/REMOVE + Diel based on hour or day
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
# Make sure to include + Diel in function formula before running
# Filling in the acf hourly residual table with data from each site
hrResidACF <- ResidACF("EI", EI_hr, hrResidACF)
hrResidACF <- ResidACF("KGI", KGI_hr, hrResidACF)
hrResidACF <- ResidACF("CI", CI_hr, hrResidACF)

# Make sure to remove + Diel in function formula before running
# Filling in the acf daily residual table with data from each site
dayResidACF <- ResidACF("EI", EI_day, dayResidACF)
dayResidACF <- ResidACF("KGI", KGI_day, dayResidACF)
dayResidACF <- ResidACF("CI", CI_day, dayResidACF)


# ------------ Step 2c: ACF table (model data) --------------
# Function to build acf table on model data
ModACF <- function(site, site_presence, acf_table) {
  row_index <- which(acf_table$site == site) # Row index for the site
  # Running loop for each species
  for (x in 1:6) {
    # Creating relevant variables
    species <- species_list[[x]]
    species_data <- site_presence[[species]]
    species_sum <- sum(species_data)
    
    # Base code to generate acf plots/values can be found in Natalie's GitHub:
    # SeasonalityAnalysis/GAM_GEEs/WAT/General/SiteSpecific_GAMGEE_knots4_WAT.R
    
    if (species_sum == 0) { # If no data for this species, return 0 for acf_value
      acf_table[row_index, species] <- 0 }
    else { 
      formula <- as.formula(paste(species, "~ bs(Day, k=4)")) # ADD/REMOVE + Diel based on hour or day
      BlockMod <- glm(formula, data = site_presence, family = binomial)
      ACF = acf(fitted(BlockMod), lag.max = 1500) 
      CI = ggfortify:::confint.acf(ACF)
      ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
      ACFval = ACFidx[1]
      acf_table[row_index, species] <- ACFval
    }
  }
  return(acf_table)
}
# Make sure to include + Diel in function formula before running
# Filling in the acf hourly residual table with data from each site
hrModACF <- ModACF("EI", EI_hr, hrModACF)
hrModACF <- ModACF("KGI", KGI_hr, hrModACF)
hrModACF <- ModACF("CI", CI_hr, hrModACF)

# Make sure to remove + Diel in function formula before running
# Filling in the acf daily residual table with data from each site
dayModACF <- ModACF("EI", EI_day, dayModACF)
dayModACF <- ModACF("KGI", KGI_day, dayModACF)
dayModACF <- ModACF("CI", CI_day, dayModACF)

# Saving acf data
write.csv(dayResidACF, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Autocorrelation/acf_table.csv", row.names = F)
write.csv(dayDetection, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Autocorrelation/binned_detections.csv", row.names = F)

# # ------------ NOT USING: Step 2c: PACF table --------------
# # Function to build pacf table
# FindPACF <- function(site, site_presence, pacf_table) {
#   row_index <- which(pacf_table$site == site) # Row index for the site
#   # Running loop for each species
#   for (x in 1:6) {
#     # Creating relevant variables
#     species <- species_list[[x]]
#     species_data <- site_presence[[species]]
#     species_sum <- sum(species_data)
#     
#     # Base code to generate acf plots/values can be found in Natalie's GitHub:
#     # SeasonalityAnalysis/GAM_GEEs/WAT/General/SiteSpecific_GAMGEE_knots4_WAT.R
#     
#     if (species_sum == 0) { # If no data for this species, return 0 for acf_value
#       pacf_table[row_index, species] <- 0 }
#     else { 
#       formula <- as.formula(paste(species, "~ bs(Day, k=4) + Diel"))
#       BlockMod <- glm(formula, data = site_presence, family = binomial)
#       PACF = pacf(residuals(BlockMod), lag.max = 100)
#       CI = ggfortify:::confint.acf(PACF)
#       PACFidx = which(PACF[["acf"]] < CI, arr.ind=TRUE)
#       PACFval = PACFidx[1]
#       pacf_table[row_index, species] <- PACFval
#     }
#   }
#   return(pacf_table)
# }
# # Filling in the pacf table with data from each site
# pacf_table <- FindPACF("EI", EI_presence, pacf_table)
# pacf_table <- FindPACF("KGI", KGI_presence, pacf_table)
# pacf_table <- FindPACF("CI", CI_presence, pacf_table)
# # ----------------Step 3: Exporting tables---------------
# directory <- "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Autocorrelation"
# write.csv(acf_table, paste(directory, "/acf_table.csv", sep = ""), row.names = F)
# write.csv(pacf_table, paste(directory, "/pacf_table.csv", sep = ""), row.names = F)
# write.csv(detection_table, paste(directory, "/binned_detections.csv", sep = ""), row.names = F)
