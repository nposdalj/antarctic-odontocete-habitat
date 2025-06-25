library(tidyverse)
library(ggfortify)
library(splines)

# Loading data and correctly formatting date/time 
odontocete <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/Antarc_Odontocetes.csv")
odontocete$Start.time <- as.POSIXct(odontocete$Start.time, format ="%Y-%m-%d %H:%M:%S")
odontocete$End.time <- as.POSIXct(odontocete$End.time, format = "%Y-%m-%d %H:%M:%S")
# Filling in cells that were not correctly reformatted
odontocete[1174,4] <- mdy_hms("2/28/2015 0:00:00") 
odontocete[1374,4] <- mdy_hms("1/15/2016 0:00:00")

# Function to create a minute-by-minute species-specific dataframe
presence_df <- function(data, site) {
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
  
  # Creating a dataframe with the minutes and Julian day for the chosen site
  minutes <- data.frame(Minute = seq(from = t1, to = t2, by = "1 min"))
  minutes$Julian.Day <- yday(minutes$Minute)
  
  # Expand each call into several minute-by-minute presences
  expanded <- filtered %>%
    rowwise() %>%
    mutate(Minute = list(seq(floor_date(Start.time, "minute"), 
                             floor_date(End.time, "minute"), by = "1 min"))) %>%
    unnest(cols = c(Minute)) %>%
    select(Species.Code, Minute) %>%
    mutate(Presence = 1)
  
  # Takes the dataframe from the previous step and splits it up by species columns 
  presence <- expanded %>%
    distinct() %>%
    pivot_wider(names_from = Species.Code, values_from = Presence, values_fill = 0)
  
  # Combine with the empty minutes dataframe and fill minutes without whale presence as 0
  minbymin <- minutes %>%
    left_join(presence, by = "Minute") %>%
    replace(is.na(.), 0)  # fill missing with 0s (absence)
  
  return(minbymin)
}

# Creating minute-by-minute dataframe for each site
EI_presence <- presence_df(odontocete, "EI")
KGI_presence <- presence_df(odontocete, "KGI")
CI_presence <- presence_df(odontocete, "CI")

site <- c("EI", "KGI", "CI")
BW29 <- c(0, 0, 0)
BW37 <- c(0, 0, 0)
BW58 <- c(0, 0, 0)
Oo <- c(0, 0, 0)
Gm <- c(0, 0, 0)
Pm <- c(0, 0, 0)
species_list <- list("BW29", "BW37", "BW58", "Oo", "Gm", "Pm")
acf_table <- data.frame(site, BW29, BW37, BW58, Oo, Gm, Pm)

FindCorrelation <- function(site_presence) {
  for (x in 1:6) {
    species <- species_list[[x]]
    formula <- as.formula(paste(species, "~ bs(Julian.Day, k=4)"))
    BlockMod <- glm(formula, data = site_presence, family = binomial)
    ACF = acf(residuals(BlockMod),lag.max = 5000)
    CI = ggfortify:::confint.acf(ACF)
    ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
    ACFval = ACFidx[1]
    print(paste("Species:", species, "ACFval:", ACFval))
  }
}
FindCorrelation(EI_presence)
# [1] "Species: BW29 ACFval: 3899"
# [1] "Species: BW37 ACFval: 21"
# [1] "Species: BW58 ACFval: 2412"
# [1] "Species: Oo ACFval: 52"
# [1] "Species: Gm ACFval: NA"
# [1] "Species: Pm ACFval: 4498"
# Warning messages:
#   1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 3: glm.fit: fitted probabilities numerically 0 or 1 occurred 



# BlockMod<-glm(PreAbs~
#                 bs(Julian,k=4)+
#                 TimeLost+
#                 as.factor(Year)
#               ,data=SiteHourTable,family=binomial)
# 
# ACF = acf(residuals(BlockMod),lag.max = 100)
# CI = ggfortify:::confint.acf(ACF)
# ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
# if (site == 'NFC'){
#   ACFval = 81
# } else 
#   if (site == 'HAT_B'){
#     ACFval = 100
#   } else { 
#     ACFval = ACFidx[1] }