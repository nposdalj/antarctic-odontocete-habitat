# Script to add all environmental data to daily presence dataframe for final modeling
# Can also plot any custom site/species/variable combinations
library(tidyverse)
library(gridExtra)
library(patchwork)

# ---------------- Step 0: Create base final dataframe-------------
dailyDetection <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/dailyDetections.csv")
allData <- dailyDetection
allData$date <- allData$Day
allData <- allData %>% subset(select = -Day)
# Choose variables + sites + species to plot
environmental_vars <- c('AAO','FSLE','SSH','EKE','temperature','salinity', 'chla','o2','productivity',
                        'mixed_layer','ice_conc','ice_thickness','ice_diff') 
# options: AAO, FSLE, SSH, EKE, temperature, salinity, mixed_layer, chla (chlorophyll), o2 (oxygen),
#   productivity (primary production), ice_conc (sea ice concentration), ice_thickness (sea ice thickness),
#   ice_diff (daily difference in ice concentration)
species <- c('Gm') # options: BW29, BW37, BW58, Oo, Pm, Gm
# BW29 = Southern bottlenose whale, BW37 & BW58 = Gray's and strap-toothed whales
# Oo = Killer whale, Pm = Sperm Whale, Gm = Long-finned pilot whale
sites <- c('KGI', 'EI','CI') # options: EI, KGI, CI

# --------------- Step 1: Format/Add Antarctic Oscillation Index ----------------
AAO <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Daily_AAO.csv")

# Make date column
# Add zero in front of single digit dates
addZero <- function(number) {
  if(number == '1') { return('01')
  } else if(number == '2') { return('02')
  } else if(number == '3') { return('03')
  } else if(number == '4') { return('04')
  } else if(number == '5') { return('05')
  } else if(number == '6') { return('06')
  } else if(number == '7') { return('07')
  } else if(number == '8') { return('08')
  } else if(number == '9') { return('09')}
}
# Parse through relevant years to add zero
AAO <- filter(AAO, year >= 2014 & year <= 2016)
AAO$month <- as.character(AAO$month)
AAO$day <- as.character(AAO$day)
for(x in 1:1096) {
  if(AAO$month[x] == '1' | AAO$month[x] == '2' | AAO$month[x] == '3' | AAO$month[x] == '4' | AAO$month[x] == '5' | 
     AAO$month[x] == '6' | AAO$month[x] == '7' | AAO$month[x] == '8' | AAO$month[x] == '9') {
    AAO$month[x] <- addZero(AAO$month[x])
  }
  if(AAO$day[x] == '1' | AAO$day[x] == '2' | AAO$day[x] == '3' | AAO$day[x] == '4' | AAO$day[x] == '5' | 
     AAO$day[x] == '6' | AAO$day[x] == '7' | AAO$day[x] == '8' | AAO$day[x] == '9') {
    AAO$day[x] <- addZero(AAO$day[x])
  }
}
# Create date column with year, month, day
AAO$date <- paste(AAO$year, AAO$month, AAO$day, sep = "")
AAO$date <- as.Date(AAO$date, "%Y%m%d")
# Keep date column & AAO index, only for relevant date ranges
AAO <- AAO %>% subset(select = c(aao_index_cdas, date)) %>%
  filter((date >= '2014-03-05' & date <= '2014-07-17') | 
                        (date >= '2015-02-10' & date <= '2016-01-29') | 
                        (date >= '2016-02-04' & date <= '2016-12-02'))
# Add AAO to final dataset
allData <- merge(allData, AAO, by=intersect(names(allData), names(AAO)))
allData <- rename(allData, AAO = aao_index_cdas)


# -------------------- Step 2: Add Copernicus Data--------
# sst, salinity, depth variables, eke, ssh, etc.
# ice variables, chlorophyll, oxygen
copernicus <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/copernicus_40km.csv")
allData <- merge(allData, copernicus, by=intersect(names(allData), names(copernicus)))
allData <- allData %>% rename(SSH=ssh) %>% rename(temperature=temp) %>%
  rename(ice_conc=sice_conc) %>% rename(ice_thickness=sice_thick)
# replacing NA in ice variables with 0
allData[is.na(allData)] <- 0
# Adding ice difference column
allData <- allData %>% mutate(ice_diff = ice_conc - lag(ice_conc, n = 1))
allData <- allData %>% mutate(ice_diff = if_else(date %in% as.Date(
  c('2014-03-05', '2015-02-10', '2016-02-04')), 0, ice_diff))

# -------------------- Step 3: Format/Add FSLEs------------
EI_fsle <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/AVISO/EI_fsle_40km")
KGI_fsle <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/AVISO/KGI_fsle_40km")
CI_fsle <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/AVISO/CI_fsle_40km")

EI_fsle$Site <- "EI"
KGI_fsle$Site <- "KGI"
CI_fsle$Site <- "CI"

all_fsle <- rbind(EI_fsle, KGI_fsle, CI_fsle)
all_fsle$X <- NULL

allData <- merge(allData, all_fsle, by=intersect(names(allData), names(all_fsle)))
allData <- rename(allData, FSLE=fsle)

# -------------------- Step 4: Save Final Dataframe -------------------
# Pivoting data to wide format, so there is one column for temperature, salinity, EKE at each depth
# Removing unnecessary columns (might need north/east sea ice velocity later for convergence calculations?)
allData <- allData %>% subset(select=-c(X,n_velocity,e_velocity,sice_e_veloc,sice_n_veloc))

# Setting variable categories for current dataframe
depth_varying <- c("temperature", "salinity", 'EKE','chla','productivity','o2')
surf_vars <- c("AAO",'SSH','mixed_layer','ice_conc','ice_thickness','FSLE','ice_diff')
species_vars <- c('BW29','BW37','BW58','Gm',"Pm", "Oo")
grouping_vars <- c("date", "depth", "Site")

# Pivoting depth-varying variables to a wide format (one column for each depth)
depth_wide <- allData %>% select(all_of(grouping_vars), all_of(depth_varying)) %>%
  pivot_longer(cols = all_of(depth_varying), names_to = "variable", values_to = "value") %>%
  mutate(variable = paste0(variable, "_", round(depth))) %>%
  select(-depth) %>% pivot_wider(names_from = variable, values_from = value)  

# Making separate dataframes for non depth-varying and species presence variables
surf_wide <- allData %>% select(date, Site, all_of(surf_vars)) %>% distinct()
species_wide <- allData %>% select(date, Site, all_of(species_vars)) %>% distinct()

# Joining all dataframes into wide format and saving file
allData_wide <- depth_wide %>% left_join(surf_wide, by = c("date", "Site")) %>%
  left_join(species_wide, by = c("date", "Site")) %>% distinct()
allData_wide$julian_day <- yday(allData_wide$date)

write.csv(allData_wide, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/allData.csv")



# -------------------- Step 5: Make Requested Plots -------------------
timeseriesPlots <- function(sites, species, vars) {
  # Keeping only requested variables, species, and sites
  filtered <- allData %>% select(any_of(c('date','Site','depth',species, vars))) %>%
    subset(Site %in% sites)
  filtered$date <- as.Date(filtered$date,'%Y-%m-%d')

  # Run through all sites
  for(si in sites) {
    
    for(sp in species){ # Make plots for each species at this site
      # Setting and filtering by species dive depths
      if(sp =='BW29') {
        depths <- c(0.5, 768.0)}
      else if(sp=='BW37' | sp=='BW58') {
        depths <- c(0.5, 67.0, 920.0) }
      else if(sp=='Oo') {
        depths <- c(0.5, 11.0, 455.0) }
      else if(sp=='Pm') {
        depths <- c(0.5, 375.0, 1665.0)}
      else if(sp=='Gm') {
        depths <- c(0.5, 16.0, 635.0) }
      else {
        print("Species code not valid. Check inputs.") }
      species_specific <- filtered %>% select(any_of(c('date','Site','depth',sp, vars))) %>%
        subset(depth %in% depths) %>% filter(Site == si)
      
      # Making final figure
      aggregatePlot(species_specific, vars, sp, depths,si)
      print(paste(si, " plot for ",sp, " done."))
    }
  }
}


# Function to combine all plots for one site/species combination into one figure
aggregatePlot <- function(data, vars, plot_species, depths,site) {

  # Creating list of plots
  all_plots <- list()

  # Making plots and adding to list of plots
  for(v in vars) {
    v_idx <- match(v,vars)
    plot <- makePlot(data, v, depths)
    all_plots[[length(all_plots)+1]] <- plot
  }
  
  # Adding timeseries for species presence to list of plots
  species_ts <- ggplot(data = data, mapping = aes(x = date, y = get(plot_species))) + 
    geom_col(width = 1, color = "purple") + scale_x_date(date_labels = "%b %Y") +
    labs(y = NULL, x = paste(plot_species,' Presence')) +
    theme(plot.margin = unit(c(0, 0.5, 0, 0.5),units = "line"), 
          axis.title = element_text(size=10))
  all_plots[[length(all_plots)+1]] <- species_ts
  
  # Arranging the final figure
  final_plot <- wrap_plots(all_plots, ncol = 1, guides = "collect") &
    plot_annotation(title = paste(name(plot_species)," at ", name(site))) &
    # Making margins small to fit plots with many variables
    theme(legend.position = "bottom", legend.margin = margin(t=0, b=2, unit="pt"),
          legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.spacing = unit(1, "pt"), legend.key.size = unit(8, "pt"),
          legend.title = element_text(size=10))
  # Showing plot
  print(final_plot)
  return(final_plot)
}

# Function to make an individual timeseries plot for one environmental variable
makePlot <- function(data, var, depths) {
  # Setting units for chosen variable
  if(var=='SSH'){
    label <- 'SSH (m)'
    color <- 'darkred'
  } else if(var=='mixed_layer'){
    label <- 'Mixed Layer Depth (m)'
    color <- 'darkcyan'
  } else if(var=='ice_thickness'){
    label <- 'Sea Ice Thickness (m)'
    color <- 'navy'
  } else if(var=='temperature') {
    label <- 'Temperature (C)'
  } else if(var=='salinity') {
    label <- 'Salinity (psu)'
  } else if(var=='ice_conc') {
    label <- 'Sea Ice Concentration'
    color <- 'mediumslateblue'
  } else if(var=='EKE') {
    label <- "EKE (cm^2/s^2)"
  } else if(var=='AAO') {
    label <- 'AAO'
    color <- 'mediumvioletred'
  } else if(var == 'FSLE') {
      label <- 'FSLE'
      color <- 'orange'
  } else if(var == 'o2') {
    label <- 'Oxygen Concentration (mmol/m3)'
    #color <- 'saddlebrown'
  } else if(var == 'chla') {
    label <- 'Chlorophyll (mg/m3)'
    #color <- 'darkolivegreen'
  } else if(var == 'productivity') {
    label <- "Net Primary Production (mg/m3/day carbon)"
    #color <- 'dimgray'
  } else if(var == 'ice_diff') {
    label <- "Daily Change in Sea Ice concentration"
    color <- 'saddlebrown'
  }
  
  # Plotting across depths for relevant variables
  if (var=='temperature' | var == 'EKE' | var == 'salinity' | var=='o2' | var=='productivity' | 
      var == 'chla') {
    # Renaming the surface depth
    data <- data %>% mutate(depth = ifelse(depth == 0.5, 'surface', as.character(depth)))
    # Returning plot colored by depth
    return(ggplot(data = data, aes(x=date, y=get(var), color=factor(depth))) + 
             geom_line() + 
             labs(y = NULL, color = "Depth (m)", x = label) +  scale_x_date(labels = NULL) +
             theme(plot.margin = unit(c(0, 0.5, 0, 0.5),units = "line"), 
                   axis.title = element_text(size=10)))
  } else { # Plotting surface values for variables without depth
    return(
      ggplot(data = data, mapping = aes(x = date, y = get(var))) + geom_line(color = color) +
        xlab(label) + ylab(NULL) +  scale_x_date(labels = NULL) +
        theme(plot.margin = unit(c(0, 0.5, 0, 0.5),units = "line"), 
              axis.title = element_text(size=10))
    )
  }
}

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

timeseriesPlots(sites, species, environmental_vars)
