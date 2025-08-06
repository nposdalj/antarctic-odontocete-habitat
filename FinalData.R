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
environmental_vars <- as.vector(c('AAO','FSLE','SSH','EKE','temperature','salinity', 'chla','o2',
                                  'productivity','mixed_layer','ice_conc','ice_thickness')) 
# options: AAO, FSLE, SSH, EKE, temperature, salinity, mixed_layer, chla (chlorophyll), o2 (oxygen),
#   productivity (primary production), ice_conc (sea ice concentration), ice_thickness (sea ice thickness),
#   ice_diff (daily difference in ice concentration), fsle_orient (orientation of fsle vector)
species <- as.vector(c('none')) # options: BW29, BW37, BW58, Oo, Pm, Gm, none (for surface variables)
# BW29 = Southern bottlenose whale, BW37 & BW58 = Gray's and strap-toothed whales
# Oo = Killer whale, Pm = Sperm Whale, Gm = Long-finned pilot whale
sites <- as.vector(c('KGI')) # options: EI, KGI, CI

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
allData <- allData %>% subset(select=-c(n_velocity_mean,e_velocity_mean,sice_e_veloc_mean,
                                        sice_n_veloc_mean,n_velocity_mad,e_velocity_mad))
allData <- allData %>% rename(SSH=ssh_mean, temperature=temp_mean, ice_conc=sice_conc_mean, 
                              ice_thickness=sice_thick_mean, salinity=salinity_mean, 
                              mixed_layer=mixed_layer_mean, chla=chla_mean, o2=o2_mean, productivity=productivity_mean,
                              EKE = EKE_mean)

# replacing NA in ice variables with 0
allData[is.na(allData)] <- 0

# Adding ice difference column
ice_diffData <- allData %>% select(date, ice_conc) %>% distinct() %>% arrange(date) %>%
  mutate(ice_diff = ice_conc - lag(ice_conc, n=1))
ice_diffData <- ice_diffData %>% mutate(ice_diff = if_else(date %in% as.Date(
  c('2014-03-05', '2015-02-10', '2016-02-04')), 0, ice_diff))
allData <- allData %>%
  left_join(ice_diffData %>% select(date, ice_diff), by = "date")
allData <- unique(allData)

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
allData <- rename(allData, FSLE=fsle_mean, fsle_orient = fsle_orientmean,
                  fsle_orient_sd = fsle_orientsd)

# -------------------- Step 4: Save Final Dataframe -------------------
# Pivoting data to wide format, so there is one column for temperature, salinity, EKE at each depth
# Removing unnecessary columns (might need north/east sea ice velocity later for convergence calculations?)
# adding ice regime and site depths

# only run this line if there is a column called X (remnant title on .csv import)
# allData <- allData %>% subset(select=-X)

# Setting variable categories for current dataframe
depth_varying <- c("temperature", "salinity", 'EKE','chla','productivity','o2',
                   'temp_sd','salinity_sd','EKE_mad','chla_sd','productivity_sd','o2_sd')
surf_vars <- c("AAO",'SSH','mixed_layer','ice_conc','ice_thickness','FSLE','fsle_orient','ice_diff',
               'ssh_sd','mixed_layer_sd','fsle_sd','fsle_orient_sd')
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

# adding in different ice regimes and site depth
allData_wide$ice_regime <- 'blank'
allData_wide$bathymetry <- 0
for (x in 1:nrow(allData_wide)) {
  if(allData_wide[x,'ice_diff'] == 0) {
    allData_wide[x,'ice_regime'] <- 'none'
  } else if(allData_wide[x,'ice_diff'] <= -0.01) {
    allData_wide[x,'ice_regime'] <- 'decreasing'
  } else if(allData_wide[x,'ice_diff'] >= 0.01) {
    allData_wide[x,'ice_regime'] <- 'increasing'
  } else
    allData_wide[x,'ice_regime'] <- 'stable'
  
  if(allData_wide[x,'Site'] == 'EI' | allData_wide[x,'Site'] == 'KGI'){
    allData_wide[x,'bathymetry'] <- 760
  } else
    allData_wide[x,'bathymetry'] <- 1030
}

write.csv(allData_wide, "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/allData.csv")



# -------------------- Step 5: Make Requested Plots -------------------
timeseriesPlots <- function(sites, species, vars) {
  # Keeping only requested variables, species, and sites
  if(species == 'none') {
    filtered <- allData %>% select(any_of(c('date','Site','depth', vars))) %>%
      subset(Site %in% sites)
  } else {
    filtered <- allData %>% select(any_of(c('date','Site','depth',species, vars))) %>%
      subset(Site %in% sites)
  }
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
      else if(sp=='none') {
        depths <- c(0.5) }
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
    plot <- makePlot(data, v, depths, plot_species)
    all_plots[[length(all_plots)+1]] <- plot
  }
  
  if(plot_species == 'none') 
    title <- paste0('Environmental Variables at ',name(site))
  else {
    # Adding timeseries for species presence to list of plots
    species_ts <- ggplot(data = data, mapping = aes(x = date, y = .data[[plot_species]])) + 
      geom_col(width = 1, color = "purple") + scale_x_date(date_labels = "%b %Y") +
      labs(y = NULL, x = paste(plot_species,' Presence')) +
      theme(plot.margin = unit(c(0, 0.5, 0, 0.5),units = "line"), 
            axis.title = element_text(size=10))  
    all_plots[[length(all_plots)+1]] <- species_ts
    
    title <- paste0(name(plot_species)," at ", name(site))
  }
  
  # Arranging the final figure
  final_plot <- wrap_plots(all_plots, ncol = 1, guides = "collect") &
    plot_annotation(title = title) &
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
makePlot <- function(data, var, depths, species) {
  # Setting labels for chosen variable
  if(var=='SSH'){
    label <- 'Sea Surface Height (m)'
    col <- 'darkred'
  } else if(var=='mixed_layer'){
    label <- 'Mixed Layer Depth (m)'
    col <- 'darkcyan'
  } else if(var=='ice_thickness'){
    label <- 'Sea Ice Thickness (m)'
    col <- 'navy'
  } else if(var=='temperature') {
    label <- 'Temperature (°C)'
  } else if(var=='salinity') {
    label <- 'Salinity (psu)'
  } else if(var=='ice_conc') {
    label <- 'Sea Ice Concentration'
    col <- 'mediumslateblue'
  } else if(var=='EKE') {
    label <- "Eddy Kinetic Energy (cm^2/s^2)"
  } else if(var=='AAO') {
    label <- 'Antarctic Oscillation Index'
    col <- 'mediumvioletred'
  } else if(var == 'FSLE') {
    label <- 'FSLE'
    col <- 'orange'
  } else if(var == 'o2') {
    label <- 'Oxygen Concentration (mmol/m3)'
  } else if(var == 'chla') {
    label <- 'Chlorophyll (mg/m3)'
  } else if(var == 'productivity') {
    label <- "Net Primary Production (mg/m3/day carbon)"
  } else if(var == 'ice_diff') {
    label <- "Daily Change in Sea Ice concentration"
    col <- 'saddlebrown'
  } else if(var == 'fsle_orient') {
    label <- "Orientation of FSLE Vector"
    col <- 'dimgray'
  }
  
  
  # Making surface plots for all variables if species is not specified
  if(species == 'none') {
    # setting colors for surface variables if not plotting across depths
    if(var == 'o2') {
      label <- 'Oxygen Concentration (mmol/m^3'
      col <- 'darkviolet'
    } else if(var == 'EKE') {
      label <- 'Eddy Kinetic Energy (cm^2/s^2)'
      col <- 'darkorange'
    } else if(var == 'temperature') {
      label <- 'Sea Surface Temperature (°C)'
      col <- 'maroon'
    } else if(var == 'salinity') {
      label <- 'Sea Surface Salinity (psu)'
      col <- 'darkgoldenrod'
    } else if(var == 'productivity') {
      label <- 'Net Primary Production (mg/m^3/day carbon)'
      col <- 'seagreen'
    } else if(var == 'chla') {
      label <- 'Chlorophyll-a (mg/m^3)'
      col <- 'darkolivegreen'
    }
    return(
      ggplot(data = data, mapping = aes(x = date, y = .data[[var]])) + geom_line(color = col) +
        xlab(label) + ylab(NULL) +  scale_x_date(labels = NULL) +
        theme(plot.margin = unit(c(0, 0.5, 0, 0.5),units = "line"), 
              axis.title = element_text(size=10))  
    )
  } else {
    # Plotting across depths for relevant variables
    # Only if species is specified
    if (var=='temperature' | var == 'EKE' | var == 'salinity' | var=='o2' | var=='productivity' | 
        var == 'chla') {
      # Renaming the surface depth
      data <- data %>% mutate(depth = ifelse(depth == 0.5, 'surface', as.character(depth)))
      # Returning plot colored by depth
      return(ggplot(data = data, aes(x=date, y=.data[[var]], color=factor(depth))) + 
               geom_line() + 
               labs(y = NULL, color = "Depth (m)", x = label) +  scale_x_date(labels = NULL) +
               theme(plot.margin = unit(c(0, 0.5, 0, 0.5),units = "line"), 
                     axis.title = element_text(size=10)))  
    } else { # Plotting surface values for variables without depth
      return(
        ggplot(data = data, mapping = aes(x = date, y = .data[[var]])) + geom_line(color = col) +
          xlab(label) + ylab(NULL) +  scale_x_date(labels = NULL) +
          theme(plot.margin = unit(c(0, 0.5, 0, 0.5),units = "line"), 
                axis.title = element_text(size=10)) 
      )
    }
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
