library(tidyverse)
library(gridExtra)
library(lubridate)

# -----------------Step 1: Prepping data---------------------
# Load data from EIE/CI site
CI <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/Antarc_EIE_01_Odontocetes.csv")
CI <- CI[-c(1:2,7)]
# Load data from SSI/KGI site
KGI <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/Antarc_SSI_01_Odontocetes.csv")
KGI <- KGI[-c(1:2,7)]
# Load data from EI site
EI <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/Antarc_EI_01_Odontocetes.csv")
EI <- EI[-c(1:2,7:12)]
# Join all data into one dataframe and save it
CI$Site <- "CI"
KGI$Site <- "KGI"
EI$Site <- "EI"
odontocete <- rbind(CI, KGI, EI)

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
odontocete <- odontocete %>% # Reformat dataframe
  mutate(
    Start.time = mdy_hms(Start.time), # Change from character to time format
    End.time = mdy_hms(End.time),
    Date = as.Date(Start.time), # New column that saves the date of call
    Call.time = as.numeric(difftime(End.time, Start.time, units = "secs"))
    # New column that has length of call in seconds
  )
odontocete$Call.time[odontocete$Call.time == 0] <- 1 # If a call is 0 seconds, correct to 1
write.csv(odontocete,"C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/Antarc_Odontocetes.csv")



# --------------Step 2: Cumulative daily hours timeseries---------------
hrTimeseries <- function(site, species) { # Function to create a timeseries plot
  # Setting date range bounds for each site
  if (site == "EI") {
    b1 <- as.Date("2014-03-05")
    b2 <- as.Date("2014-07-17")
  } else if (site == "KGI") {
    b1 <- as.Date("2015-02-10")
    b2 <- as.Date("2016-01-29")
  } else if (site == "CI") {
    b1 <- as.Date("2016-02-04")
    b2 <- as.Date("2016-12-02")
  } else {
    b1 <- min(odontocete$Date, na.rm = TRUE)
    b2 <- max(odontocete$Date, na.rm = TRUE)
  }
  # Filtering dataframe by the relevant site and species
  filtered <- filter(odontocete, Site == site & Species.Code == species)
  # Creating new dataframe that stores the total call duration (hours) for each day
  daily <- filtered %>% group_by(Date, Call) %>%
    summarise(Duration = sum(Call.time, na.rm = TRUE), .groups = "drop")
  daily$Duration <- daily$Duration / 3600
  # Creating the site- and species-specific timeseries
  ggplot(data = daily, mapping = aes(x = Date, y = Duration, fill = Call)) + geom_col(width = 0.9) + 
    scale_x_date(limits = c(b1, b2), date_labels = "%b %Y")+ 
    labs(subtitle = name(species), y = "Duration (hrs)") + 
    guides(fill = guide_legend(keywidth = .6, keyheight = 0.4)) +
    theme(
      legend.position = "bottom", legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 8), 
      plot.subtitle = element_text(size = 9, face = "bold"), axis.title = element_text(size = 9),
      axis.text = element_text(size = 7), legend.margin = margin(c(0, 0, 0.2, 0), unit = 'cm'), 
      plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), units = "line")
      )
}

hrSite_ts <- function(site) { # Function to aggregate all plots from a particular site
  # Joining all species-specific timeseries into one figure
  windows(width = 20, height = 30)
  final_plot <- grid.arrange(hrTimeseries(site, "BW29"), hrTimeseries(site, "BW37"), 
                             hrTimeseries(site, "BW58"), hrTimeseries(site, "Gm"), 
                             hrTimeseries(site, "Oo"), hrTimeseries(site, "Pm"),
                             nrow = 6, top = name(site))
  return(final_plot)
}

# Generating timeseries plots for Elephant Island, King George Island, and Clarence Island
EI_ts <- hrSite_ts("EI")
KGI_ts <- hrSite_ts("KGI")
CI_ts <- hrSite_ts("CI")


# ----------------Step 3: Daily Binned Timeseries------------
dailyDetections <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/dailyDetections.csv")

dayTimeseries <- function(site, species) { # Function to create a timeseries plot
  # Filtering dataframe by the relevant site
  filtered <- filter(dailyDetections, Site == site)
  filtered$Day <- as.Date(filtered$Day, "%Y-%m-%d")
  
  # Setting date range bounds
  b1 <- min(filtered$Day, na.rm = TRUE)
  b2 <- max(filtered$Day, na.rm = TRUE)
  
  # Adding a week start date variable to group data
  filtered$week_start <- floor_date(filtered$Day, unit = "week")
  
  # Creating new dataframe that stores the total call days for each week, only for the relevant species
  weekly <- filtered %>% group_by(week_start) %>%
    summarise(total_days = sum(.data[[species]], na.rm = TRUE), .groups = "drop")
  # Creating the site- and species-specific timeseries
  ggplot(data = weekly, mapping = aes(x = week_start, y = total_days)) + geom_col(width = 1, color = "mediumvioletred") +
    scale_x_date(limits = c(b1, b2), date_labels = "%b %Y")+
    labs(subtitle = name(species), y = NULL, x = NULL) + ylim(0,7) +
    theme(
      plot.subtitle = element_text(size = 9, face = "bold"), 
      plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), units = "line")
   )
}

daySite_ts <- function(site) { # Function to aggregate all plots from a particular site
  # Joining all species-specific timeseries into one figure
  windows(width = 20, height = 30)
  final_plot <- grid.arrange(dayTimeseries(site, "BW29"), dayTimeseries(site, "BW37"), 
                             dayTimeseries(site, "BW58"), dayTimeseries(site, "Gm"), 
                             dayTimeseries(site, "Oo"), dayTimeseries(site, "Pm"),
                             nrow = 6, top = name(site), left = "Number of Days Detected",
                             bottom = "Week")
  return(final_plot)
}

# Generating timeseries plots for Elephant Island, King George Island, and Clarence Island
EI_week <- daySite_ts("EI")
KGI_week <- daySite_ts("KGI")
CI_week <- daySite_ts("CI")