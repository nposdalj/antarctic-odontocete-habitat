library(tidyverse)
library(gridExtra)

# Load data from EIE/CI site
CI <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Antarc_EIE_01_Odontocetes.csv")
CI <- CI[-c(1:2,7)]
# Load data from SSI/KGI site
KGI <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Antarc_SSI_01_Odontocetes.csv")
KGI <- KGI[-c(1:2,7)]
# Load data from EI site
EI <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Antarc_EI_01_Odontocetes.csv")
EI <- EI[-c(1:2,7:12)]
# Join all data into one dataframe and save it
CI$Site <- "CI"
KGI$Site <- "KGI"
EI$Site <- "EI"
odontocete <- rbind(CI, KGI, EI)
write.csv(odontocete,"C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Antarc_Odontocetes.csv")

# Timeseries plots
odontocete <- odontocete %>% # Reformat dataframe
  mutate(
    Start.time = mdy_hms(Start.time), # Change from character to time format
    End.time = mdy_hms(End.time),
    Date = as.Date(Start.time), # New column that saves the date of call
    Call.time = as.numeric(difftime(End.time, Start.time, units = "secs"))
    # New column that has length of call in seconds
  )
odontocete$Call.time[odontocete$Call.time == 0] <- 1 # If a call is 0 seconds, correct to 1

timeseries <- function(site, species) { # Function to create a timeseries plot
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
  daily <- filtered %>% group_by(Date) %>% 
    summarise(Duration = sum(Call.time, na.rm = TRUE))
  daily$Duration <- daily$Duration / 3600
  # Creating the site- and species-specific timeseries
  ggplot(data = daily, mapping = aes(x = Date, y = Duration)) + geom_col() + 
    scale_x_date(limits = c(b1, b2), date_labels = "%b %Y")+ 
    labs(subtitle = species, y = "Duration (hrs)")
}

bySite_ts <- function(site) { # Function to aggregate all plots from a particular site
  # Joining all species-specific timeseries into one figure
  final_plot <- grid.arrange(timeseries(site, "BW29"), timeseries(site, "BW37"), 
                             timeseries(site, "BW58"), timeseries(site, "Gm"), 
                             timeseries(site, "Oo"), timeseries(site, "Pm"),
                             nrow = 6, top = site)
  return(final_plot)
}

# Generating timeseries plots for Elephant Island, King George Island, and Clarence Island
EI_ts <- bySite_ts("EI")
KGI_ts <- bySite_ts("KGI")
CI_ts <- bySite_ts("CI")
