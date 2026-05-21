
## Read in Daily Detections .csv file
dailyDetections <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/dailyDetections.csv") # replace this folder path with where you download the dailyDetections.csv
dailyDetections$Day <- as.Date(as.character(dailyDetections$Day), format = "%m/%d/%Y")


# CLARENCE ISLAND - CI
# Pm
species <- "Pm" # change to "BW29", "Gm", etc.
site <- "CI" # change to CI, KGI, EI

## Subset the data for one site
site_data <- dailyDetections[dailyDetections$Site == site, ]

# Generate plot
plot(site_data$Day, site_data[[species]], type = "h",
     xlab = "Date", ylab = "Presence (0/1)",
     main = paste("Presence of", species, "at", site))

# Oo
species <- "Oo" # change to "BW29", "Gm", etc.
site <- "CI" # change to CI, KGI, EI

## Subset the data for one site
site_data <- dailyDetections[dailyDetections$Site == site, ]

# Generate plot
plot(site_data$Day, site_data[[species]], type = "h",
     xlab = "Date", ylab = "Presence (0/1)",
     main = paste("Presence of", species, "at", site))

