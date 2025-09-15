library(tidyverse)
library(gridExtra)
library(lubridate)
library(purrr)

# -----------------Step 1: Prepping data---------------------
# Load data from EIE/CI site
CI <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/Antarc_EIE_01_Odontocetes.csv")
CI <- CI[-c(1:2,7)]
# Load data from SSI/KGI site
KGI <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/Antarc_SSI_01_Odontocetes.csv")
KGI <- KGI[-c(1:2,7)]
# Load data from EI site
EI <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/Antarc_EI_01_Odontocetes.csv")
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
write.csv(odontocete,"/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/Antarc_Odontocetes.csv")



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
  #windows(width = 20, height = 30)
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


# ----------------Step 3: Daily Binned Timeseries (By Site) ------------
dailyDetections <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/dailyDetections.csv")

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
  #windows(width = 20, height = 30)
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

# ----------------Step 3: Daily Binned Timeseries (By Species) ------------
combined_ts <- function(species) { # Function to create & aggregate all plots from a particular species
  # Prepping data
  dailyDetections$Day <- as.Date(dailyDetections$Day, "%Y-%m-%d")
  # Adding a week start date variable to group data
  dailyDetections$week_start <- floor_date(dailyDetections$Day, unit = "week")
  
  # Creating dataframe that stores the total call days for each week, only for the relevant species
  week <- dailyDetections %>% group_by(week_start) %>% 
    reframe(total_days = sum(.data[[species]],na.rm=TRUE),Site = Site,week_start = week_start) %>% unique()
  
  
  # Dataframe of weeks with partial recording effort
  # Use in a geom_point layer that indicates weeks with partial effort
  # Weeks with partial effort:
  # 2014-03-02: starting 3/05 (4 days on)
  # 2014-07-13: ending 7/17 (5 days on)
  # 2015-02-08: starting 2/10 (5 days on)
  # 2016-01-31: ending 1/29 (5 days on)
  partialWeeks <- c('2014-03-02','2014-07-13','2015-02-08','2016-01-31')
  partialWeeks <- as.Date(partialWeeks, format = '%Y-%m-%d')
  daysRecorded <- c(4.2, 5.2, 5.2, 5.2) # .1 to nudge point away from columns
  partialEffort <- data.frame(partialWeeks,daysRecorded)
  
  # timeseries across all site
  ggplot(data = week, mapping = aes(x = week_start, y = total_days)) + 
    # adding species presence
    geom_col(width = 3, mapping = aes(fill = Site)) +
    
    # adding weeks with partial recording effort
    geom_point(data=partialEffort, aes(x=partialWeeks,y=daysRecorded), color ='slateblue', alpha = 0.7) +
    
    # adding gray rectangle where no recording effort was made
    annotate("rect", xmin=as.Date('2014-7-20'), xmax=as.Date('2015-02-01'), 
             ymin=-Inf, ymax=Inf, alpha=0.5, fill="gray") +
    
    # scaling, theming, labeling plot
    scale_x_date(breaks = seq(as.Date('2014-03-05'), as.Date('2016-12-02'), by = '4 months'),
                 date_labels = "%b %Y") +
    labs(title = paste0(name(species),' Acoustic Presence'), y = 'Number of Days Detected', 
         x = 'Week') + 
    scale_y_continuous(limits=c(0, 7.5), expand = c(0, 0)) +
    scale_fill_manual(values = c("EI" = "deeppink", 
                                 "KGI" = "darkmagenta", "CI" = "royalblue")) +
    theme_bw() + theme(plot.subtitle = element_text(size = 9, face = "bold"),
                       plot.margin = unit(c(0.5, 2, 0.5, 1), units='line'),
                       legend.position = 'right', 
                       legend.margin = margin(t = 0, unit='cm'),
                       legend.title = element_text(size=9,face='bold'),
                       legend.text = element_text(size=8),
                       legend.key.size = unit(.8,'line'),
                       plot.title =element_text(face='bold'))
}

Gm_week <- combined_ts('Gm')
Oo_week <- combined_ts('Oo')
Pm_week <- combined_ts('Pm')
BW37_week <- combined_ts('BW37')
BW29_week <- combined_ts('BW29')
BW58_week <- combined_ts('BW58')

# ----------------- Step 4: 5-min binary presence (dailyDetections-like) -----------------
# Site bounds (use the same windows you used for plotting)
# Site bounds (use the same windows you used for plotting)
site_bounds <- function(site) {
  if (site == "EI") {
    b1 <- as.Date("2014-03-05"); b2 <- as.Date("2014-07-17")
  } else if (site == "KGI") {
    b1 <- as.Date("2015-02-10"); b2 <- as.Date("2016-01-29")
  } else if (site == "CI") {
    b1 <- as.Date("2016-02-04"); b2 <- as.Date("2016-12-02")
  } else {
    b1 <- min(odontocete$Date, na.rm = TRUE)
    b2 <- max(odontocete$Date, na.rm = TRUE)
  }
  # Cover entire days: start at 00:00 of b1, end at 23:55 of b2
  start_dt <- as.POSIXct(b1, tz = "UTC")
  end_dt   <- as.POSIXct(b2 + 1, tz = "UTC") - minutes(5)  # <-- changed
  list(start = start_dt, end = end_dt)
}


# Build full 5-min grid for a site
make_site_bins <- function(site) {
  rng <- site_bounds(site)
  tibble(
    Bin  = seq(from = rng$start, to = rng$end, by = "5 min"),
    Site = site
  )
}

# For a given site+species, return just the bins that are present (1’s)
present_bins_for <- function(site, species) {
  df <- odontocete %>%
    filter(Site == site, Species.Code == species) %>%
    transmute(Start = Start.time, End = End.time) %>%
    drop_na(Start, End)
  
  if (nrow(df) == 0) {
    return(tibble(Bin = as.POSIXct(character()), present = integer()))
  }
  
  # Snap each detection to 5-min bin sequence it touches
  df %>%
    mutate(
      start_bin = floor_date(Start, "5 minutes"),
      end_bin   = floor_date(End - seconds(1), "5 minutes")
    ) %>%
    rowwise() %>%
    mutate(Bin = list(seq(from = start_bin, to = end_bin, by = "5 min"))) %>%
    unnest(Bin) %>%
    distinct(Bin) %>%
    mutate(present = 1L) %>%
    select(Bin, present)
}

# Build the wide 5-min table that mirrors dailyDetections layout
build_5minDetections <- function(
    sites = c("EI","KGI","CI"),
    species_order = c("Pm","BW29","Oo","Gm","BW37","BW58")
) {
  # Start with a full bin grid per site
  all_bins <- map_df(sites, make_site_bins)
  
  # Add each species as a column of 0/1
  out <- all_bins
  for (sp in species_order) {
    sp_bins <- map_df(sites, ~present_bins_for(.x, sp) %>% mutate(Site = .x))
    out <- out %>%
      left_join(sp_bins, by = c("Bin", "Site")) %>%
      mutate(present = replace_na(present, 0L)) %>%
      rename(!!sp := present)
  }
  
  # Order & select columns to match your daily file (first col = time bin instead of Day)
  out %>%
    arrange(Site, Bin) %>%
    select(Bin, all_of(species_order), Site)
}

# ---- run & save ----
fiveMinDetections <- build_5minDetections()
# set the species columns in the order you use elsewhere
species_order <- c("Pm","BW29","Oo","Gm","BW37","BW58")

daily_minutes <- fiveMinDetections %>%
  mutate(Day = as.Date(Bin)) %>%
  group_by(Site, Day) %>%
  summarise(
    across(all_of(species_order), ~ sum(.x, na.rm = TRUE) * 5),  # 5 min per bin
    .groups = "drop"
  ) %>%
  select(Day, all_of(species_order), Site) %>%   # match your daily layout
  arrange(Site, Day)
write.csv(
  daily_minutes,
  "/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/Antarc_Odontocetes_daily_minutes.csv",
  row.names = FALSE
)
