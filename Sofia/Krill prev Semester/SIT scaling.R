dailyDetections <- read.csv("F://dailyDetections.csv")
dailyDetections$Day <- as.Date(dailyDetections$Day, format="%m/%d/%Y")

species <- "Pm"    # detection column
site <- "CI"       # CI, EI, KGI

# Subset data for the site
df <- dailyDetections[dailyDetections$Site == site, ]

# Week-start date (previous Sunday)
df$WeekStart <- df$Day - as.numeric(format(df$Day, "%w"))

# --- Weekly detection counts (0–7) ---
weekly_presence <- tapply(df[[species]], df$WeekStart, sum, na.rm = TRUE)

# --- Weekly sea-ice thickness (column 8 assumed) ---
seaice_column <- 8
weekly_icethick <- tapply(df[[seaice_column]], df$WeekStart, mean, na.rm = TRUE)

# Maximum thickness for scaling
max_thick <- max(weekly_icethick, na.rm = TRUE)

####### PLOT WITH SEA ICE THICKNESS #######

par(mar = c(7, 4, 4, 4) + 0.1)

# --- Barplot of weekly presence ---
bp <- barplot(
  weekly_presence,
  names.arg = NA,
  col = "lightpink1",
  ylim = c(0, 7),
  ylab = "Weekly Presence (0–7)",
  main = paste("Weekly Presence of", species, "and Sea Ice Thickness at", site)
)

# --- X-axis labels (week-start dates) ---
week_dates <- as.Date(names(weekly_presence))
label_index <- seq(1, length(week_dates), by = 2)
axis(1, at = bp[label_index], labels = format(week_dates[label_index], "%b %d"),
     las = 2, cex.axis = 0.8)
mtext("Date (Month Day)", side = 1, line = 5)

# --- Scale thickness to match 0–7 y-axis ---
icethick_scaled <- weekly_icethick / max_thick * 7

# Overlay sea ice thickness line
par(new = TRUE)
plot(bp, icethick_scaled, type = "l", col = "deeppink", lwd = 2,
     axes = FALSE, xlab = "", ylab = "", ylim = c(0, 7))

# --- Right Y-axis for actual thickness values ---
axis(4,
     at = seq(0, 7, length.out = 6),
     labels = round(seq(0, max_thick, length.out = 6), 2),
     las = 2)
mtext("Sea Ice Thickness (m)", side = 4, line = 2.5)

