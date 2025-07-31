# ------------------------------------------------------------------------------
#                                  TO DEFINE
# ------------------------------------------------------------------------------
# Define the paths for scripts and data
scripts_path <- "C:/Users/Morgane/Desktop/Morgane_GOM/GAMs/GOMCetaceanMapper"   # Path where functions are located
data_path <- file.path(scripts_path,"0.Data")   # Data path
results_path <- "C:/Users/Morgane/Desktop/Morgane_GOM/GAMs/GOMCetaceanMapper_Morgane_Results"   # Results path, your hard disk

# Define the species to model
species <- "densUD"

# ------------------------------------------------------------------------------
#                                 LOAD DATA
# ------------------------------------------------------------------------------
# Load the dataset
data <- read.csv(file.path(data_path,
                           "AllSites_GOM_2020_2024_data_habitat_modeling_GLBy.csv"), 
                 header = TRUE)

# ------------------------------------------------------------------------------
#                         IMPORT INTERNAL FUNCTIONS
# ------------------------------------------------------------------------------
# Source internal functions
source(file.path(scripts_path, "4.Modeling", "check-functions.R"))
source(file.path(scripts_path, "4.Modeling", "split-functions.R"))
source(file.path(scripts_path, "4.Modeling", "calibration-functions.R"))
source(file.path(scripts_path, "4.Modeling", "evaluation-functions.R"))
source(file.path(scripts_path, "4.Modeling", "better_check.R"))

# ------------------------------------------------------------------------------
#                         INSTALL REQUIRED PACKAGES
# ------------------------------------------------------------------------------
required_packages <- c(
  "dismo", "gbm",                                       # BRT models and model selection
  "ggpubr", "ggplot2", "corrplot",                      # Plotting
  "parallel", "doParallel", "foreach", "tcltk",         # Parallel
  "pROC", "Hmisc", "dplyr", "Metrics", "PresenceAbsence","spdep",# Metrics: AUC, correlation
  "nortest",                                            # Extra Normality test
  "mgcv", "car",                                        # Multicollinearity check
  "EWSmethods",                                         # Decompose variables
  "lindesaysh/MRSea", "MuMIn",                          # SALSA1D
  "ncdf4", "raster", "sf",                              # General packages
  "thomasp85/scico", "tesselle/khroma", "hafen/stlplus",# GitHub repo for plotting
  "statmod","tweedie","gridBase","grid", "itsadug"      # GAM check for tweedie
)

# Install missing packages
.install_packages_if_missing(required_packages)

# ------------------------------------------------------------------------------
#                             CREATE FOLDERS
# ------------------------------------------------------------------------------
# Create folder if it doesn't exist to store results
if (!dir.exists(results_path)) {dir.create(results_path)}

# Create species subfolder
species_results_path <- file.path(results_path, species)
if (!dir.exists(species_results_path)) {dir.create(species_results_path)}

# ------------------------------------------------------------------------------
#                             DATA FORMATTING
# ------------------------------------------------------------------------------
# Convert data to the right format
data <- data %>%
  dplyr::mutate(across(everything(), ~ifelse(is.nan(.), NA, .))) %>%  # Replace nans with NA
  dplyr::mutate(
    Time = as.Date(Time, format = '%d-%b-%Y'),   # Convert 'Time' to Date format
    Site = as.factor(Site),                      # Convert 'Site' to a factor
    # typeEddy = as.factor(typeEddy),              # Convert 'typeEddy' to a factor # NONE IN THIS DATASET
    DepthRange = cut(Depth,                       
                     breaks = c(-Inf, -2500, -1500, -1000, -500, -200),
                     labels = c(">2500", "1500-2500", "1000-1500", "500-1000", "200-500"),
                     right = FALSE) %>% 
      factor(ordered = TRUE), # Treat as ordered factor
    SlopeRange = cut(Slope,                       
                     breaks = c(0, 1, 3, Inf),
                     labels = c("0-1", "1-3", ">3"),
                     right = FALSE) %>% 
      factor(ordered = TRUE)  # Treat as ordered factor
  )

# Remove rows where the specified species density is NA
dataSp <- data %>% 
  dplyr::filter(!is.na(.data[[species]]))

#Create new plotting window and create plot
dataSp <- data %>% 
  dplyr::filter(!is.na(.data[[species]]))
{windows(width = 8, height = 6)
  plot.new()
  dev.new(width = 5, height = 4, noRStudioGD = TRUE)
  D1 = ggplot2::ggplot(dataSp,aes(x=.data[[species]]))+
    ggplot2::geom_histogram(position="identity", fill="#6ABEDB",binwidth = 0.1,color="black")+
    ggplot2::ggtitle("Actual distribution")+
    theme_classic()
  
  #Output plot, saves as pdf, and closes graphics device
  print(D1)
  plot.File = file.path(species_results_path, paste(species,"Response_Transformation.pdf",sep="_"))
  dev.print(pdf, plot.File) 
  dev.off()
  }

# ------------------------------------------------------------------------------
#                     PLOT ORIGINAL COVARIATE DISTRIBUTIONS
# ------------------------------------------------------------------------------
# Create dataSp_New dataset with only densUD and with RowID added (row numbers 
# to the existing data to keep track once the data is reshuffled and manipulated 
# in log transformations and deseasoning)

# Identify all dens* columns except for "densUD"
cols_to_drop <- grep("^dens", names(dataSp), value = TRUE)
cols_to_drop <- setdiff(cols_to_drop, "densUD")  # keep densUD

# Apply transformation
dataSp_New_UD <- dataSp %>%
  dplyr::mutate(RowID = dplyr::row_number()) %>%
  dplyr::select(-dplyr::all_of(cols_to_drop)) %>%  # drop all dens* except densUD
  dplyr::filter(Time < as.Date("2024-06-01")) # remove all data after June 2024

pdf_file <- file.path(species_results_path, paste0(species, "_Covariate_Distributions_Original.pdf"))
pdf(pdf_file, width = 10, height = 8)
par(mfrow = c(5, 4), mar = c(4, 4, 2, 1))
for (v in names(dataSp_New_UD)) {
  if (v !="RowID") {
    plot(x = dataSp_New_UD$RowID,
         y = dataSp_New_UD[[v]],
         main = v, xlab = "RowID", ylab = v,
         pch = 19, cex = 0.3)
  }
}
dev.off()
# ------------------------------------------------------------------------------
#                           LOG TRANSFORM VARS
# ------------------------------------------------------------------------------
dataSp_New_UD <- dataSp_New_UD %>%
  dplyr::mutate(Chla_log = log10(Chla))

# ------------------------------------------------------------------------------
#                      SEASONAL TREND DECOMPOSITION
# ------------------------------------------------------------------------------
# This method of decomposition to remove seasonality will be used on SST, SSS, 
# and MXL. A daily-mean time series will be created for every variable - the 
# seasonality will be counted based on this GOM wide time series for each 
# variable and then will be removed from the original variable datasets. This 
# method is meant to prevent the addition of NAs in the deseasoned datasets.

.deseason_variable_average <- function(data, var, n.p = 365) {
  
  # calculate daily average across all sites
  global_daily_avg <- data %>%
    dplyr::group_by(Time) %>%
    dplyr::summarise(var_avg = mean(.data[[var]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(Time)
  
  # STL decomposition on global mean time series
  stl_result <- tryCatch({
    stlplus::stlplus(
      x = global_daily_avg$var_avg,
      t = global_daily_avg$Time,
      s.window = "periodic",
      n.p = n.p
    )
  }, error = function(e) NULL)
  if (is.null(stl_result)) {
    stop(paste("STL decomposition failed for", var))
  }
  
  # Create seasonal component table
  seasonality_df <- data.frame(Time = global_daily_avg$Time)
  seasonality_df[[paste0(var, "_seasonal_avg")]] <- stl_result$data$seasonal
  seasonality_df[[paste0(var, "_trend_avg")]]    <- stl_result$data$trend
  seasonality_df[[paste0(var, "_remainder_avg")]] <- stl_result$data$remainder
  
  # keeps he original average values for plotting
  seasonality_df[[paste0(var, "_avg")]] <- global_daily_avg$var_avg
  
  # Join seasonal components back to full data
  data_augmented <- data %>%
    dplyr::left_join(seasonality_df, by = "Time")
  
  # original variable - seasonal avg calculation 
  data_augmented[[paste0(var, "_deseasoned")]] <- 
    data_augmented[[var]] - data_augmented[[paste0(var, "_seasonal_avg")]]
  
  return(list(
    data = data_augmented,
    seasonal_table = seasonality_df
  ))
}

# applied to all three variables separately
SST_out <- .deseason_variable_average(dataSp_New_UD, "SST", n.p = 365)
MXL_out <- .deseason_variable_average(SST_out$data, "MXL", n.p = 365)
SSS_out <- .deseason_variable_average(MXL_out$data, "SSS", n.p = 365)

dataSp_New_UD <- SSS_out$data
SST_daily_avg <- SST_out$seasonal_table
MXL_daily_avg <- MXL_out$seasonal_table
SSS_daily_avg <- SSS_out$seasonal_table

# Reorder by RowID to restore original order
dataSp_New_UD <- dataSp_New_UD %>%
  dplyr::arrange(RowID)

# Set output PDF path
pdf_file <- file.path(species_results_path, paste0(species, "_SST_SSS_MXL_Decomposition_Comparison.pdf"))
pdf(pdf_file, width = 15, height = 9)  # 5 cols x 3 rows layout
par(mfrow = c(3, 5), mar = c(4, 4, 2, 1))  # 3 rows, 5 columns
# Helper function for plotting with consistent format
plot_with_rowid <- function(rowid, y, label, col = "black") {
  plot(rowid, y, type = "l", col = col,
       xlab = "RowID", ylab = label,
       main = label)
}
# Plot for SST (row 1)
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$SST, "SST (Original)")
plot_with_rowid(dataSp_New_UD$RowID, SST_daily_avg$SST_avg[match(dataSp_New_UD$Time, SST_daily_avg$Time)], "SST_avg (Daily Mean)", col = "black")
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$SST_seasonal_avg, "SST_seasonal_avg", col = "black")
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$SST_trend_avg, "SST_trend_avg", col = "black")
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$SST_deseasoned, "SST_deseasoned", col = "black")
# Plot for SSS (row 2)
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$SSS, "SSS (Original)")
plot_with_rowid(dataSp_New_UD$RowID, SSS_daily_avg$SSS_avg[match(dataSp_New_UD$Time, SSS_daily_avg$Time)], "SSS_avg (Daily Mean)", col = "black")
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$SSS_seasonal_avg, "SSS_seasonal_avg", col = "black")
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$SSS_trend_avg, "SSS_trend_avg", col = "black")
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$SSS_deseasoned, "SSS_deseasoned", col = "black")
# Plot for MXL (row 3)
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$MXL, "MXL (Original)")
plot_with_rowid(dataSp_New_UD$RowID, MXL_daily_avg$MXL_avg[match(dataSp_New_UD$Time, MXL_daily_avg$Time)], "MXL_avg (Daily Mean)", col = "black")
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$MXL_seasonal_avg, "MXL_seasonal_avg", col = "black")
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$MXL_trend_avg, "MXL_trend_avg", col = "black")
plot_with_rowid(dataSp_New_UD$RowID, dataSp_New_UD$MXL_deseasoned, "MXL_deseasoned", col = "black")
dev.off()

# ------------------------------------------------------------------------------
# ADDED AS TEST TO SEE IF THIS WOULD IMPROVE PREDICTIONS (PRIOR TO SPLIT BATHY/MESO): 
# FOR EACH SITE THAT IS NOT AT THE DEEPEST DEPTHRANGE, CHANGED NAS TO ZEROS FOR 
# ALL KE AND VORT VARIABLES THAT ARE DEEPER THAN THE DEPTH OF THE SITE
# -----------------------------------------------------------------------------
zero_fill_by_depthrange <- function(data) {
  # Map from DepthRange to numeric depth in meters
  depthrange_map <- c(
    "200-500" = 500,
    "500-1000" = 1000,
    "1000-1500" = 1500,
    "1500-2500" = 2500,
    ">2500" = 3000
  )
  
  # Depths for each KE and Vort variable
  depth_by_var <- list(
    "KE_z300_600"     = 600,
    "KE_z700_1250"    = 1250,
    "KE_z1500_3000"   = 3000,
    "Vort_z0_250"     = 250,
    "Vort_z300_600"   = 600,
    "Vort_z700_1250"  = 1250,
    "Vort_z1500_3000" = 3000
  )
  
  # Loop through each site
  unique_sites <- unique(data$Site)
  
  for (site in unique_sites) {
    site_data <- data[data$Site == site, ]
    site_depthrange <- unique(site_data$DepthRange)
    
    # convert to character to avoid factor key mismatch
    site_max_depth <- depthrange_map[as.character(site_depthrange)]
    if (length(site_depthrange) != 1 || is.na(site_max_depth)) next
    
    for (var in names(depth_by_var)) {
      if (!var %in% names(data)) next
      var_depth <- depth_by_var[[var]]
      
      # Replace with 0 only if var is deeper and all values are NA
      if (var_depth > site_max_depth) {
        data[data$Site == site & is.na(data[[var]]), var] <- 0
        message(sprintf("Set %s to 0 for site %s (depth range: %s)", var, site, site_depthrange))
      }
    }
  }
  
  return(data)
}

# applies the function to dataSp_New_UD
dataSp_New_UD <- zero_fill_by_depthrange(dataSp_New_UD)

# ------------------------------------------------------------------------------
#                         SPLIT DATA INTO TRAIN/TEST
#                  (Temporal and Spatial Block Cross-Validation)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# TEMPORAL SPLITTING: Split the data into train and test sets for each site
# while ensuring that the test segments meet the minimum and maximum
# week criteria of continuous data.
#
# Constants:
# - ratioTest: Proportion of the data to allocate to the test set
# - minWeeks: Minimum number of continuous weeks a segment must span to be
#             considered valid for testing.
# - maxWeeks: Maximum number of continuous weeks a test segment can span.
# ------------------------------------------------------------------------------
# Constants for splitting data
yearTest = 2024
sitesTest = c("Y4A","Y4B","Y4D")
split_result <- .train_test_split_year_time_series(dataSp_New_UD, yearTest, sitesTest)
# Create summary table of proportions allocated to "Train" and "Test" per site
summary_set_splits <- split_result %>%
  dplyr::group_by(Site, Set) %>%
  dplyr::summarise(SetTotalDays = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(Percentage = (SetTotalDays / sum(SetTotalDays)) * 100) %>%
  dplyr::ungroup()
# Update main data frame with 'Set' column
dataSp_New_UD$Set = split_result$Set
# Separate data into train and test sets
train_data <- dataSp_New_UD %>% dplyr::filter(Set == "Train")
test_data <- dataSp_New_UD %>% dplyr::filter(Set == "Test")
# PLOT: Visualize train and test set distribution per site
{plot.dims = c(16,9)
  plot.new()
  dev.new(width = plot.dims[1], height = plot.dims[2], noRStudioGD = TRUE)
  Split.plot = ggplot2::ggplot(dataSp_New_UD, ggplot2::aes(x = Time, y=.data[[species]], color = Set)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::facet_wrap(~ Site, scales = "free_y") +  # Facet by site
    ggplot2::scale_color_manual(values = c("Train" = "#6ABEDB", "Test" = "#E39400")) +
    ggplot2::labs(title = "Train and Test Sets Across Sites",
                  x = NULL,
                  y = "Species Density",
                  color = NULL) +
    ggplot2::scale_x_date(date_labels = "%b'%y") +  # Format x-axis labels as Month'YY
    ggpubr::theme_pubr() +
    ggplot2::theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 9, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 9),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA)
    )
  print(Split.plot)
  plot.File = file.path(species_results_path, paste(species,"Train_Test_Site_Splits.pdf",sep="_"))
  dev.print(pdf, plot.File)
  dev.off()
}

# ------------------------------------------------------------------------------
#                           PREDICTORS - deseasoned
#                    Select predictors with no correlation, 
#             based on Spearman's rank correlation coefficient
# ------------------------------------------------------------------------------
# Select all variables not related to density, which correspond to enviro vars
numeric_columns <- train_data %>%
  dplyr::select(-SST_trend_avg, -SST_remainder_avg, -SST_avg,
                -SSS_seasonal_avg, -SSS_trend_avg, -SSS_remainder_avg, -SSS_avg,
                -MXL_seasonal_avg, -MXL_trend_avg, -MXL_remainder_avg, -MXL_avg,
                -SST, -SSS, -MXL, -RowID, -Latitude, -Longitude, -Depth, -Slope, 
                -Chla) %>%
  dplyr::select_if(is.numeric) %>%
  dplyr::select(-matches("^dens")) %>%
  dplyr::filter(if_all(everything(), ~ !is.na(.) & is.finite(.)))

# Compute the Spearman correlation matrix and corresponding p-values
spearman_corr <- Hmisc::rcorr(as.matrix(numeric_columns), type = "spearman")

# PLOT: Visualize correlation between predictor variables
{plot.new()
  dev.new(width = 8,height = 8,noRStudioGD = TRUE)
  corrplot::corrplot(spearman_corr$r, # correlation matrix
                     method = "ellipse",
                     order = "original",
                     addgrid.col = "darkgray",
                     col = colorRampPalette(c("blue","cyan","gray", "white","gray","yellow", "red"))(10),
                     tl.cex = 0.8,  # text label size
                     cl.pos = "r",  # legend position (bottom)
                     tl.col = "black",  # text label color
                     tl.srt = 90,  # Text label rotation
                     p.mat = spearman_corr$P, # matrix of p-values for significance testing
                     sig.level = 0.05,  # p-value significance
                     insig = "pch",  # mark insignificant correlations with crosses
                     pch.col = "black",  # cross color
                     pch.cex = 1,# cross size
                     diag = F, type = "upper")
  plot.File = file.path(species_results_path, paste(species,"Correlation_Predictor_Variables_Deseasoned.pdf",sep="_"))
  dev.print(pdf, plot.File) 
  dev.off()
}

selected.predictors <- c(
  "DepthRange", "SlopeRange",                         # Group 1; "Latitude", "Longitude", "Depth", "Slope"
  "MoonIllu",                                         # Group 2: "MoonPhase",
  #                                                   # Group 3 (all excluded): "LCext", "LCmaxN", "LCmaxW" 
  "SST_deseasoned", "SSS_deseasoned", "SST_seasonal_avg", "Chla_log", "MXL_deseasoned", "UI", # Group 4
  #                                                   # Group 5 (all excluded): "distEddy", "ampEddy", "vortEddy",  "MESI","diverEddy", "EKE", 
  "Sigma22Z", "Sigma26Z",                             # Group 6
  "KE_z300_600", "KE_z700_1250", "KE_z1500_3000",       # Group 7: "KE_z0_250", 
  "Vort_z0_250", "Vort_z300_600", "Vort_z700_1250", "Vort_z1500_3000"  # Group 8
)

# ------------------------------------------------------------------------------
#                            GAM PRE-CHECKS
# ------------------------------------------------------------------------------

# pre-GAM modeling checks:
# - check for normal distribution of response variable
# - check autocorrelation (spatial and temporal)
# - check multicollinearity and pairwise correlation
# - select smooths
# - (if not normal) find appropriate distribution (family) - tweedie distribution

# ------------------------------------------------------------------------------
#                       TEST FOR NORMAL DISTRIBUTION
# ------------------------------------------------------------------------------
# testing the distribution of our response variable
# Shapiro-Wilk's test:
# Ho: data is normal
# Halt: data isn't normal
#shapiro_test <- stats::shapiro.test(dataSp$densUD) # length of dataset is 9983, too large for shapiro

# random sample of 5000 for shapiro test
shapiro_test <- stats::shapiro.test(sample(train_data$densUD,5000))
print(shapiro_test) 
# result of shapiro-wilk test with random sample -> not normally distributed therefore tweedie distribution justified
# data:  sample(train_data$densUD, 5000)
# W = 0.64786, p-value < 2.2e-16

# Anderson-Darling normality test (for really large sample size): 
# NOTE: if we do have a really large dataset and want to use Shapiro-Wilk's, you can take a subset of your data and run SW
anderson_darling_test <- nortest::ad.test(train_data$densUD)
print(anderson_darling_test)
# result of anderson-darling test -> tweedie distribution justified, not normal distribution
# data:  train_data$densUD
# A = 1416.4, p-value < 2.2e-16

# ------------------------------------------------------------------------------
#                         CALCULATE TEMPORAL AUTOCORRELATION
# ------------------------------------------------------------------------------
# Estimate ACF by site for predictor diagnostics
pdf_file <- file.path(species_results_path, paste0(species, "_ACF_densUD_by_Site.pdf"))
pdf(pdf_file, width = 12, height = 10)
sites <- as.character(unique(train_data$Site))
n_sites <- length(sites)
n_cols <- 7
n_rows <- ceiling(n_sites / n_cols)
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

for (s in sites) {
  site_data <- train_data$densUD[train_data$Site == s]
  if (sum(!is.na(site_data)) > 1) {
    acf(site_data, main = paste(s),
        xlab = "Lag", ylab = "ACF", na.action = na.pass)
  } else {
    plot.new()
    title(main = paste("Site:", s, "- insufficient data"))
  }
}
dev.off()
# ------------------------------------------------------------------------------
#                         TEST SPATIAL AUTOCORRELATION
#                                 Moran's I
# ------------------------------------------------------------------------------
# Define spatial coordinates 
dataSp_spatialcorr_UD <- train_data %>%
  dplyr::filter(!is.na(densUD) & !is.na(Longitude) & !is.na(Latitude))
coords <- cbind(dataSp_spatialcorr_UD$Longitude, dataSp_spatialcorr_UD$Latitude)

# Check for duplicate coordinates
if (sum(duplicated(coords)) > 0) {
  warning("Duplicate coordinates detected! Switching to distance-based neighbors.")
  
  # Use distance-based neighbors instead of k-nearest neighbors
  nb <- spdep::dnearneigh(coords, d1 = 0, d2 = 1, longlat = TRUE)
} else {
  # Create a k-nearest neighbors list
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = 5))
}
# Output: Duplicate coordinates detected! Switching to distance-based neighbors.

# Convert neighbors list to spatial weights (lw)
lw <- spdep::nb2listw(nb, style = "W") # style W means that the weights are row-standardized

# Ensure densUD and lw have the same length
if (length(dataSp_spatialcorr_UD$densUD) != length(lw$neighbours)) {
  stop("Error: Response variable and spatial weights have different lengths. Check for missing values in coordinates.")
}

# Run Moran’s I test on response variable
moran_test <- spdep::moran.test(dataSp_spatialcorr_UD$densUD, lw)

# Print Moran’s I test results
print(moran_test)

# Moran I test under randomisation
# 
# data:  dataSp_spatialcorr_UD$densUD  
# weights: lw    
# 
# Moran I statistic standard deviate = 224.24, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
#       1.105456e-01     -8.237232e-05      2.433965e-07 

#   POSITIVE AUTOCORRELATION since Moran I > 0, Moran I = 0.111
#   look at the p-value, <0.05 = spatial autocorrelation
#   p-value = 2.2e-16, statistically significant!
#   Variance is small, therefore, test result is stable and reliable
#   Therefore, strong spatial autocorrelation, statistically significant

# ------------------------------------------------------------------------------
#                    SPLIT train_data INTO BATHY AND MESO
# ------------------------------------------------------------------------------
# Plot depths of each site to see spread
plot(train_data$Site,train_data$Depth) # decided on 2 groupings: bathy (>=1500m), meso (<1500)

# bathy site grouping
bathySites <- train_data %>%
  dplyr::filter(abs(Depth) >= 1500)

# meso site grouping
mesoSites <- train_data %>%
  dplyr::filter(abs(Depth) < 1500) %>%
  dplyr::select(-dplyr::any_of(c("KE_z1500_3000", "Vort_z1500_3000"))) # Drop deep KE/Vort variables from mesoSites (all deeper than 1500m)

plot(bathySites$Site,bathySites$Depth) # plot of bathysites vs depth
plot(mesoSites$Site,mesoSites$Depth)  # plot of mesosites vs depth

# Function to preprocess data for multicollinearity checks and GAMs 
# - deseasoned version
.prepare_data <- function(data) {
  
  # Remove non-numeric and unselected variables (including all that start with 
  # dens except densUD) "data" is used here to be able to run bathy and meso 
  # separately after the function is defined
  # Keep densUD and selected predictors
  filtered_data_with_factors <- data %>%
    dplyr::select(any_of(c("Site", "Time", "DepthRange", "SlopeRange", "densUD", selected.predictors)))
  
  # Remove columns with all NAs
  cols_with_all_na <- names(which(colSums(is.na(filtered_data_with_factors)) == nrow(filtered_data_with_factors)))
  filtered_data_with_factors <- filtered_data_with_factors %>%
    dplyr::select(-all_of(cols_with_all_na))
  
  # Calculates the number of of NA values for each column in the dataset
  filtered_data_with_factors %>% dplyr::summarise(across(everything(), ~ sum(is.na(.))))
  
  # Remove rows with any remaining NAs
  filtered_data_with_factors <-filtered_data_with_factors %>%
    tidyr::drop_na() 
  
  # Make dataset without NAs and factors (only numeric) for GVIF
  filtered_datasets_for_GVIF <-filtered_data_with_factors %>%
    dplyr::select_if(is.numeric) %>%
    dplyr::select(-dplyr::any_of("densUD")) # Removes densUD
  
  # Return both datasets (one for GAM, one for GVIF)
  return(list(filtered_data_with_factors, filtered_datasets_for_GVIF))
}

# Run the function for bathy and meso datasets
bathy_results <- .prepare_data(bathySites)
meso_results <- .prepare_data(mesoSites)

# Store cleaned datasets for GAMs
bathy_filtered_data_with_factors <- bathy_results[[1]]
meso_filtered_data_with_factors <- meso_results[[1]]

# Store cleaned numeric-only datasets for GVIF
bathy_filtered_datasets_for_GVIF <- bathy_results[[2]]
meso_filtered_datasets_for_GVIF <- meso_results[[2]]

# ------------------------------------------------------------------------------
#                        PLOT COVARIATE DISTRIBUTIONS
# ------------------------------------------------------------------------------
# BATHY DESEASONED
pdf_file <- file.path(species_results_path, paste0(species, "_Covariate_Distributions_Bathy_deseasoned.pdf"))
pdf(pdf_file, width = 10, height = 8)
par(mfrow = c(5, 4), mar = c(4, 4, 2, 1))
for (v in names(bathy_filtered_datasets_for_GVIF)) {
  plot(x = seq_len(nrow(bathy_filtered_datasets_for_GVIF)),
       y = bathy_filtered_datasets_for_GVIF[[v]],
       main = v, xlab = "Row Number", ylab = v,
       pch = 19, cex = 0.3)
}
dev.off()

# MESO DESEASONED
pdf_file <- file.path(species_results_path, paste0(species, "_Covariate_Distributions_Meso_deseasoned.pdf"))
pdf(pdf_file, width = 10, height = 8)
par(mfrow = c(5, 4), mar = c(4, 4, 2, 1))
for (v in names(meso_filtered_datasets_for_GVIF)) {
  plot(x = seq_len(nrow(meso_filtered_datasets_for_GVIF)),
       y = meso_filtered_datasets_for_GVIF[[v]],
       main = v, xlab = "Row Number", ylab = v,
       pch = 19, cex = 0.3)
}
dev.off()

# ------------------------------------------------------------------------------
#                         CHECK MULTICOLINEARITY
# ------------------------------------------------------------------------------
#GVIF/VIF - use both because categorical and continuous predictor variables
# NOTE: now we want to see how do each explanatory variables correlate with all 
# of the other explanatory variables, perhaps in individual pairs we get low 
# colinearity, but if a variable is slightly correlated with a lot of our 
# variables that could be problematic. here we test for multicolinearity

# CODE FOR GVIF/VIF FUNCTIONS:
# Github for code: https://gist.github.com/jfroeschke/a7a900b41c7348b3c0b30a36766223e5
#Library files for courses provided by: Highland Statistics Ltd.
#To cite these functions, use:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.
#Copyright Highland Statistics LTD.

# Compute the Spearman correlation matrix and corresponding p-values
.compute_correlation <- function(data, label) {
  data <- data %>% dplyr::select(-dplyr::any_of(c("densUD", "RowID"))) # Remove densUD from corrplots, since it is the response variable
  spearman_corr <- Hmisc::rcorr(as.matrix(data), type = "spearman")
  
  # PLOT: Visualize correlation between predictor variables
  plot.new()
  dev.new(width = 8,height = 8,noRStudioGD = TRUE)
  corrplot::corrplot(spearman_corr$r, # correlation matrix
                     method = "ellipse",
                     order = "original",
                     addgrid.col = "darkgray",
                     col = colorRampPalette(c("blue","cyan","gray", "white","gray","yellow", "red"))(10),
                     tl.cex = 0.8,  # text label size
                     cl.pos = "r",  # legend position (bottom)
                     tl.col = "black",  # text label color
                     tl.srt = 90,  # Text label rotation
                     p.mat = spearman_corr$P, # matrix of p-values for significance testing
                     sig.level = 0.05,  # p-value significance
                     insig = "pch",  # mark insignificant correlations with crosses
                     pch.col = "black",  # cross color
                     pch.cex = 1,# cross size
                     diag = F, type = "upper",
                     addCoef.col = "black", # Display correlation coefficients
                     number.cex = 0.5)  # Adjust text size for correlation values
  title(main = paste("Correlation Plot for", label, "Sites"))
  
  # Save plots
  plot_file <- file.path(species_results_path, paste(species, "Correlation_Plot_for", label, "Sites.pdf", sep = "_"))
  dev.print(pdf, plot_file)
  dev.off()
}

# Generate correlation plots
if (!is.null(bathy_filtered_datasets_for_GVIF)) .compute_correlation(bathy_filtered_datasets_for_GVIF, "bathy") # deseasoned vars
if (!is.null(meso_filtered_datasets_for_GVIF)) .compute_correlation(meso_filtered_datasets_for_GVIF, "meso") # deseasoned vars

# Function for computing GVIF/VIF
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: VIFs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1]) {
    diag(tmp_cor) <- 0
    if (any(tmp_cor == 1.0)) {
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  # If all predictors have Df = 1, then only GVIF should be returned
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF = result[, 1])  # Only show GVIF
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2])) # Compute GVIF^(1/(2Df))
  }
  print("VIF Results from myvif():")
  print(result)
  
  return(result)
}

# Function to run multicollinearity check
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  # Create the formula dynamically for the model
  form <- formula(paste("fooy ~ ", paste(names(dataz), collapse = " + ")))
  
  # Create a new dataset with a dummy response variable 'fooy'
  dataz <- data.frame(fooy = 1, dataz)
  
  # Fit the linear model
  lm_mod <- lm(form, dataz)
  
  vif_results <- myvif(lm_mod)
}

# Run VIF/Multicollinearity check for deseasoned variables
message(paste("Running VIF check for bathy"))
tryCatch({
  vif_results_bathy <- corvif(bathy_filtered_datasets_for_GVIF) 
  if (any(is.na(vif_results_bathy))) {
    message("VIF computation failed due to singularities or insufficient data.")
  }
}, error = function(e) {
  message(paste("VIF check failed for bathy with error:", e$message))
})
message(paste("Running VIF check for meso"))
tryCatch({
  vif_results_meso <- corvif(meso_filtered_datasets_for_GVIF)
  if (any(is.na(vif_results_meso))) {
    message("VIF computation failed due to singularities or insufficient data.")
  }
}, error = function(e) {
  message(paste("VIF check failed for depth range meso with error:", e$message))
})

# ------------------------------------------------------------------------------
#                 LOOP TO REMOVE HIGHLY COLLINEAR VARIABLES
#             Need to run this loop even if no variables have GVIF >5
# ------------------------------------------------------------------------------
# starting with random high number (this variable would contain the highest vif 
# in each run, so this value here needs to be anything above 5 to start the loop)
filter_gvif <- function(data) {
  vifMax <- 5000
  threshVIF = 5
  selectedCovariates = setdiff(names(data), "RowID")
  
  while (vifMax > threshVIF){
    
    # compute vif for the selected predictors
    viftemp = corvif(data[selectedCovariates])
    #print(data.frame(viftemp),digits=3)
    
    # find max vif 
    vifMaxIndx = match (max(viftemp[,1]), viftemp[,1])
    vifMax = viftemp[vifMaxIndx, 1]
    if (vifMax>threshVIF){   
      selectedCovariates <- selectedCovariates[-vifMaxIndx]
    }
  }
  # Return the final selected variables
  return(selectedCovariates)
}

# Run GVIF filtering separately for bathy and meso sites, stores list of 
# final vars that survived the GVIF filtering
selectedCovariates_bathy <- filter_gvif(bathy_filtered_datasets_for_GVIF)
selectedCovariates_meso <- filter_gvif(meso_filtered_datasets_for_GVIF)

# ------------------------------------------------------------------------------
#              SPLIT TEST_DATA INTO BATHY AND MESO BEFORE PREDICT.BAM
# ------------------------------------------------------------------------------
# Bathy test data (Depth >= 1500m)
test_data_bathy <- test_data %>%
  dplyr::filter(abs(Depth) >= 1500)

# Meso test data (Depth < 1500m)
test_data_meso <- test_data %>%
  dplyr::filter(abs(Depth) < 1500) %>%
  dplyr::select(-dplyr::any_of(c("KE_z1500_3000", "Vort_z1500_3000")))

# ------------------------------------------------------------------------------
#         ensure that the factcors have same levels in test vs train data
# ------------------------------------------------------------------------------
test_data_bathy <- test_data_bathy %>%
  dplyr::mutate(
    DepthRange = factor(DepthRange, levels = levels(train_data$DepthRange)),
    SlopeRange = factor(SlopeRange, levels = levels(train_data$SlopeRange))
  ) 

test_data_meso <- test_data_meso %>%
  dplyr::mutate(
    DepthRange = factor(DepthRange, levels = levels(train_data$DepthRange)),
    SlopeRange = factor(SlopeRange, levels = levels(train_data$SlopeRange))
  ) 

# ------------------------------------------------------------------------------
#                        PLOT COVARIATE DISTRIBUTIONS 
#                                   TEST DATA
# ------------------------------------------------------------------------------
# BATHY DESEASONED
# Variables to skip manually
manual_skip <- c("SST", "Chla", "SSS","MXL","Latitude", "Longitude", "Depth", 
                 "densUD",  "Slope", "RowID", "SST_remainder_avg", "SST_avg", 
                 "MXL_seasonal_avg", "MXL_trend_avg", "MXL_remainder_avg", 
                 "MXL_avg", "SSS_seasonal_avg", "SSS_trend_avg", "SSS_remainder_avg",
                 "SSS_avg", "MoonPhase")  # Add any vars you want to ignore

pdf_file <- file.path(species_results_path, paste0(species, "_Covariate_Distributions_Bathy_test.pdf"))
grDevices::pdf(pdf_file, width = 10, height = 8)
graphics::par(mfrow = c(5, 4), mar = c(4, 4, 2, 1))

# Initialize trackers
skipped_vars <- c()
manually_skipped <- c()

for (v in names(test_data_bathy)) {
  if (v %in% manual_skip) {
    manually_skipped <- c(manually_skipped, v)
    next
  }
  
  vec <- test_data_bathy[[v]]
  
  # Only plot if numeric and has finite values
  if (is.numeric(vec) && any(is.finite(vec))) {
    graphics::plot(
      x = seq_len(length(vec)),
      y = vec,
      main = v, xlab = "Row Number", ylab = v,
      pch = 19, cex = 0.3
    )
  } else {
    skipped_vars <- c(skipped_vars, v)
  }
}

grDevices::dev.off()

# Print results
if (length(skipped_vars) > 0) {
  message("Skipped due to all NA or non-numeric: ", paste(skipped_vars, collapse = ", "))
} else {
  message("All variables plotted successfully.")
}


# MESO DESEASONED
# Variables to skip manually
manual_skip <- c("SST", "Chla", "SSS","MXL","Latitude", "Longitude", "Depth", 
                 "densUD", "Slope", "RowID", "SST_remainder_avg", "SST_avg", 
                 "MXL_seasonal_avg", "MXL_trend_avg", "MXL_remainder_avg", 
                 "MXL_avg", "SSS_seasonal_avg", "SSS_trend_avg", "SSS_remainder_avg",
                 "SSS_avg", "MoonPhase")  # Add any vars you want to ignore

pdf_file <- file.path(species_results_path, paste0(species, "_Covariate_Distributions_Meso_test.pdf"))
grDevices::pdf(pdf_file, width = 10, height = 8)
graphics::par(mfrow = c(5, 4), mar = c(4, 4, 2, 1))

# Initialize trackers
skipped_vars <- c()
manually_skipped <- c()

for (v in names(test_data_meso)) {
  if (v %in% manual_skip) {
    manually_skipped <- c(manually_skipped, v)
    next
  }
  
  vec <- test_data_meso[[v]]
  
  # Only plot if numeric and has finite values
  if (is.numeric(vec) && any(is.finite(vec))) {
    graphics::plot(
      x = seq_len(length(vec)),
      y = vec,
      main = v, xlab = "Row Number", ylab = v,
      pch = 19, cex = 0.3
    )
  } else {
    skipped_vars <- c(skipped_vars, v)
  }
}

grDevices::dev.off()

# Print results
if (length(skipped_vars) > 0) {
  message("Skipped due to all NA or non-numeric: ", paste(skipped_vars, collapse = ", "))
} else {
  message("All variables plotted successfully.")
}

# FITTING MODELS AND THEN PREDICTING FOR DENSUD!! ONE MODEL AT A TIME
# ------------------------------------------------------------------------------
#                     1. BAM MODEL (TWEEDIE) for BATHY SITES
#                   baseline model without autocorrelation
#                    with DepthRange and SlopeRange
# ------------------------------------------------------------------------------
# When running gamm(), use corAR1(form= ~ Time | Site) is best 
# since the data is hierarchical corAR1 is designed for grouped data, modeling 
# correlation of residuals across time within groups, and works directly into GAMMs using REML
# Instead, running bam() because it supports tweedie with AR1 structure

gc() # Cleans memory before running the model iterations, optional

# Prep bathy dataset
bathy_gam_data <- bathy_filtered_data_with_factors %>%
  dplyr::select(Site, Time, DepthRange, SlopeRange, densUD, SST_seasonal_avg, all_of(selectedCovariates_bathy)) %>% 
  stats::na.omit()

kVal <- 9
smooth_terms <- paste0(
  "s(", selectedCovariates_bathy[selectedCovariates_bathy != "SST_seasonal_avg"], ", bs = 'ts', k= ", kVal, ")")
# add SST seasonal avg with cyclic spline
seasonal_term <- "s(SST_seasonal_avg, bs = 'cc', k = 9)"
all_smooths <- c(smooth_terms, seasonal_term)

# construct the formula
formula_text <- paste(
  "densUD ~ DepthRange + SlopeRange +",
  paste(all_smooths, collapse = "+")
)
formulaModel <- as.formula(formula_text)

# Fit bam Tweedie model for bathy sites, with site-level random effects
FullModel_bathy <- mgcv::bam(
  formula = formulaModel,
  family = mgcv::tw(),
  data = bathy_gam_data,
  method = "fREML",
  select = TRUE,
  AR.start = NULL
)

namesCovariates_bathy = selectedCovariates_bathy
iVarList = seq_along (selectedCovariates_bathy) # for numbering when going through selections
NVars = length(selectedCovariates_bathy)

# based on baseline model, determine the value of lag 1 = Rho to try in next model
r1 <- itsadug::start_value_rho(FullModel_bathy, plot=TRUE)
acf(resid(FullModel_bathy), plot=FALSE)$acf[2]
# [1] 0.01323486 = 0.013

# ------------------------------------------------------------------------------
#                     2. BAM MODEL (TWEEDIE) for BATHY SITES
#             fit with autocorrelation AR1 style using Rho + AR.start
#                               with factors
# ------------------------------------------------------------------------------
# Prep bathy dataset
bathy_gam_data <- bathy_filtered_data_with_factors %>%
  dplyr::select(Site, Time, DepthRange, SlopeRange, densUD, SST_seasonal_avg, all_of(selectedCovariates_bathy)) %>% 
  stats::na.omit() %>%
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

kVal <- 9
smooth_terms <- paste0(
  "s(", selectedCovariates_bathy[selectedCovariates_bathy != "SST_seasonal_avg"], ", bs = 'ts', k= ", kVal, ")")
# add SST seasonal avg with cyclic spline
seasonal_term <- "s(SST_seasonal_avg, bs = 'cc', k = 9)"
all_smooths <- c(smooth_terms, seasonal_term)

# Model formula
formula_text <- paste(
  "densUD ~ DepthRange + SlopeRange +",
  paste(all_smooths, collapse = "+")
)
formulaModel <- as.formula(formula_text)

# Fit bam Tweedie model for bathy sites with AR1 autocorrelation
FullModel_bathy_Rho0.013 <- mgcv::bam(
  formula = formulaModel,
  family = mgcv::tw(),
  data = bathy_gam_data,
  method = "fREML",
  select = TRUE,
  discrete = TRUE,
  rho = 0.013,
  AR.start = bathy_gam_data$AR.start
)

# Check if AR.start is working right, TRUE at start of each site by checking bathy_gam_data added column
# good

namesCovariates_bathy = selectedCovariates_bathy
iVarList = seq_along (selectedCovariates_bathy) # for numbering when going through selections
NVars = length(selectedCovariates_bathy)

summary(FullModel_bathy_Rho0.013)
# Family: Tweedie(p=1.387) 
# Link function: log 
# 
# Formula:
#   densUD ~ DepthRange + SlopeRange + s(MoonIllu, bs = "ts", k = 9) + 
#   s(SST_deseasoned, bs = "ts", k = 9) + s(SSS_deseasoned, bs = "ts", 
#                                           k = 9) + s(Chla_log, bs = "ts", k = 9) + s(MXL_deseasoned, 
#                                                                                      bs = "ts", k = 9) + s(UI, bs = "ts", k = 9) + s(Sigma22Z, 
#                                                                                                                                      bs = "ts", k = 9) + s(Sigma26Z, bs = "ts", k = 9) + s(KE_z300_600, 
#                                                                                                                                                                                            bs = "ts", k = 9) + s(KE_z700_1250, bs = "ts", k = 9) + s(KE_z1500_3000, 
#                                                                                                                                                                                                                                                      bs = "ts", k = 9) + s(Vort_z0_250, bs = "ts", k = 9) + s(Vort_z300_600, 
#                                                                                                                                                                                                                                                                                                               bs = "ts", k = 9) + s(Vort_z700_1250, bs = "ts", k = 9) + 
#   s(Vort_z1500_3000, bs = "ts", k = 9) + s(SST_seasonal_avg, 
#                                            bs = "cc", k = 9)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    5.4426     0.1050  51.857  < 2e-16 ***
#   DepthRange.L   0.0000     0.0000     NaN      NaN    
# SlopeRange.L  -0.5019     0.0667  -7.526 8.54e-14 ***
#   SlopeRange.Q   0.4812     0.1644   2.926  0.00348 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
#   s(MoonIllu)         1.426e-05      8 0.000  0.8916    
#   s(SST_deseasoned)   5.254e-01      8 0.136  0.1436    
#   s(SSS_deseasoned)   4.418e-01      8 0.098  0.1776    
#   s(Chla_log)         3.866e+00      8 5.367  <2e-16 ***
#   s(MXL_deseasoned)   8.229e-01      8 0.510  0.0225 *  
#   s(UI)               2.242e-05      8 0.000  0.5666    
#   s(Sigma22Z)         5.174e-05      8 0.000  0.4475    
#   s(Sigma26Z)         3.363e+00      8 4.960  <2e-16 ***
#   s(KE_z300_600)      3.307e+00      8 4.071  <2e-16 ***
#   s(KE_z700_1250)     2.627e-05      8 0.000  0.7934    
#   s(KE_z1500_3000)    1.703e+00      8 0.535  0.0666 .  
#   s(Vort_z0_250)      7.135e-01      8 0.300  0.0614 .  
#   s(Vort_z300_600)    1.818e-05      8 0.000  0.9252    
#   s(Vort_z700_1250)   1.540e-05      8 0.000  0.6228    
#   s(Vort_z1500_3000)  3.440e-01      8 0.065  0.2148    
#   s(SST_seasonal_avg) 1.626e-05      7 0.000  0.9318    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 130/131
# R-sq.(adj) =  0.139   Deviance explained = 15.6%
# fREML = 6009.7  Scale est. = 67.368    n = 1684

AIC(FullModel_bathy_Rho0.013)

# shrinkage for model 2
# Drop the first non-significant smooth term (p > 0.2)
vars_to_drop <- c("Chla_log", "SST_deseasoned", "MoonIllu", "UI", "KE_z700_1250", "KE_z1500_3000",
                  "Vort_z300_600", "Vort_z700_1250",  "Vort_z1500_3000")

pattern <- paste(vars_to_drop, collapse = "|")

all_smooths_reduced <- all_smooths[!grepl(pattern, all_smooths)]

formula_text_reduced <- paste(
  "densUD ~ DepthRange + SlopeRange +",
  paste(all_smooths_reduced, collapse = "+")
)
formulaModel_reduced <- as.formula(formula_text_reduced)

# Refit the model
ReducedModel_bathy_Rho0.013 <- mgcv::bam(
  formula = formulaModel_reduced,
  family = mgcv::tw(),
  data = bathy_gam_data,
  method = "fREML",
  select = TRUE,
  discrete = TRUE,
  rho = 0.013,
  AR.start = bathy_gam_data$AR.start
)

summary(ReducedModel_bathy_Rho0.013)
# Family: Tweedie(p=1.393) 
# Link function: log 
# 
# Formula:
#   densUD ~ DepthRange + SlopeRange + s(SSS_deseasoned, bs = "ts", 
#                                        k = 9) + s(MXL_deseasoned, bs = "ts", k = 9) + s(Sigma22Z, 
#                                                                                         bs = "ts", k = 9) + s(Sigma26Z, bs = "ts", k = 9) + s(KE_z300_600, 
#                                                                                                                                               bs = "ts", k = 9) + s(Vort_z0_250, bs = "ts", k = 9) + s(SST_seasonal_avg, 
#                                                                                                                                                                                                        bs = "cc", k = 9)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.971e+00  2.859e-01  17.389   <2e-16 ***
#   DepthRange.L -4.850e-01  0.000e+00    -Inf   <2e-16 ***
#   SlopeRange.L -5.782e-08  6.483e-02   0.000    1.000    
# SlopeRange.Q  1.536e-01  1.639e-01   0.937    0.349    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
#   s(SSS_deseasoned)   1.0063      8 2.278 5.68e-06 ***
#   s(MXL_deseasoned)   0.7211      8 0.307  0.06240 .  
#   s(Sigma22Z)         0.8869      8 0.761  0.00613 ** 
#   s(Sigma26Z)         3.2339      8 4.597  < 2e-16 ***
#   s(KE_z300_600)      3.0383      8 3.196 2.77e-06 ***
#   s(Vort_z0_250)      0.7887      8 0.436  0.03181 *  
#   s(SST_seasonal_avg) 1.3548      7 0.418  0.09184 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 58/59
# R-sq.(adj) =   0.12   Deviance explained = 13.7%
# fREML = 5992.2  Scale est. = 66.651    n = 1684

# Compare AIC
AIC(FullModel_bathy_Rho0.013, ReducedModel_bathy_Rho0.013)

# Smooth term plots
windows(width = 10, height = 10)
par(mfrow = c(4, 5))  
plot(ReducedModel_bathy_Rho0.013, all.terms = TRUE, shade = TRUE)
# Diagnostics using better_check
windows(width = 8, height = 8)
par(mfrow = c(2, 2))
better_check(ReducedModel_bathy_Rho0.013) # using better_check instead of gam.check for tweedie
# check acf of residuals
acf(resid(ReducedModel_bathy_Rho0.013), lag.max = 15)

# Final covariates based on ReducedModel_bathy_Rho0.035
finalCovariates_bathy_factors_0.013 <- c(
  "SSS_deseasoned",
  "MXL_deseasoned",
  "Sigma22Z",
  "Sigma26Z",
  "KE_z300_600",
  "Vort_z0_250",
  "SST_seasonal_avg"
)

# PREDICT
# ------------------------------------------------------------------------------
#              Predict by Site for bathy Model - With Factors
#                      ReducedModel_bathy_Rho0.013
#                             test_data
# ------------------------------------------------------------------------------
# keeps all rows with NA and sorts them
test_data_bathy <- test_data_bathy %>%
  dplyr::select(Site, Time, DepthRange, SlopeRange, densUD, SST_seasonal_avg, all_of(finalCovariates_bathy_factors_0.013)) %>% 
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Unique sites in the test data
bathy_sites <- unique(test_data_bathy$Site)

# Initialize list to store results
site_predictions_bathy_0.013 <- list()

# Loop through each site and predict - Without Factors
for (site in bathy_sites) {
  
  cat("\nProcessing Site:", site, "\n")
  
  # Subset data by Site
  site_data_bathy_0.013 <- test_data_bathy %>% dplyr::filter(Site == site)
  
  # Predict only if data is available for the Site
  if (nrow(site_data_bathy_0.013) > 0) {
    
    # Predict using the bathy model without factors
    prediction_bathy_0.013 <- mgcv::predict.bam(
      object = ReducedModel_bathy_Rho0.013,
      newdata = site_data_bathy_0.013,
      type = 'response',
      se.fit = TRUE, # returns standard errors of prediction in addition to the predictions (FALSE only returns predictors)
      terms = NULL, # default meaning all terms from model will be included
      exclude = NULL, # default therefore no terms will be excluded
      block.size = 50000, # default size of data chunks for computation time and memory
      newdata.guaranteed = FALSE, # checks whether predictors are in correct structure - TRUE if sure that they already are
      na.action = na.pass, # keeps NA values and produces NA in the output, doesn't remove (na.omit if want them removed)
      cluster = NULL, # default is not specifying a computing cluster to run on
      discrete = TRUE, # since model was run with this, keep TRUE
      n.threads = 1, # runs parallel computations if higher than 1
      gc.level = 0 # garbage collection for memory usage - if set to 1 or 2 can clear unused memory
    )
    
    # Combine predictions with the site data
    site_data_bathy_0.013 <- site_data_bathy_0.013 %>%
      dplyr::mutate(
        predicted_response_bathy_0.013 = prediction_bathy_0.013$fit,
        predicted_se_bathy_0.013 = prediction_bathy_0.013$se.fit,
        residuals_bathy_0.013 = densUD - prediction_bathy_0.013$fit
      )
    
    # Store in the list
    site_predictions_bathy_0.013[[site]] <- site_data_bathy_0.013
  }
}

# Combine all site-specific predictions into a single data frame
predictions_bathy_by_site_no_factors <- do.call(rbind, site_predictions_bathy_0.013)

# Plot by Site - Without Factors
windows()
ggplot(predictions_bathy_by_site_no_factors, aes(x = Time)) +
  geom_line(aes(y = densUD, color = "Observed"), linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = predicted_response_bathy_0.013, color = "Predicted"), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Observed" = "grey", "Predicted" = "red")) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b %Y") +
  facet_wrap(~ Site, scales = "free_y", ncol = 2) +
  labs(
    title = "Observed vs. Predicted Density by Site (bathy Sites - With Factors) Test Data",
    x = "Time",
    y = "Density (densUD)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines")
  )


# ------------------------------------------------------------------------------
#             Predict by Site for bathy Model with factors 
#                      ReducedModel_bathy_Rho0.013
#                           Training Data
# ------------------------------------------------------------------------------
# keeps all rows with NA and sorts them
train_data_bathy <- bathy_filtered_data_with_factors %>%
  dplyr::select(Site, Time, DepthRange, SlopeRange, densUD, SST_seasonal_avg, all_of(finalCovariates_bathy_factors_0.013)) %>% 
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Unique sites in the training data
train_bathy_sites <- unique(train_data_bathy$Site)

# Initialize list to store results for training data
train_site_predictions_bathy_0.013 <- list()

# Loop through each site and predict - Training Data
for (site in train_bathy_sites) {
  
  cat("\nProcessing Site (Training Data):", site, "\n")
  
  # Subset data by Site
  site_train_data_bathy_0.013 <- train_data_bathy %>% dplyr::filter(Site == site)
  
  # Predict only if data is available for the Site
  if (nrow(site_train_data_bathy_0.013) > 0) {
    
    # Predict using the bathy model with factors
    prediction_train_bathy_0.013 <- mgcv::predict.bam(
      object = ReducedModel_bathy_Rho0.013,
      newdata = site_train_data_bathy_0.013,
      type = 'response',
      se.fit = TRUE,
      na.action = na.pass,
      discrete = TRUE
    )
    
    # Combine predictions with the site data
    site_train_data_bathy_0.013 <- site_train_data_bathy_0.013 %>%
      dplyr::mutate(
        predicted_response_train_bathy_0.013 = prediction_train_bathy_0.013$fit,
        predicted_se_train_bathy_0.013 = prediction_train_bathy_0.013$se.fit,
        residuals_train_bathy_0.013 = densUD - prediction_train_bathy_0.013$fit
      )
    
    # Store in the list
    train_site_predictions_bathy_0.013[[site]] <- site_train_data_bathy_0.013
  }
}

# Combine all site-specific predictions into a single data frame for training data
predictions_train_bathy_by_site_factors <- do.call(rbind, train_site_predictions_bathy_0.013)

# Plot by Site - Training Data with Factors
windows()
ggplot(predictions_train_bathy_by_site_factors, aes(x = Time)) +
  geom_line(aes(y = densUD, color = "Observed"), linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = predicted_response_train_bathy_0.013, color = "Predicted"), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Observed" = "grey", "Predicted" = "red")) +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  facet_wrap(~ Site, scales = "free_y") +
  labs(
    title = "Observed vs. Predicted Density by Site (Training Data - bathy Sites - With Factors)",
    x = "Time",
    y = "Density (densUD)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# ------------------------------------------------------------------------------
#                     3. BAM MODEL (TWEEDIE) for BATHY SITES
#                   baseline model without autocorrelation
#                    without FACTORS: DepthRange and SlopeRange
# ------------------------------------------------------------------------------

gc() # Cleans memory before running the model iterations, optional

# Prep bathy dataset
bathy_gam_data_wout_DepthSlope <- bathy_filtered_data_with_factors %>%
  dplyr::select(Site, Time, densUD, SST_seasonal_avg, all_of(selectedCovariates_bathy)) %>% # removed factors here
  stats::na.omit()

kVal <- 9
smooth_terms <- paste0(
  "s(", selectedCovariates_bathy[selectedCovariates_bathy != "SST_seasonal_avg"], ", bs = 'ts', k= ", kVal, ")")
# add SST seasonal avg with cyclic spline
seasonal_term <- "s(SST_seasonal_avg, bs = 'cc', k = 9)"
all_smooths <- c(smooth_terms, seasonal_term)

# construct the formula
formula_text <- paste(
  "densUD ~ ", # remove factors here too
  paste(all_smooths, collapse = "+")
)
formulaModel <- as.formula(formula_text)

# Fit bam Tweedie model for bathy sites, with site-level random effects
FullModel_bathy_woutF <- mgcv::bam(
  formula = formulaModel,
  family = mgcv::tw(),
  data = bathy_gam_data,
  method = "fREML",
  select = TRUE,
  AR.start = NULL
)

namesCovariates_bathy = selectedCovariates_bathy
iVarList = seq_along (selectedCovariates_bathy) # for numbering when going through selections
NVars = length(selectedCovariates_bathy)

# based on baseline model, determine the value of lag 1 = Rho to try in next model
r1 <- itsadug::start_value_rho(FullModel_bathy_woutF, plot=TRUE)
acf(resid(FullModel_bathy_woutF), plot=FALSE)$acf[2]
# [1] 0.1660331 = 0.166
# ------------------------------------------------------------------------------
#                     4. BAM MODEL (TWEEDIE) for BATHY SITES
#             fit with autocorrelation AR1 style using Rho + AR.start
#                             without factors
# ------------------------------------------------------------------------------
# Prep bathy dataset
bathy_gam_data_wout_DepthSlope <- bathy_filtered_data_with_factors %>%
  dplyr::select(Site, Time, densUD, SST_seasonal_avg, all_of(selectedCovariates_bathy)) %>% 
  stats::na.omit() %>%
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    AR.start = dplyr::row_number() == 1 
  ) %>%
  dplyr::ungroup()

kVal <- 9
smooth_terms <- paste0(
  "s(", selectedCovariates_bathy[selectedCovariates_bathy != "SST_seasonal_avg"], ", bs = 'ts', k= ", kVal, ")")
# add SST seasonal avg with cyclic spline
seasonal_term <- "s(SST_seasonal_avg, bs = 'cc', k = 9)"
all_smooths <- c(smooth_terms, seasonal_term)

# Model formula
formula_text <- paste(
  "densUD ~ ",
  paste(all_smooths, collapse = "+")
)
formulaModel <- as.formula(formula_text)

# Fit bam Tweedie model for bathy sites with AR1 autocorrelation
FullModel_bathy_Rho0.166 <- mgcv::bam(
  formula = formulaModel,
  family = mgcv::tw(),
  data = bathy_gam_data,
  method = "fREML",
  select = TRUE,
  discrete = TRUE,
  rho = 0.166,
  AR.start = bathy_gam_data$AR.start
)

namesCovariates_bathy = selectedCovariates_bathy
iVarList = seq_along (selectedCovariates_bathy) # for numbering when going through selections
NVars = length(selectedCovariates_bathy)

summary(FullModel_bathy_Rho0.166)
# Family: Tweedie(p=1.388) 
# Link function: log 
# 
# Formula:
#   densUD ~ s(MoonIllu, bs = "ts", k = 9) + s(SST_deseasoned, bs = "ts", 
#                                              k = 9) + s(SSS_deseasoned, bs = "ts", k = 9) + s(Chla_log, 
#                                                                                               bs = "ts", k = 9) + s(MXL_deseasoned, bs = "ts", k = 9) + 
#   s(UI, bs = "ts", k = 9) + s(Sigma22Z, bs = "ts", k = 9) + 
#   s(Sigma26Z, bs = "ts", k = 9) + s(KE_z300_600, bs = "ts", 
#                                     k = 9) + s(KE_z700_1250, bs = "ts", k = 9) + s(KE_z1500_3000, 
#                                                                                    bs = "ts", k = 9) + s(Vort_z0_250, bs = "ts", k = 9) + s(Vort_z300_600, 
#                                                                                                                                             bs = "ts", k = 9) + s(Vort_z700_1250, bs = "ts", k = 9) + 
#   s(Vort_z1500_3000, bs = "ts", k = 9) + s(SST_seasonal_avg, 
#                                            bs = "cc", k = 9)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   5.4676     0.1933   28.29   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
#   s(MoonIllu)         1.994e-05      8 0.000 0.756326    
#   s(SST_deseasoned)   8.023e-01      8 0.453 0.029412 *  
#   s(SSS_deseasoned)   8.652e-01      8 0.697 0.006820 ** 
#   s(Chla_log)         3.285e+00      8 3.420 1.99e-06 ***
#   s(MXL_deseasoned)   6.230e-01      8 0.201 0.103021    
#   s(UI)               5.141e+00      8 2.699 0.000312 ***
#   s(Sigma22Z)         3.093e-01      8 0.056 0.212972    
#   s(Sigma26Z)         3.604e+00      8 4.816  < 2e-16 ***
#   s(KE_z300_600)      3.162e+00      8 2.798 1.64e-05 ***
#   s(KE_z700_1250)     5.203e-05      8 0.000 0.482116    
#   s(KE_z1500_3000)    2.203e+00      8 0.971 0.014982 *  
#   s(Vort_z0_250)      7.397e-01      8 0.336 0.047270 *  
#   s(Vort_z300_600)    1.764e-05      8 0.000 0.574002    
#   s(Vort_z700_1250)   6.521e-01      8 0.226 0.085651 .  
#   s(Vort_z1500_3000)  2.913e-05      8 0.000 0.602520    
#   s(SST_seasonal_avg) 2.885e-05      7 0.000 0.746186    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.134   Deviance explained = 14.5%
# fREML = 6004.5  Scale est. = 67.444    n = 1684


# shrinking model 3 without factors
# Drop the first non-significant smooth term (p > 0.2)
vars_to_drop <- c("Chla_log", "MoonIllu", "KE_z700_1250", "Vort_z300_600", 
                  "Vort_z1500_3000")

pattern <- paste(vars_to_drop, collapse = "|")

all_smooths_reduced <- all_smooths[!grepl(pattern, all_smooths)]

formula_text_reduced <- paste(
  "densUD ~ ",
  paste(all_smooths_reduced, collapse = "+")
)
formulaModel_reduced <- as.formula(formula_text_reduced)

# Refit the model
ReducedModel_bathy_Rho0.166_woutF <- mgcv::bam(
  formula = formulaModel_reduced,
  family = mgcv::tw(),
  data = bathy_gam_data,
  method = "fREML",
  select = TRUE,
  discrete = TRUE,
  rho = 0.166,
  AR.start = bathy_gam_data$AR.start
)

summary(ReducedModel_bathy_Rho0.166_woutF)
# Family: Tweedie(p=1.394) 
# Link function: log 
# 
# Formula:
#   densUD ~ s(SST_deseasoned, bs = "ts", k = 9) + s(SSS_deseasoned, 
#                                                    bs = "ts", k = 9) + s(MXL_deseasoned, bs = "ts", k = 9) + 
#   s(UI, bs = "ts", k = 9) + s(Sigma22Z, bs = "ts", k = 9) + 
#   s(Sigma26Z, bs = "ts", k = 9) + s(KE_z300_600, bs = "ts", 
#                                     k = 9) + s(KE_z1500_3000, bs = "ts", k = 9) + s(Vort_z0_250, 
#                                                                                     bs = "ts", k = 9) + s(Vort_z700_1250, bs = "ts", k = 9) + 
#   s(SST_seasonal_avg, bs = "cc", k = 9)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   5.0153     0.3214   15.61   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
#   s(SST_deseasoned)   0.7592      8 0.363  0.04390 *  
#   s(SSS_deseasoned)   1.0403      8 3.001 6.40e-07 ***
#   s(MXL_deseasoned)   0.4328      8 0.095  0.18180    
#   s(UI)               4.0929      8 1.776  0.00372 ** 
#   s(Sigma22Z)         0.8530      8 0.616  0.01196 *  
#   s(Sigma26Z)         3.4970      8 4.467  < 2e-16 ***
#   s(KE_z300_600)      2.9886      8 2.381 7.74e-05 ***
#   s(KE_z1500_3000)    2.1257      8 0.866  0.02222 *  
#   s(Vort_z0_250)      0.7974      8 0.447  0.02738 *  
#   s(Vort_z700_1250)   0.7104      8 0.293  0.06184 .  
#   s(SST_seasonal_avg) 1.6737      7 0.669  0.02501 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.114   Deviance explained = 12.8%


AIC(FullModel_bathy_Rho0.166, ReducedModel_bathy_Rho0.166_woutF)

# Final covariates based on ReducedModel_bathy_Rho0.166_woutF
finalCovariates_bathy_0.166 <- c(
  "SST_deseasoned",
  "SSS_deseasoned",
  "MXL_deseasoned",
  "UI",
  "Sigma22Z",
  "Sigma26Z",
  "KE_z300_600",
  "KE_z1500_3000",
  "Vort_z0_250",
  "Vort_z700_1250",
  "SST_seasonal_avg"
)

# Smooth term plots
windows(width = 10, height = 10)
par(mfrow = c(4, 5))  
plot(ReducedModel_bathy_Rho0.166_woutF, all.terms = TRUE, shade = TRUE)
# Diagnostics using better_check
windows(width = 8, height = 8)
par(mfrow = c(2, 2))
better_check(ReducedModel_bathy_Rho0.166_woutF) # using better_check instead of gam.check for tweedie
# check acf of residuals
acf(resid(ReducedModel_bathy_Rho0.166_woutF), lag.max = 15)

# ------------------------------------------------------------------------------
#              Predict by Site for bathy Model - Without Factors
#                      ReducedModel_bathy_Rho0.166_woutF
#                             test_data
# ------------------------------------------------------------------------------
# may have to refresh bathy test data before running this:
test_data_bathy <- test_data %>%
  dplyr::filter(abs(Depth) >= 1500)

test_data_bathy <- test_data_bathy %>%
  dplyr::mutate(
    DepthRange = factor(DepthRange, levels = levels(train_data$DepthRange)),
    SlopeRange = factor(SlopeRange, levels = levels(train_data$SlopeRange))
  ) 

# ok now predicting code
# keeps all rows with NA and sorts them
test_data_bathy <- test_data_bathy %>%
  dplyr::select(Site, Time, densUD, SST_seasonal_avg, all_of(finalCovariates_bathy_0.166)) %>% 
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Unique sites in the test data
bathy_sites <- unique(test_data_bathy$Site)

# Initialize list to store results
site_predictions_bathy_0.166 <- list()

# Loop through each site and predict - Without Factors
for (site in bathy_sites) {
  
  cat("\nProcessing Site:", site, "\n")
  
  # Subset data by Site
  site_data_bathy_0.166 <- test_data_bathy %>% dplyr::filter(Site == site)
  
  # Predict only if data is available for the Site
  if (nrow(site_data_bathy_0.166) > 0) {
    
    # Predict using the bathy model without factors
    prediction_bathy_0.166 <- mgcv::predict.bam(
      object = ReducedModel_bathy_Rho0.166_woutF,
      newdata = site_data_bathy_0.166,
      type = 'response',
      se.fit = TRUE, # returns standard errors of prediction in addition to the predictions (FALSE only returns predictors)
      terms = NULL, # default meaning all terms from model will be included
      exclude = NULL, # default therefore no terms will be excluded
      block.size = 50000, # default size of data chunks for computation time and memory
      newdata.guaranteed = FALSE, # checks whether predictors are in correct structure - TRUE if sure that they already are
      na.action = na.pass, # keeps NA values and produces NA in the output, doesn't remove (na.omit if want them removed)
      cluster = NULL, # default is not specifying a computing cluster to run on
      discrete = TRUE, # since model was run with this, keep TRUE
      n.threads = 1, # runs parallel computations if higher than 1
      gc.level = 0 # garbage collection for memory usage - if set to 1 or 2 can clear unused memory
    )
    
    # Combine predictions with the site data
    site_data_bathy_0.166 <- site_data_bathy_0.166 %>%
      dplyr::mutate(
        predicted_response_bathy_0.166 = prediction_bathy_0.166$fit,
        predicted_se_bathy_0.166 = prediction_bathy_0.166$se.fit,
        residuals_bathy_0.166 = densUD - prediction_bathy_0.166$fit
      )
    
    # Store in the list
    site_predictions_bathy_0.166[[site]] <- site_data_bathy_0.166
  }
}

# Combine all site-specific predictions into a single data frame
predictions_bathy_by_site_no_factors <- do.call(rbind, site_predictions_bathy_0.166)

# Plot by Site - Without Factors
windows()
ggplot(predictions_bathy_by_site_no_factors, aes(x = Time)) +
  geom_line(aes(y = densUD, color = "Observed"), linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = predicted_response_bathy_0.166, color = "Predicted"), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Observed" = "grey", "Predicted" = "red")) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b %Y") +
  facet_wrap(~ Site, scales = "free_y", ncol = 2) +
  labs(
    title = "Observed vs. Predicted Density by Site (Testing Data - bathy Sites - Without Factors)",
    x = "Time",
    y = "Density (densUD)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines")
  )


# ------------------------------------------------------------------------------
#             Predict by Site for bathy Model without factors 
#                      ReducedModel_bathy2_Rho0.166_woutF
#                           Training Data
# ------------------------------------------------------------------------------
# keeps all rows with NA and sorts them
train_data_bathy <- bathy_filtered_data_with_factors %>%
  dplyr::select(Site, Time, densUD, SST_seasonal_avg, all_of(finalCovariates_bathy_0.166)) %>% 
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Unique sites in the training data
train_bathy_sites <- unique(train_data_bathy$Site)

# Initialize list to store results for training data
train_site_predictions_bathy_0.166 <- list()

# Loop through each site and predict - Training Data
for (site in train_bathy_sites) {
  
  cat("\nProcessing Site (Training Data):", site, "\n")
  
  # Subset data by Site
  site_train_data_bathy_0.166 <- train_data_bathy %>% dplyr::filter(Site == site)
  
  # Predict only if data is available for the Site
  if (nrow(site_train_data_bathy_0.166) > 0) {
    
    # Predict using the bathy model with factors
    prediction_train_bathy_0.166 <- mgcv::predict.bam(
      object = ReducedModel_bathy_Rho0.166_woutF,
      newdata = site_train_data_bathy_0.166,
      type = 'response',
      se.fit = TRUE,
      na.action = na.pass,
      discrete = TRUE
    )
    
    # Combine predictions with the site data
    site_train_data_bathy_0.166 <- site_train_data_bathy_0.166 %>%
      dplyr::mutate(
        predicted_response_train_bathy_0.166 = prediction_train_bathy_0.166$fit,
        predicted_se_train_bathy_0.166 = prediction_train_bathy_0.166$se.fit,
        residuals_train_bathy_0.166 = densUD - prediction_train_bathy_0.166$fit
      )
    
    # Store in the list
    train_site_predictions_bathy_0.166[[site]] <- site_train_data_bathy_0.166
  }
}

# Combine all site-specific predictions into a single data frame for training data
predictions_train_bathy_by_site_factors <- do.call(rbind, train_site_predictions_bathy_0.166)

# Plot by Site - Training Data with Factors
windows()
ggplot(predictions_train_bathy_by_site_factors, aes(x = Time)) +
  geom_line(aes(y = densUD, color = "Observed"), linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = predicted_response_train_bathy_0.166, color = "Predicted"), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Observed" = "grey", "Predicted" = "red")) +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  facet_wrap(~ Site, scales = "free_y") +
  labs(
    title = "Observed vs. Predicted Density by Site (Training Data - bathy Sites - Without Factors)",
    x = "Time",
    y = "Density (densUD)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# ------------------------------------------------------------------------------
#                     1. BAM MODEL (TWEEDIE) for MESO SITES
#                   baseline model without autocorrelation
#                    with DepthRange and SlopeRange
# ------------------------------------------------------------------------------
# Instead, running bam() because it supports tweedie with AR1 structure

gc() # Cleans memory before running the model iterations, optional

# Prep meso dataset
meso_gam_data <- meso_filtered_data_with_factors %>%
  dplyr::select(Site, Time, DepthRange, SlopeRange, densUD, SST_seasonal_avg, all_of(selectedCovariates_meso)) %>% 
  stats::na.omit()

kVal <- 9
smooth_terms <- paste0(
  "s(", selectedCovariates_meso[selectedCovariates_meso != "SST_seasonal_avg"], ", bs = 'ts', k= ", kVal, ")")
# add SST seasonal avg with cyclic spline
seasonal_term <- "s(SST_seasonal_avg, bs = 'cc', k = 9)"
all_smooths <- c(smooth_terms, seasonal_term)

# construct the formula
formula_text <- paste(
  "densUD ~ DepthRange + SlopeRange +",
  paste(all_smooths, collapse = "+")
)
formulaModel <- as.formula(formula_text)

# Fit bam Tweedie model for bathy sites, with site-level random effects
FullModel_meso <- mgcv::bam(
  formula = formulaModel,
  family = mgcv::tw(),
  data = meso_gam_data,
  method = "fREML",
  select = TRUE,
  AR.start = NULL
)

namesCovariates_meso = selectedCovariates_meso
iVarList = seq_along (selectedCovariates_meso) # for numbering when going through selections
NVars = length(selectedCovariates_meso)

# based on baseline model, determine the value of lag 1 = Rho to try in next model
r1 <- itsadug::start_value_rho(FullModel_meso, plot=TRUE)
acf(resid(FullModel_meso), plot=FALSE)$acf[2]
# [1] 0.0241708 - 0.024

# ------------------------------------------------------------------------------
#                     2. BAM MODEL (TWEEDIE) for MESO SITES
#             fit with autocorrelation AR1 style using Rho + AR.start
#                               with factors
# ------------------------------------------------------------------------------
# Prep meso dataset
meso_gam_data <- meso_filtered_data_with_factors %>%
  dplyr::select(Site, Time, DepthRange, SlopeRange, densUD, SST_seasonal_avg, all_of(selectedCovariates_meso)) %>% 
  stats::na.omit() %>%
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

kVal <- 9
smooth_terms <- paste0(
  "s(", selectedCovariates_meso[selectedCovariates_meso != "SST_seasonal_avg"], ", bs = 'ts', k= ", kVal, ")")
# add SST seasonal avg with cyclic spline
seasonal_term <- "s(SST_seasonal_avg, bs = 'cc', k = 9)"
all_smooths <- c(smooth_terms, seasonal_term)

# Model formula
formula_text <- paste(
  "densUD ~ DepthRange + SlopeRange +",
  paste(all_smooths, collapse = "+")
)
formulaModel <- as.formula(formula_text)

# Fit bam Tweedie model for meso sites with AR1 autocorrelation
FullModel_meso_Rho0.024 <- mgcv::bam(
  formula = formulaModel,
  family = mgcv::tw(),
  data = meso_gam_data,
  method = "fREML",
  select = TRUE,
  discrete = TRUE,
  rho = 0.024,
  AR.start = meso_gam_data$AR.start
)

namesCovariates_meso = selectedCovariates_meso
iVarList = seq_along (selectedCovariates_meso) # for numbering when going through selections
NVars = length(selectedCovariates_meso)

summary(FullModel_meso_Rho0.024)
# Family: Tweedie(p=1.406) 
# Link function: log 
# 
# Formula:
#   densUD ~ DepthRange + SlopeRange + s(MoonIllu, bs = "ts", k = 9) + 
#   s(SST_deseasoned, bs = "ts", k = 9) + s(SSS_deseasoned, bs = "ts", 
#                                           k = 9) + s(Chla_log, bs = "ts", k = 9) + s(MXL_deseasoned, 
#                                                                                      bs = "ts", k = 9) + s(UI, bs = "ts", k = 9) + s(Sigma22Z, 
#                                                                                                                                      bs = "ts", k = 9) + s(Sigma26Z, bs = "ts", k = 9) + s(KE_z300_600, 
#                                                                                                                                                                                            bs = "ts", k = 9) + s(KE_z700_1250, bs = "ts", k = 9) + s(Vort_z0_250, 
#                                                                                                                                                                                                                                                      bs = "ts", k = 9) + s(Vort_z300_600, bs = "ts", k = 9) + 
#   s(Vort_z700_1250, bs = "ts", k = 9) + s(SST_seasonal_avg, 
#                                           bs = "cc", k = 9)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   6.69453    0.11834  56.569  < 2e-16 ***
#   DepthRange.L  0.23693    0.06756   3.507 0.000462 ***
#   DepthRange.Q -0.76205    0.11195  -6.807 1.26e-11 ***
#   SlopeRange.L  0.70870    0.06158  11.509  < 2e-16 ***
#   SlopeRange.Q  0.95209    0.09793   9.722  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
#   s(MoonIllu)         5.677e-01      8 0.159 0.133282    
#   s(SST_deseasoned)   1.934e-04      8 0.000 0.344414    
#   s(SSS_deseasoned)   4.752e-01      8 0.111 0.162497    
#   s(Chla_log)         9.718e-01      8 1.340 0.000514 ***
#   s(MXL_deseasoned)   9.382e-01      8 1.194 0.001067 ** 
#   s(UI)               3.114e+00      8 3.242 3.40e-06 ***
#   s(Sigma22Z)         2.242e+00      8 3.359  < 2e-16 ***
#   s(Sigma26Z)         1.019e+00      8 1.817 7.46e-05 ***
#   s(KE_z300_600)      1.036e+00      8 1.759 8.93e-05 ***
#   s(KE_z700_1250)     1.220e+00      8 5.448  < 2e-16 ***
#   s(Vort_z0_250)      2.241e-05      8 0.000 0.684584    
#   s(Vort_z300_600)    8.978e-01      8 0.887 0.004410 ** 
#   s(Vort_z700_1250)   2.593e-05      8 0.000 0.566896    
#   s(SST_seasonal_avg) 2.658e+00      7 2.539 7.03e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.168   Deviance explained = 16.6%
# fREML = 8120.7  Scale est. = 43.637    n = 2389

AIC(FullModel_meso_Rho0.024)

# shrinkage for model 2
# Drop the first non-significant smooth term (p > 0.2)
vars_to_drop <- c("Chla_log", "SST_deseasoned", "SSS_deseasoned", "Vort_z0_250", "Vort_z700_1250")

pattern <- paste(vars_to_drop, collapse = "|")

all_smooths_reduced <- all_smooths[!grepl(pattern, all_smooths)]

formula_text_reduced <- paste(
  "densUD ~ DepthRange + SlopeRange +",
  paste(all_smooths_reduced, collapse = "+")
)
formulaModel_reduced <- as.formula(formula_text_reduced)

# Refit the model
ReducedModel_meso_Rho0.024 <- mgcv::bam(
  formula = formulaModel_reduced,
  family = mgcv::tw(),
  data = meso_gam_data,
  method = "fREML",
  select = TRUE,
  discrete = TRUE,
  rho = 0.024,
  AR.start = meso_gam_data$AR.start
)

summary(ReducedModel_meso_Rho0.024)
# Family: Tweedie(p=1.41) 
# Link function: log 
# 
# Formula:
#   densUD ~ DepthRange + SlopeRange + s(MoonIllu, bs = "ts", k = 9) + 
#   s(MXL_deseasoned, bs = "ts", k = 9) + s(UI, bs = "ts", k = 9) + 
#   s(Sigma22Z, bs = "ts", k = 9) + s(Sigma26Z, bs = "ts", k = 9) + 
#   s(KE_z300_600, bs = "ts", k = 9) + s(KE_z700_1250, bs = "ts", 
#                                        k = 9) + s(Vort_z300_600, bs = "ts", k = 9) + s(SST_seasonal_avg, 
#                                                                                        bs = "cc", k = 9)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   6.72143    0.11842  56.758  < 2e-16 ***
#   DepthRange.L  0.28522    0.06680   4.270 2.03e-05 ***
#   DepthRange.Q -0.78765    0.11078  -7.110 1.53e-12 ***
#   SlopeRange.L  0.68666    0.06124  11.212  < 2e-16 ***
#   SlopeRange.Q  0.92974    0.09765   9.522  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
#   s(MoonIllu)         1.9130      8 0.929 0.014424 *  
#   s(MXL_deseasoned)   0.9452      8 1.281 0.000786 ***
#   s(UI)               3.1633      8 3.420 1.82e-06 ***
#   s(Sigma22Z)         2.2056      8 6.557  < 2e-16 ***
#   s(Sigma26Z)         1.0649      8 2.646 2.40e-06 ***
#   s(KE_z300_600)      1.0629      8 1.961 3.79e-05 ***
#   s(KE_z700_1250)     1.2744      8 5.719  < 2e-16 ***
#   s(Vort_z300_600)    0.8756      8 0.750 0.008142 ** 
#   s(SST_seasonal_avg) 2.4692      7 1.740 0.001379 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.156   Deviance explained = 16.4%
# fREML = 8120.6  Scale est. = 43.549    n = 2389

AIC(FullModel_meso_Rho0.024, ReducedModel_meso_Rho0.024)

# Smooth term plots
windows(width = 10, height = 10)
par(mfrow = c(4, 4))  
plot(ReducedModel_meso_Rho0.024, all.terms = TRUE, shade = TRUE)
# Diagnostics using better_check
windows(width = 8, height = 8)
par(mfrow = c(2, 2))
better_check(ReducedModel_meso_Rho0.024) # using better_check instead of gam.check for tweedie
# check acf of residuals
acf(resid(ReducedModel_meso_Rho0.024), lag.max = 30)

# Final covariates based on ReducedModel_meso_Rho0.024
finalCovariates_meso_0.024 <- c(
  "MoonIllu",
  "MXL_deseasoned",
  "UI",
  "Sigma22Z",
  "Sigma26Z",
  "KE_z300_600",
  "KE_z700_1250",
  "Vort_z300_600",
  "SST_seasonal_avg"
)

# ------------------------------------------------------------------------------
#              Predict by Site for Meso Model - With Factors
#                      ReducedModel_meso_Rho0.024
#                               test_data
# ------------------------------------------------------------------------------

# keeps all rows with NA and sorts them
test_data_meso <- test_data_meso %>%
  dplyr::select(Site, Time, DepthRange, SlopeRange, densUD, SST_seasonal_avg, all_of(finalCovariates_meso_0.024)) %>% 
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Unique sites in the test data
meso_sites <- unique(test_data_meso$Site)

# Initialize list to store results
site_predictions_meso_0.024 <- list()

# Loop through each site and predict - Without Factors
for (site in meso_sites) {
  
  cat("\nProcessing Site:", site, "\n")
  
  # Subset data by Site
  site_data_meso_0.024 <- test_data_meso %>% dplyr::filter(Site == site)
  
  # Predict only if data is available for the Site
  if (nrow(site_data_meso_0.024) > 0) {
    
    # Predict using the meso model without factors
    prediction_meso_0.024 <- mgcv::predict.bam(
      object = ReducedModel_meso_Rho0.024,
      newdata = site_data_meso_0.024,
      type = 'response',
      se.fit = TRUE,
      na.action = na.pass,
      discrete = TRUE
    )
    
    # Combine predictions with the site data
    site_data_meso_0.024 <- site_data_meso_0.024 %>%
      dplyr::mutate(
        predicted_response_meso_0.024 = prediction_meso_0.024$fit,
        predicted_se_meso_0.024 = prediction_meso_0.024$se.fit,
        residuals_meso_0.024 = densUD - prediction_meso_0.024$fit
      )
    
    # Store in the list
    site_predictions_meso_0.024[[site]] <- site_data_meso_0.024
  }
}

# Combine all site-specific predictions into a single data frame
predictions_meso_by_site_no_factors <- do.call(rbind, site_predictions_meso_0.024)

# Plot by Site - Without Factors
windows()
ggplot(predictions_meso_by_site_no_factors, aes(x = Time)) +
  geom_line(aes(y = densUD, color = "Observed"), linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = predicted_response_meso_0.024, color = "Predicted"), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Observed" = "grey", "Predicted" = "red")) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b %Y") +
  facet_wrap(~ Site, scales = "free_y", ncol = 2) +
  labs(
    title = "Observed vs. Predicted Density by Site (Meso Model - Test Data - With Factors)",
    x = "Time",
    y = "Density (densUD)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines")
  )

# ------------------------------------------------------------------------------
#              Predict by Site for Meso Model - with factors
#                      ReducedModel_meso_Rho0.024
#                         Training Data
# ------------------------------------------------------------------------------
# keeps all rows with NA and sorts them
train_data_meso <- meso_filtered_data_with_factors %>%
  dplyr::select(Site, Time, DepthRange, SlopeRange, densUD, SST_seasonal_avg, all_of(finalCovariates_meso_0.024)) %>% 
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Unique sites in the training data
train_meso_sites <- unique(train_data_meso$Site)

# Initialize list to store results for training data
train_site_predictions_meso_0.024 <- list()

# Loop through each site and predict - Training Data
for (site in train_meso_sites) {
  
  cat("\nProcessing Site (Training Data):", site, "\n")
  
  # Subset data by Site
  site_train_data_meso_0.024 <- train_data_meso %>% dplyr::filter(Site == site)
  
  # Predict only if data is available for the Site
  if (nrow(site_train_data_meso_0.024) > 0) {
    
    # Predict using the meso model with factors
    prediction_train_meso_0.024 <- mgcv::predict.bam(
      object = ReducedModel_meso_Rho0.024,
      newdata = site_train_data_meso_0.024,
      type = 'response',
      se.fit = TRUE,
      na.action = na.pass,
      discrete = TRUE
    )
    
    # Combine predictions with the site data
    site_train_data_meso_0.024 <- site_train_data_meso_0.024 %>%
      dplyr::mutate(
        predicted_response_train_meso_0.024 = prediction_train_meso_0.024$fit,
        predicted_se_train_meso_0.024 = prediction_train_meso_0.024$se.fit,
        residuals_train_meso_0.024 = densUD - prediction_train_meso_0.024$fit
      )
    
    # Store in the list
    train_site_predictions_meso_0.024[[site]] <- site_train_data_meso_0.024
  }
}

# Combine all site-specific predictions into a single data frame for training data
predictions_train_meso_by_site_factors <- do.call(rbind, train_site_predictions_meso_0.024)

# Plot by Site - Training Data with Factors
windows()
ggplot(predictions_train_meso_by_site_factors, aes(x = Time)) +
  geom_line(aes(y = densUD, color = "Observed"), linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = predicted_response_train_meso_0.024, color = "Predicted"), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Observed" = "grey", "Predicted" = "red")) +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  facet_wrap(~ Site, scales = "free_y") +
  labs(
    title = "Observed vs. Predicted Density by Site (Training Data - Meso Sites - With Factors)",
    x = "Time",
    y = "Density (densUD)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ------------------------------------------------------------------------------
#                     3. BAM MODEL (TWEEDIE) for MESO SITES
#                   baseline model without autocorrelation
#                    without FACTORS: DepthRange and SlopeRange
# ------------------------------------------------------------------------------

gc() # Cleans memory before running the model iterations, optional

# Prep meso dataset
meso_gam_data_wout_DepthSlope <- meso_filtered_data_with_factors %>%
  dplyr::select(Site, Time, densUD, SST_seasonal_avg, all_of(selectedCovariates_meso)) %>% # removed factors here
  stats::na.omit()

kVal <- 9
smooth_terms <- paste0(
  "s(", selectedCovariates_meso[selectedCovariates_meso != "SST_seasonal_avg"], ", bs = 'ts', k= ", kVal, ")")
# add SST seasonal avg with cyclic spline
seasonal_term <- "s(SST_seasonal_avg, bs = 'cc', k = 9)"
all_smooths <- c(smooth_terms, seasonal_term)

# construct the formula
formula_text <- paste(
  "densUD ~ ", # remove factors here too
  paste(all_smooths, collapse = "+")
)
formulaModel <- as.formula(formula_text)

# Fit bam Tweedie model for meso sites, with site-level random effects
FullModel_meso_woutF <- mgcv::bam(
  formula = formulaModel,
  family = mgcv::tw(),
  data = meso_gam_data_wout_DepthSlope,
  method = "fREML",
  select = TRUE,
  AR.start = NULL
)

namesCovariates_meso = selectedCovariates_meso
iVarList = seq_along (selectedCovariates_meso) # for numbering when going through selections
NVars = length(selectedCovariates_meso)

# based on baseline model, determine the value of lag 1 = Rho to try in next model
r1 <- itsadug::start_value_rho(FullModel_meso_woutF, plot=TRUE)
acf(resid(FullModel_meso_woutF), plot=FALSE)$acf[2]
# [1] 0.02539918 = 0.025

# ------------------------------------------------------------------------------
#                     4. BAM MODEL (TWEEDIE) for meso SITES
#             fit with autocorrelation AR1 style using Rho + AR.start
#                             without factors
# ------------------------------------------------------------------------------
# Prep meso dataset
meso_gam_data_wout_DepthSlope <- meso_filtered_data_with_factors %>%
  dplyr::select(Site, Time, densUD, SST_seasonal_avg, all_of(selectedCovariates_meso)) %>% 
  stats::na.omit() %>%
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup() 

kVal <- 9
smooth_terms <- paste0(
  "s(", selectedCovariates_meso[selectedCovariates_meso != "SST_seasonal_avg"], ", bs = 'ts', k= ", kVal, ")")
# add SST seasonal avg with cyclic spline
seasonal_term <- "s(SST_seasonal_avg, bs = 'cc', k = 9)"
all_smooths <- c(smooth_terms, seasonal_term)

# Model formula
formula_text <- paste(
  "densUD ~ ",
  paste(all_smooths, collapse = "+")
)
formulaModel <- as.formula(formula_text)

# Fit bam Tweedie model for meso sites with AR1 autocorrelation
FullModel_meso_Rho0.025 <- mgcv::bam(
  formula = formulaModel,
  family = mgcv::tw(),
  data = meso_gam_data_wout_DepthSlope,
  method = "fREML",
  select = TRUE,
  discrete = TRUE,
  rho = 0.025,
  AR.start = meso_gam_data_wout_DepthSlope$AR.start
)

namesCovariates_meso = selectedCovariates_meso
iVarList = seq_along (selectedCovariates_meso) # for numbering when going through selections
NVars = length(selectedCovariates_meso)

summary(FullModel_meso_Rho0.025)
# Family: Tweedie(p=1.47) 
# Link function: log 
# 
# Formula:
#   densUD ~ s(MoonIllu, bs = "ts", k = 9) + s(SST_deseasoned, bs = "ts", 
#                                              k = 9) + s(SSS_deseasoned, bs = "ts", k = 9) + s(Chla_log, 
#                                                                                               bs = "ts", k = 9) + s(MXL_deseasoned, bs = "ts", k = 9) + 
#   s(UI, bs = "ts", k = 9) + s(Sigma22Z, bs = "ts", k = 9) + 
#   s(Sigma26Z, bs = "ts", k = 9) + s(KE_z300_600, bs = "ts", 
#                                     k = 9) + s(KE_z700_1250, bs = "ts", k = 9) + s(Vort_z0_250, 
#                                                                                    bs = "ts", k = 9) + s(Vort_z300_600, bs = "ts", k = 9) + 
#   s(Vort_z700_1250, bs = "ts", k = 9) + s(SST_seasonal_avg, 
#                                           bs = "cc", k = 9)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   6.1812     0.1196    51.7   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
#   s(MoonIllu)         3.224e-01      8  0.059 0.224301    
#   s(SST_deseasoned)   5.958e-01      8  0.177 0.116057    
#   s(SSS_deseasoned)   4.453e-01      8  0.098 0.168381    
#   s(Chla_log)         6.088e+00      8  4.803  < 2e-16 ***
#   s(MXL_deseasoned)   2.128e-05      8  0.000 0.904740    
#   s(UI)               4.191e+00      8  5.605  < 2e-16 ***
#   s(Sigma22Z)         9.661e-01      8  1.303 0.000124 ***
#   s(Sigma26Z)         5.946e+00      8 10.709  < 2e-16 ***
#   s(KE_z300_600)      9.091e-01      8  0.831 0.005092 ** 
#   s(KE_z700_1250)     1.843e+00      8  2.419 1.05e-05 ***
#   s(Vort_z0_250)      4.561e-01      8  0.103 0.173377    
#   s(Vort_z300_600)    2.667e+00      8  2.630 5.35e-06 ***
#   s(Vort_z700_1250)   1.046e+00      8  2.043 2.67e-05 ***
#   s(SST_seasonal_avg) 2.517e+00      7  1.932 0.000609 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.136   Deviance explained = 13.3%
# fREML = 8132.5  Scale est. = 42.892    n = 2389

# shrinkage for model 4
# Drop the first non-significant smooth term (p > 0.2)
vars_to_drop <- c("Chla_log", "MoonIllu", "SST_deseasoned", "MXL_deseasoned")

pattern <- paste(vars_to_drop, collapse = "|")

all_smooths_reduced <- all_smooths[!grepl(pattern, all_smooths)]

formula_text_reduced <- paste(
  "densUD ~ ",
  paste(all_smooths_reduced, collapse = "+")
)
formulaModel_reduced <- as.formula(formula_text_reduced)

# Refit the model
ReducedModel_meso_Rho0.025 <- mgcv::bam(
  formula = formulaModel_reduced,
  family = mgcv::tw(),
  data = meso_gam_data_wout_DepthSlope,
  method = "fREML",
  select = TRUE,
  discrete = TRUE,
  rho = 0.025,
  AR.start = meso_gam_data_wout_DepthSlope$AR.start
)

summary(ReducedModel_meso_Rho0.025)
# Family: Tweedie(p=1.421) 
# Link function: log 
# 
# Formula:
#   densUD ~ s(SSS_deseasoned, bs = "ts", k = 9) + s(UI, bs = "ts", 
#                                                    k = 9) + s(Sigma22Z, bs = "ts", k = 9) + s(Sigma26Z, bs = "ts", 
#                                                                                               k = 9) + s(KE_z300_600, bs = "ts", k = 9) + s(KE_z700_1250, 
#                                                                                                                                             bs = "ts", k = 9) + s(Vort_z0_250, bs = "ts", k = 9) + s(Vort_z300_600, 
#                                                                                                                                                                                                      bs = "ts", k = 9) + s(Vort_z700_1250, bs = "ts", k = 9) + 
#   s(SST_seasonal_avg, bs = "cc", k = 9)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   6.1615     0.1162   53.02   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
# s(SSS_deseasoned)   0.9144      8  0.781 0.004698 ** 
#   s(UI)               4.0915      8  5.407  < 2e-16 ***
#   s(Sigma22Z)         0.8487      8  0.585 0.011783 *  
#   s(Sigma26Z)         6.0352      8 10.498  < 2e-16 ***
#   s(KE_z300_600)      0.8793      8  0.710 0.008752 ** 
#   s(KE_z700_1250)     2.1519      8  1.581 0.000610 ***
#   s(Vort_z0_250)      4.6999      8  2.360 0.000516 ***
#   s(Vort_z300_600)    0.7702      8  0.364 0.047796 *  
#   s(Vort_z700_1250)   3.2398      8  3.343 1.48e-06 ***
#   s(SST_seasonal_avg) 2.7729      7  2.865 2.53e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.117   Deviance explained = 12.6%
# fREML = 8124.6  Scale est. = 42.716    n = 2389

AIC(FullModel_meso_Rho0.025, ReducedModel_meso_Rho0.025)
# df      AIC
# FullModel_meso_Rho0.025    42.12473 8963.717
# ReducedModel_meso_Rho0.025 36.09505 8972.545

# Smooth term plots
windows(width = 10, height = 10)
par(mfrow = c(3, 4))  
plot(ReducedModel_meso_Rho0.025, all.terms = TRUE, shade = TRUE)
# Diagnostics using better_check
windows(width = 8, height = 8)
par(mfrow = c(2, 2))
better_check(ReducedModel_meso_Rho0.025) # using better_check instead of gam.check for tweedie
# check acf of residuals
windows()
acf(resid(ReducedModel_meso_Rho0.025), lag.max = 50)

# Final covariates based on ReducedModel_bathy_Rho0.025_woutF
finalCovariates_meso_0.025 <- c( 
  "SSS_deseasoned",
  "UI", 
  "Sigma22Z",
  "Sigma26Z",
  "KE_z300_600",
  "KE_z700_1250", 
  "Vort_z0_250",
  "Vort_z300_600",
  "Vort_z700_1250",
  "SST_seasonal_avg" 
)

# ------------------------------------------------------------------------------
#              Predict by Site for Meso Model - Without Factors
#                       ReducedModel_meso_Rho0.025
#                           Test_data
# ------------------------------------------------------------------------------
# Meso test data (Depth < 1500m) # running this again to refresh from last
test_data_meso <- test_data %>%
  dplyr::filter(abs(Depth) < 1500) %>%
  dplyr::select(-dplyr::any_of(c("KE_z1500_3000", "Vort_z1500_3000")))


# keeps all rows with NA and sorts them
test_data_meso <- test_data_meso %>%
  dplyr::select(Site, Time, densUD, SST_seasonal_avg, all_of(finalCovariates_meso_0.025)) %>% 
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Unique sites in the test data
meso_sites <- unique(test_data_meso$Site)

# Initialize list to store results
site_predictions_meso_0.025 <- list()

# Loop through each site and predict - Without Factors
for (site in meso_sites) {
  
  cat("\nProcessing Site:", site, "\n")
  
  # Subset data by Site
  site_data_meso_0.025 <- test_data_meso %>% dplyr::filter(Site == site)
  
  # Predict only if data is available for the Site
  if (nrow(site_data_meso_0.025) > 0) {
    
    # Predict using the meso model without factors
    prediction_meso_0.025 <- mgcv::predict.bam(
      object = ReducedModel_meso_Rho0.025,
      newdata = site_data_meso_0.025,
      type = 'response',
      se.fit = TRUE,
      na.action = na.pass,
      discrete = TRUE
    )
    
    # Combine predictions with the site data
    site_data_meso_0.025 <- site_data_meso_0.025 %>%
      dplyr::mutate(
        predicted_response_meso_0.025 = prediction_meso_0.025$fit,
        predicted_se_meso_0.025 = prediction_meso_0.025$se.fit,
        residuals_meso_0.025 = densUD - prediction_meso_0.025$fit
      )
    
    # Store in the list
    site_predictions_meso_0.025[[site]] <- site_data_meso_0.025
  }
}

# Combine all site-specific predictions into a single data frame
predictions_meso_by_site_no_factors <- do.call(rbind, site_predictions_meso_0.025)

# Plot by Site - Without Factors
windows()
ggplot(predictions_meso_by_site_no_factors, aes(x = Time)) +
  geom_line(aes(y = densUD, color = "Observed"), linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = predicted_response_meso_0.025, color = "Predicted"), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Observed" = "grey", "Predicted" = "red")) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b %Y") +
  facet_wrap(~ Site, scales = "free_y", ncol = 2) +
  labs(
    title = "Observed vs. Predicted Density by Site (Test Data Meso Model - Without Factors)",
    x = "Time",
    y = "Density (densUD)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines")
  )

# ------------------------------------------------------------------------------
#              Predict by Site for Meso Model - Without Factors
#                      ReducedModel_meso_Rho0.025
#                             Training Data
# ------------------------------------------------------------------------------
# keeps all rows with NA and sorts them
train_data_meso <- meso_filtered_data_with_factors %>%
  dplyr::select(Site, Time, densUD, SST_seasonal_avg, all_of(finalCovariates_meso_0.025)) %>% 
  dplyr::arrange(Site, Time) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Unique sites in the training data
train_meso_sites <- unique(train_data_meso$Site)

# Initialize list to store results for training data
train_site_predictions_meso_0.025 <- list()

# Loop through each site and predict - Training Data without Factors
for (site in train_meso_sites) {
  
  cat("\nProcessing Site (Training Data):", site, "\n")
  
  # Subset data by Site
  site_train_data_meso_0.025 <- train_data_meso %>% dplyr::filter(Site == site)
  
  # Predict only if data is available for the Site
  if (nrow(site_train_data_meso_0.025) > 0) {
    
    # Predict using the meso model without factors
    prediction_train_meso_0.025 <- mgcv::predict.bam(
      object = ReducedModel_meso_Rho0.025,
      newdata = site_train_data_meso_0.025,
      type = 'response',
      se.fit = TRUE,
      na.action = na.pass,
      discrete = TRUE
    )
    
    # Combine predictions with the site data and explicitly set to NA where predictors are NA
    site_train_data_meso_0.025 <- site_train_data_meso_0.025 %>%
      dplyr::mutate(
        predicted_response_train_meso_0.025 = prediction_train_meso_0.025$fit,
        predicted_se_train_meso_0.025 = prediction_train_meso_0.025$se.fit,
        residuals_train_meso_0.025 = densUD - prediction_train_meso_0.025$fit
      )
    
    # Store in the list
    train_site_predictions_meso_0.025[[site]] <- site_train_data_meso_0.025
  }
}

# Combine all site-specific predictions into a single data frame for training data
predictions_train_meso_by_site_no_factors <- do.call(rbind, train_site_predictions_meso_0.025)

# Plot by Site - Training Data Without Factors
windows()
ggplot(predictions_train_meso_by_site_no_factors, aes(x = Time)) +
  geom_line(aes(y = densUD, color = "Observed"), linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = predicted_response_train_meso_0.025, color = "Predicted"), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Observed" = "grey", "Predicted" = "red")) +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  facet_wrap(~ Site, scales = "free_y", ncol = 2) +
  labs(
    title = "Observed vs. Predicted Density by Site (Training Data - Meso Sites - Without Factors)",
    x = "Time",
    y = "Density (densUD)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



# ------------------------------------------------------------------------------
#         CHECK MULTICOLLINEARITY FOR ONE SITE AT A TIME
# ------------------------------------------------------------------------------

# Function to compute and plot correlation matrix
.compute_correlation <- function(data, site_name) {
  data <- data %>% dplyr::select(-dplyr::any_of(c("densUD", "RowID")))  # Remove response + ID
  spearman_corr <- Hmisc::rcorr(as.matrix(data), type = "spearman")
  
  # Plot
  plot.new()
  dev.new(width = 8, height = 8, noRStudioGD = TRUE)
  corrplot::corrplot(spearman_corr$r,
                     method = "ellipse",
                     order = "original",
                     addgrid.col = "darkgray",
                     col = colorRampPalette(c("blue", "cyan", "gray", "white", "gray", "yellow", "red"))(10),
                     tl.cex = 0.8,
                     cl.pos = "r",
                     tl.col = "black",
                     tl.srt = 90,
                     p.mat = spearman_corr$P,
                     sig.level = 0.05,
                     insig = "pch",
                     pch.col = "black",
                     pch.cex = 1,
                     diag = FALSE,
                     type = "upper",
                     addCoef.col = "black",
                     number.cex = 0.5)
  title(main = paste("Correlation Plot for", site_name, "Site"))
  
  # Save PDF
  plot_file <- file.path(species_results_path, paste(species, "Correlation_Plot_for", site_name, "Site.pdf", sep = "_"))
  dev.print(pdf, plot_file)
  dev.off()
}

# Loop through all unique sites and generate correlation plots
for (site_name in unique(dataSp_New_UD$Site)) {
  site_data <- dataSp_New_UD %>%
    dplyr::filter(Site == site_name) %>%
    dplyr::select(all_of(c("densUD", "RowID", selected.predictors))) %>%
    dplyr::select_if(is.numeric) %>%
    tidyr::drop_na()
  
  if (nrow(site_data) >= 10) {  # Only plot if enough data
    message("Generating correlation plot for site: ", site_name)
    .compute_correlation(site_data, site_name)
  } else {
    message("Skipping site ", site_name, ": not enough complete cases.")
  }
}

