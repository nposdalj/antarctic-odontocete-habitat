# =============================================================================
# config.R — Shared configuration for all manuscript scripts
# -----------------------------------------------------------------------------
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Author  : Natalie Posdaljian
# Created : 2026-11-16
#
# Purpose : Single source of truth for site metadata, spatial boxes, recording
#           periods, species definitions, data-product file names, and paths.
#           Every script in final_scripts_odontocetes/ starts with:
#
#               BASE_DIR <- "."   # path to the repository root
#               source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
#
#           Changing a value here changes it everywhere, so manuscript numbers
#           and code cannot drift apart.
#
# Notes   : All times are UTC. Longitudes are in -180..180 (degrees East).
# =============================================================================

if (!exists("BASE_DIR")) BASE_DIR <- "."

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
FINAL_DIR <- file.path(BASE_DIR, "final_scripts_odontocetes")
OUT_DIR   <- file.path(FINAL_DIR, "outputs")          # all script outputs
RAW_DIR   <- file.path(BASE_DIR, "data")              # detection logs, deployment summary
ENV_DIR   <- file.path(BASE_DIR, "Environmental Data") # existing environmental inputs

# Large downloaded environmental files (Copernicus, AVISO, HYCOM, MODIS).
# Point this at an external drive if needed; scripts create subfolders.
DOWNLOAD_DIR <- file.path(OUT_DIR, "downloads")

out_path <- function(section, file) {
  d <- file.path(OUT_DIR, section)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  file.path(d, file)
}

# -----------------------------------------------------------------------------
# 2. Sites
# -----------------------------------------------------------------------------
# Coordinates, depths, and recording times from data/DataSummary.xlsx.
# NOTE: DataSummary.xlsx lists the EI Data_End_Date as 2014-07-14, but its
# Rec_Dur (134.5 d) and the detection log (detections until 2014-07-17 10:54)
# both indicate recording ended 2014-07-17 12:03:50, which is used here.
#
# Duty cycle: recording on `duty_on_min` minutes of every `duty_cycle_min`.
# CI recorded continuously.
#
# Spatial boxes:
#   box_*   : MAIN analysis box, +/- BOX_HALF_KM km around each HARP (computed
#             below). +/- 20 km (≈ 40 × 40 km) was chosen to match odontocete
#             detection ranges and the original 40 × 40 km study design.
#   box40_* : +/- 40 km box (≈ ±0.361° latitude, ≈ ±0.74° longitude at 61°S),
#             the box used for the original AVISO/Copernicus extractions. Used
#             (i) as the DOWNLOAD region, which contains the main box, and
#             (ii) for the box-size comparison reported in the supplement.
BOX_HALF_KM            <- 20
BOX_HALF_KM_COMPARISON <- 40
SITES <- data.frame(
  site           = c("EI", "KGI", "CI"),
  site_name      = c("Elephant Island", "King George Island", "Clarence Island"),
  lat            = c(-60.88690, -61.457817, -61.251867),
  lon            = c(-55.95400, -57.941917, -53.483433),
  depth_m        = c(762, 768, 1033),
  rec_start      = as.POSIXct(c("2014-03-05 00:00:00", "2015-02-10 00:00:00",
                                "2016-02-03 23:59:59"), tz = "UTC"),
  rec_end        = as.POSIXct(c("2014-07-17 12:03:50", "2016-01-29 13:15:25",
                                "2016-12-02 05:15:05"), tz = "UTC"),
  duty_on_min    = c(5, 10, NA),
  duty_cycle_min = c(6, 13.33333, NA),
  detection_log  = c("Antarc_EI_01_Odontocetes.csv",
                     "Antarc_SSI_01_Odontocetes.csv",   # original file name for KGI
                     "Antarc_EIE_01_Odontocetes.csv"),  # original file name for CI
  box40_lat_min  = c(-61.24778, -61.81870, -61.61275),
  box40_lat_max  = c(-60.52602, -61.09694, -60.89099),
  box40_lon_min  = c(-56.69255, -58.69396, -54.23054),
  box40_lon_max  = c(-55.21545, -57.18988, -52.73632),
  stringsAsFactors = FALSE
)

# Box of half-width `half_km` (km) around each site: latitude half-width =
# half_km / 111.32 km per degree; longitude half-width scaled by cos(latitude).
site_boxes <- function(half_km) {
  dlat <- half_km / 111.32
  dlon <- half_km / (111.32 * cos(SITES$lat * pi / 180))
  data.frame(site = SITES$site,
             lat_min = SITES$lat - dlat, lat_max = SITES$lat + dlat,
             lon_min = SITES$lon - dlon, lon_max = SITES$lon + dlon)
}
main_box <- site_boxes(BOX_HALF_KM)
SITES$box_lat_min <- main_box$lat_min
SITES$box_lat_max <- main_box$lat_max
SITES$box_lon_min <- main_box$lon_min
SITES$box_lon_max <- main_box$lon_max
rm(main_box)
SITE_ORDER  <- c("EI", "KGI", "CI")
SITE_COLORS <- c(EI = "#D81B60", KGI = "#7B1FA2", CI = "#1E88E5")

# Years of environmental data around each deployment used for the monthly
# climatology (deseasoning): deployment start year - 5 to deployment end year + 5,
# truncated to the availability of each product.
CLIM_HALF_WINDOW_YEARS <- 5

# -----------------------------------------------------------------------------
# 3. Species
# -----------------------------------------------------------------------------
# depths_m: nominal depths (m) at which water-column predictors are extracted
# for each species (surface + typical foraging/dive depths). Extraction scripts
# use the nearest Copernicus model level to each nominal depth.
SPECIES <- list(
  Pm   = list(name = "Sperm whale",                 depths_m = c(0, 375, 1665), model = TRUE),
  Gm   = list(name = "Long-finned pilot whale",     depths_m = c(0, 16, 635),   model = TRUE),
  Oo   = list(name = "Killer whale",                depths_m = c(0, 11, 455),   model = TRUE),
  BW29 = list(name = "Southern bottlenose whale",   depths_m = c(0, 768),       model = TRUE),
  BW37 = list(name = "BW37 beaked whale",           depths_m = c(0, 67, 920),   model = TRUE),
  BW58 = list(name = "BW58 beaked whale",           depths_m = numeric(0),      model = FALSE)  # too rare to model (5 events)
)
# BW37 and BW58 are echolocation signal types of unconfirmed species origin,
# hypothesized to be produced by Gray's and strap-toothed beaked whales,
# respectively (Trickey et al., 2015; Baumann-Pickering et al., 2015;
# Barlow et al., 2021); Shepherd's beaked whale has also been proposed
# (Barlow et al., 2021). Labels therefore name the signal type, not a species.
SPECIES_CODES <- names(SPECIES)
MODELED_SPECIES <- names(Filter(function(s) s$model, SPECIES))
ALL_DEPTHS_M <- sort(unique(unlist(lapply(SPECIES, `[[`, "depths_m"))))

# Only echolocation clicks define presence (sperm whale slow clicks and
# killer whale HFM whistles are excluded).
PRESENCE_CALLS <- c("Clicks", "regular-click")

# Presence is resolved in 5-min bins; daily presence = minutes in bins with
# >= 1 detection.
PRESENCE_BIN_MIN <- 5

# -----------------------------------------------------------------------------
# 4. Environmental products
# -----------------------------------------------------------------------------
COPERNICUS <- list(
  phy = list(product = "GLOBAL_MULTIYEAR_PHY_001_030",
             layer   = "cmems_mod_glo_phy_my_0.083deg_P1D-m",       # 1/12°, daily
             vars    = c("thetao", "so", "uo", "vo", "zos", "mlotst", "siconc", "sithick")),
  bgc = list(product = "GLOBAL_MULTIYEAR_BGC_001_029",
             layer   = "cmems_mod_glo_bgc_my_0.25deg_P1D-m",        # 1/4°, daily
             vars    = c("chl", "o2", "nppv"))
)
# Maximum depth (m) at which chl-a and primary productivity are used as
# predictors (below the euphotic zone they are ~0 and uninformative).
CHL_PROD_MAX_DEPTH_M <- 50

# Lagged bottom-up predictors: mean over the 30-day window ending
# lag * 30 days before each date.
LAG_MONTHS      <- 1:6
LAG_WINDOW_DAYS <- 30
LAG_VARS        <- c("chla", "productivity")

# -----------------------------------------------------------------------------
# 5. Data products passed between scripts (the "data contract")
# -----------------------------------------------------------------------------
# 04_acoustic_detection
F_PRESENCE_DAILY   <- out_path("04_acoustic_detection", "presence_daily.csv")
#   date, site, rec_minutes, <species>  (presence minutes per day, clicks only)
F_PRESENCE_5MIN    <- out_path("04_acoustic_detection", "presence_5min.csv")
#   bin_start, site, recorded (logical), <species> (logical)
F_DUTY_SIM         <- out_path("04_acoustic_detection", "duty_cycle_simulation.csv")

# 05_environmental_data — daily box statistics at every model level
#   (span = whatever the input NetCDFs cover; climatology window once the
#    10-yr downloads exist, otherwise deployment periods only)
F_COP_PHY_DEPTH    <- out_path("05_environmental_data", "copernicus_phy_by_depth.csv")
#   date, site, depth_m (model level), temperature, temperature_sd, salinity,
#   salinity_sd, u, v, u_mad, v_mad, coverage
F_COP_BGC_DEPTH    <- out_path("05_environmental_data", "copernicus_bgc_by_depth.csv")
#   date, site, depth_m (model level), chla, chla_sd, o2, o2_sd,
#   productivity, productivity_sd, coverage
F_COP_SURFACE      <- out_path("05_environmental_data", "copernicus_surface.csv")
#   date, site, ssh, ssh_sd, mixed_layer, mixed_layer_sd, ice_conc,
#   ice_conc_sd, ice_thickness   (ice NA -> 0: NaN = open water in source)
F_COP_COVERAGE     <- out_path("05_environmental_data", "copernicus_depth_coverage.csv")
#   site, product, depth_m, coverage (fraction of box cells with water)
F_FSLE_DAILY       <- out_path("05_environmental_data", "fsle_daily.csv")
#   date, site, FSLE, FSLE_sd, fsle_orient, fsle_orient_sd
F_ENV_DAILY        <- out_path("05_environmental_data", "env_daily_wide.csv")
#   date, site, julian_day, ice_conc, ice_thickness, SSH, mixed_layer, FSLE,
#   fsle_orient, temperature_<d>, salinity_<d>, o2_<d>, EKE_<d>, EKE_mad_<d>,
#   chla_<d>, productivity_<d> (d <= CHL_PROD_MAX_DEPTH_M only),
#   chla_lag<N>, productivity_lag<N> (surface, 30-day window means);
#   <d> = nominal depth from ALL_DEPTHS_M (nearest model level)
F_MODEL_INPUT      <- out_path("05_environmental_data", "model_input_daily.csv")
#   recording days only: env_daily_wide joined to presence_daily (+ rec_minutes)

# Helper: site row by code
site_info <- function(code) SITES[SITES$site == code, ]
