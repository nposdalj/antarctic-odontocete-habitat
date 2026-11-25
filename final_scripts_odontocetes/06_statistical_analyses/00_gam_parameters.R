# =============================================================================
# Script  : 00_gam_parameters.R
# Section : 2.6 Statistical analyses (also used by 3.2 GAM results)
# Purpose : Single parameter block for the GAM habitat-model workflow
#           (Stages 1-3 in 06_statistical_analyses/ and the result scripts in
#           08_results_gams/). Every analysis setting that appears in the
#           methods text is defined here, once.
# Inputs  : final_scripts_odontocetes/config.R (sourced by the calling script)
# Outputs : none (defines objects only)
# Author  : Natalie Posdaljian
# Date    : 2026-11-16
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Sourced by 01_stage1_automated_selection.R,
#           02_stage2_decision_table.R, 03_stage3_harmonized_models.R,
#           04_sea_ice_models.R,
#           08_results_gams/results_tables.R and 08_results_gams/response_curves.R
#           after config.R. Do not run on its own.
#
#           Environment-variable overrides (testing only; production runs use
#           the defaults from config.R):
#             ODONTOCETE_OUTPUT_ROOT  - root folder for all outputs
#             ODONTOCETE_ENV_DAILY    - path to env_daily_wide.csv
#             ODONTOCETE_MODEL_INPUT  - path to model_input_daily.csv
#             ODONTOCETE_ALLOW_SHORT_CLIMATOLOGY - "TRUE" lets a climatology
#               shorter than CLIM_MIN_YEARS run (with a loud warning); this
#               exists ONLY to test the code before the 10-yr downloads exist.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Input / output locations
# -----------------------------------------------------------------------------
# Defaults are the data-contract files in config.R. The overrides let the
# workflow be tested on a temporary copy of the data without touching
# final_scripts_odontocetes/outputs.
OUTPUT_ROOT       <- Sys.getenv("ODONTOCETE_OUTPUT_ROOT", unset = OUT_DIR)
GAM_ENV_DAILY     <- Sys.getenv("ODONTOCETE_ENV_DAILY",   unset = F_ENV_DAILY)
GAM_MODEL_INPUT   <- Sys.getenv("ODONTOCETE_MODEL_INPUT", unset = F_MODEL_INPUT)
ALLOW_SHORT_CLIMATOLOGY <- toupper(Sys.getenv("ODONTOCETE_ALLOW_SHORT_CLIMATOLOGY",
                                              unset = "FALSE")) == "TRUE"

STATS_DIR   <- file.path(OUTPUT_ROOT, "06_statistical_analyses")
RESULTS_DIR <- file.path(OUTPUT_ROOT, "08_results_gams")

# Folder for one stage's outputs ("stage1", "stage2", "stage3", "sea_ice"),
# created on use.
stage_dir <- function(stage, ...) {
  d <- file.path(STATS_DIR, stage, ...)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

# -----------------------------------------------------------------------------
# 2. Model family and smoother settings (kept as in the original pipeline)
# -----------------------------------------------------------------------------
# Tweedie with log link; the power parameter p is estimated within [1.1, 1.9]
# (between Poisson-like and gamma-like) because the response (mean daily
# presence minutes per bin) is non-negative, continuous and zero-inflated.
gam_family <- function() mgcv::tw(link = "log", a = 1.1, b = 1.9)
GAM_METHOD <- "REML"
GAM_K      <- 4      # basis dimension of every s(x): limits curves to simple shapes
P_THRESH   <- 0.05   # backward elimination: drop terms until all p <= 0.05
VIF_THRESH <- 5      # sequential VIF screening threshold (GVIF^(1/(2*Df)))

# Cyclic predictors: s(x, bs = "cc", k = GAM_K) with the ends joined at the
# given knots. FSLE orientation is an axial angle (eigenvector orientation,
# theta == theta + 180 deg; Haller 2015, Annu. Rev. Fluid Mech., Eq. 7) in
# (-90, 90] deg anticlockwise from East, so -90 and 90 deg are the same
# direction. Used in every GAM fit (theme screening, final models, Stage 3)
# and in the response curves (evaluated over the full -90..90 circle).
CYCLIC_VARS <- list(fsle_orient = c(-90, 90))

# -----------------------------------------------------------------------------
# 3. Temporal binning (residual autocorrelation)
# -----------------------------------------------------------------------------
# Block model fitted to DAILY presence minutes per species x site; the bin
# length is the first lag (>= 1 day) at which the residual ACF falls inside the
# 95% white-noise bounds. `day_index` = days since the first recording day, so
# the smooth absorbs the within-deployment seasonal trend (julian_day is not
# used because it wraps at KGI, which spans Feb 2015 - Jan 2016).
# Confirmed by Natalie (2026-11-17); the original baleen pipeline used
# site-specific environmental terms with a fixed sp = 0.1.
ACF_BLOCK_RHS   <- "s(day_index, k = 4)"
ACF_LAG_MAX     <- 30     # days examined in the block-model residual ACF
ACF_CI_LEVEL    <- 0.95
BIN_CAP_DAYS    <- 4      # cap: >= 30 bins at the shortest deployment (EI, 135 d)
BIN_MIN_DAYS    <- 1      # bins with fewer recording days are dropped [CONFIRM]

# -----------------------------------------------------------------------------
# 4. Inclusion rule and term caps (bins with detections; counted after the ice
#    filter for habitat models, over all bins for the sea-ice models)
# -----------------------------------------------------------------------------
# Testing overrides (sensitivity runs only; production uses the defaults):
#   ODONTOCETE_MIN_DETECTION_BINS, ODONTOCETE_LOW_SAMPLE_MAX_TERMS ("Inf" =
#   no low-sample cap), ODONTOCETE_MAX_TERMS_PER_BINS ("1" = effectively no
#   term cap).
env_num <- function(name, default) as.numeric(Sys.getenv(name, unset = default))
MIN_DETECTION_BINS  <- env_num("ODONTOCETE_MIN_DETECTION_BINS", 10)  # < 10 -> not modeled
LOW_SAMPLE_MAX_BINS <- 20  # 10-20 -> modeled, flagged low_sample
LOW_SAMPLE_MAX_TERMS <- env_num("ODONTOCETE_LOW_SAMPLE_MAX_TERMS", 2)  # max terms, low sample
# Habitat-model term cap: at most floor(n_bins / MAX_TERMS_PER_BINS) smooth
# terms per species x site (n_bins = all bins in the model); low-sample models
# use the smaller of this cap and LOW_SAMPLE_MAX_TERMS.
MAX_TERMS_PER_BINS <- env_num("ODONTOCETE_MAX_TERMS_PER_BINS", 10)

# -----------------------------------------------------------------------------
# 5. Deseasoning (10-year monthly climatology, periodic cubic spline)
# -----------------------------------------------------------------------------
# Climatology window per site: deployment start year - CLIM_HALF_WINDOW_YEARS
# to deployment end year + CLIM_HALF_WINDOW_YEARS (config.R). A variable must
# have data in >= CLIM_MIN_YEARS distinct years and in all 12 calendar months.
CLIM_MIN_YEARS <- 10
# Day of year of the 15th of each month (non-leap reference year): anchors of
# the 12 monthly means for the periodic spline.
MID_MONTH_YDAY <- as.integer(format(as.Date(sprintf("2001-%02d-15", 1:12)), "%j"))

# Raw OR anomaly (never both) - decision 2.6, revised 2026-11-19:
#   The sea-ice models describe the seasonal pattern of presence; the habitat
#   models test responses to conditions that are unusual for the time of
#   year. Every habitat candidate is therefore the ANOMALY (<var>_anom) from
#   the 10-yr climatology, except the mesoscale and front variables, which
#   have almost no seasonal cycle (seasonal share of daily variance <= 0.08
#   at all sites) and stay RAW: EKE_<d>, FSLE, fsle_orient.
#   (Before 2026-11-19 the anomaly was used only when the seasonal cycle
#   explained > 50% of daily variance, and chl-a/productivity were always
#   raw; raw chl-a then acted as a proxy for season.)
#   The seasonal share (1 - var(daily - spline climatology) / var(daily)) is
#   still computed and reported in seasonal_variance.csv.
ANOM_SURFACE_VARS <- c("mixed_layer")
# Layer-averaged water-column variables (decision 2026-11-19; layers in
# config.R LAYERS: surface 0-100 m, mid-water 100-500 m, deep 500-1000 m).
# Each species uses the surface layer plus the layers of its foraging depths
# (config.R species_layers()).
#   * Water mass of each layer: the temperature anomaly of the layer, one
#     theme per layer (decision 2026-11-19). Temperature is the only
#     water-mass variable: against World Ocean Database profiles the
#     reanalysis reproduced layer-mean temperature well (r = 0.75-0.90 in the
#     three layers) but not mid-water salinity (r = 0.19) or oxygen (r <= 0.15;
#     05b_hycom_comparison/06_compare_models_insitu.R). Salinity and oxygen
#     are therefore not candidates; ODONTOCETE_WATER_MASS_VARS can restore
#     them for sensitivity runs.
#   * Eddy activity: EKE of the mid and deep layers, one theme per layer.
#     Surface mesoscale structure is represented by FSLE and FSLV orientation
#     (satellite altimetry), as in the companion baleen whale analysis.
#   Temperature enters as anomalies; EKE stays raw.
# Testing override: ODONTOCETE_EXCLUDE_VARS (comma-separated candidate
# source columns, e.g. "chla_lag3") removes candidates (sensitivity runs only).
EXCLUDE_VARS <- strsplit(Sys.getenv("ODONTOCETE_EXCLUDE_VARS", unset = ""), ",")[[1]]
# Testing override: ODONTOCETE_WATER_MASS_VARS (comma-separated, e.g.
# "temperature") restricts the water-mass candidates (sensitivity runs only).
WATER_MASS_BASES  <- strsplit(Sys.getenv("ODONTOCETE_WATER_MASS_VARS",
                                         unset = "temperature"), ",")[[1]]
SUBSURFACE_ONLY   <- c("o2", "EKE")      # not used in the surface layer
LAYER_BASES       <- c(WATER_MASS_BASES, "EKE")
ANOM_LAYER_BASES  <- WATER_MASS_BASES
# Lagged chl-a is log10-transformed before deseasoning (chl-a is strongly
# right-skewed); the anomaly is then a departure in log10 chl-a.
CHL_LOG10         <- TRUE

# -----------------------------------------------------------------------------
# 6. Candidate predictors and themes
# -----------------------------------------------------------------------------
# Surface predictors used by every species. SSH anomaly and within-box flow
# variability (EKE_mad) were dropped as candidates (decision 2026-11-19).
SURFACE_VARS      <- c("mixed_layer", "FSLE", "fsle_orient")
# Sea ice is NOT a habitat-model candidate (decision 2.6, 2026-11-17): it
# defines the ice filter and is modeled on its own in 04_sea_ice_models.R
# (presence ~ s(ice_conc), all bins). ice_thickness is not used in any model.
ICE_MODEL_VAR     <- "ice_conc"
# Signed distance from the recorder to the ice edge (15%; km; < 0 = inside the
# ice; AMSR2, 03_sea_ice/amsr2_site_series.R). Second candidate in the
# sea-ice models at the sites with an ice season during recording (decision
# N. Posdaljian 2026-11-20; EI recorded only before ice advance).
ICE_EDGE_VAR      <- "dist_ice_edge_km"
ICE_EDGE_SITES    <- c("KGI", "CI")
# Productivity theme = lagged chl-a only (decisions 2026-11-19): the 30-day
# window mean of surface chl-a ending 3 months before each day, one common
# lag for all species and sites (chosen over 5 months: significant alone in
# 9 vs 6 of 11 Stage 1 models, mean deviance explained 37% vs 24%).
# Same-day chl-a is not used (odontocete prey respond to blooms with a delay)
# and net primary productivity is dropped (largely redundant with chl-a).
CHL_LAG_MONTHS    <- 3
LAG_COLUMNS       <- paste0("chla_lag", CHL_LAG_MONTHS)
# Standardized chl-a anomaly: divide the
# anomaly by its seasonal SD (smooth monthly SD of the anomaly over the
# 10-year record), so anomalies are comparable across months (Wilks 2011).
# Adopted 2026-11-20 (decision N. Posdaljian); set the env variable to FALSE
# to reproduce the unstandardized run.
CHL_STANDARDIZE   <- toupper(Sys.getenv("ODONTOCETE_CHL_STANDARDIZE", unset = "TRUE")) == "TRUE"

# One winner per theme enters VIF screening and backward elimination. Themes
# come from the base variable name and, for layer variables, the layer
# (parse_variable()): water_mass_<layer> (temperature / salinity / oxygen),
# eddy_<layer> (EKE), surface_fronts (FSLE, satellite altimetry),
# front_orientation (FSLV), stratification (mixed layer depth) and
# productivity (lagged chl-a).
THEME_OF_BASE <- c(
  FSLE         = "surface_fronts",
  fsle_orient  = "front_orientation",
  chla         = "productivity",
  mixed_layer  = "stratification"
)
THEME_ORDER <- c("water_mass_surface", "water_mass_mid", "water_mass_deep",
                 "surface_fronts", "front_orientation", "eddy_mid", "eddy_deep",
                 "stratification", "productivity")

# -----------------------------------------------------------------------------
# 7. Diagnostics and shape classification
# -----------------------------------------------------------------------------
RANDOM_SEED <- 20261116          # k-index p-values are permutation based
# Residual ACF of the final model (bin scale, ice-removed bins kept as gaps),
# plotted to RESID_ACF_MAX_LAG. Residual autocorrelation is tested jointly
# over lags 1..LJUNG_BOX_LAGS with the Ljung-Box test (decision 2026-11-19):
# "OK" when p > ACF_TEST_ALPHA. (The earlier rule, "no single lag in 2..10
# outside the 95% bounds", fails ~37% of independent series by chance.)
RESID_ACF_MAX_LAG <- 10
LJUNG_BOX_LAGS    <- 5
ACF_TEST_ALPHA    <- 0.05
# Residual ACF results are reported as a limitation (bins stay capped at 4 d
# and the selection threshold stays p <= 0.05; decision 2.6).
# Smooth shapes: terms with p > P_THRESH are "not significant"; otherwise the
# fitted centered curve (no CI band) is classified on a grid of SHAPE_GRID_N
# points spanning the observed range (classify_smooth_shape()).
SHAPE_GRID_N        <- 100
SHAPE_MONO_TOL      <- 0.10  # share of total variation allowed against the trend
SHAPE_SLOPE_NEGLIG  <- 0.05  # |slope| below this share of max |slope| is ignored

# -----------------------------------------------------------------------------
# A priori candidate sets for the habitat models (decision N. Posdaljian
# 2026-11-20). Replaces the Stage 1 union sets: each species' candidates are
# defined from the literature and its foraging layers, not from the data.
#   All species : surface fronts (FSLE, FSLV orientation), stratification
#                 (mixed layer depth anomaly)
#   + water mass (temperature anomaly) and eddy activity (EKE) in the layers
#     each species forages in; + surface temperature for the shallower
#     foragers (pilot and killer whales)
# Chl-a removed from the main sets (decision N. Posdaljian 2026-11-25, PI
# feedback): no mechanistic basis for a specific lag at 3-4 trophic levels
# above phytoplankton; the 3-month lag was chosen from the same data. Not
# reported anywhere in the paper (decision 2026-11-25); the unreported check
# is 09_chla_sensitivity.R (none / lag 0 / lag 3; no consistent signal).
# Pilot whales use the mid-water layer (was deep; decision 2026-11-25).
# Used by 03_stage3_harmonized_models.R and 07_permutation_hypotheses.R (C).
# -----------------------------------------------------------------------------
COMMON_SET <- c("FSLE", "fsle_orient", "mixed_layer_anom")
HYPOTHESIS_SETS <- list(
  Pm   = c(COMMON_SET, "temperature_mid_anom", "temperature_deep_anom", "EKE_mid", "EKE_deep"),
  Gm   = c(COMMON_SET, "temperature_surface_anom", "temperature_mid_anom", "EKE_mid"),
  Oo   = c(COMMON_SET, "temperature_surface_anom", "temperature_mid_anom", "EKE_mid"),
  BW29 = c(COMMON_SET, "temperature_deep_anom", "EKE_deep"),
  BW37 = c(COMMON_SET, "temperature_deep_anom", "EKE_deep")
)
