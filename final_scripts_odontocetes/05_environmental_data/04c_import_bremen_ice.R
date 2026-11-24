# =============================================================================
# Script  : 04c_import_bremen_ice.R
# Section : 05_environmental_data (sea ice, Section 2.3)
# Purpose : Import daily AMSR2 ASI sea-ice concentration (University of Bremen,
#           3.125 km grid) box means for each site, check them against the
#           recording periods, and write the daily sea-ice table used by the
#           models when SEA_ICE_SOURCE == "bremen" (config.R).
# Inputs  : final_scripts_odontocetes/config.R (SITES, BREMEN_ICE_DIR,
#           BREMEN_MAX_GAP_DAYS)
#           <BREMEN_ICE_DIR>/<site>/*.csv  extracted box means, one or more
#             files per site (e.g., one per year from EvaSeaIce.py). Columns:
#             a date column (Date/date: yyyy-mm-dd, yyyymmdd, or MATLAB datenum)
#             and one concentration column (bs_mean, ice_conc, sic,
#             concentration; percent 0-100 or fraction 0-1).
#           F_COP_SURFACE (optional; for the Bremen vs Copernicus check)
# Outputs : F_BREMEN_ICE  bremen_ice_daily.csv
#             date, site, ice_conc (fraction 0-1), filled, source_file
#           outputs/05_environmental_data/bremen_ice_gaps.csv
#             site, gap_start, gap_end, n_days, filled
#           outputs/05_environmental_data/bremen_vs_copernicus_ice.csv
#             site, n_days, bias (Bremen - Copernicus), rmse, r,
#             mean_bremen, mean_copernicus
#           outputs/05_environmental_data/bremen_vs_copernicus_ice.png
# Author  : Natalie Posdaljian
# Date    : 2026-11-23
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - The Bremen files are NOT downloaded or extracted here; they were
#             extracted on another computer (EvaSeaIce.py / SeaIceExtraction.m
#             in Environmental Data/Sea Ice Bremen/). [CONFIRM extraction
#             area: EvaSeaIce.py averaged over a 30 km circle (QGIS *30.shp),
#             excluding land (255) and the 120 flag; SeaIceExtraction.m used
#             the single nearest pixel. All other predictors use the +/- 20 km
#             box.]
#           - Units: values > 1.5 anywhere in a site's record are taken as
#             percent and divided by 100, so ice_conc is a fraction like the
#             Copernicus siconc it replaces.
#           - Only days inside each site's recording period (calendar days of
#             rec_start to rec_end) are kept. Missing days (e.g., no AMSR2 file
#             on the server) are linearly interpolated when the gap is
#             <= BREMEN_MAX_GAP_DAYS; longer gaps stop the script.
#           - Duplicate dates across files must agree (|difference| < 1e-6);
#             otherwise the script stops.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))

library(tidyverse)

F_GAPS    <- out_path("05_environmental_data", "bremen_ice_gaps.csv")
F_CMP_CSV <- out_path("05_environmental_data", "bremen_vs_copernicus_ice.csv")
F_CMP_PNG <- out_path("05_environmental_data", "bremen_vs_copernicus_ice.png")

VALUE_COLS <- c("bs_mean", "ice_conc", "sic", "concentration", "ice", "mean")

# -----------------------------------------------------------------------------
# 1. Read the extracted files
# -----------------------------------------------------------------------------
if (!dir.exists(BREMEN_ICE_DIR)) {
  stop("Bremen sea-ice folder not found: ", BREMEN_ICE_DIR, "\n",
       "Copy the extracted CSVs into ", BREMEN_ICE_DIR, "/<site>/ (EI, KGI, CI), ",
       "or set SEA_ICE_SOURCE <- \"copernicus\" in config.R.")
}

parse_dates <- function(x) {
  if (is.numeric(x)) {
    if (all(x > 7e5, na.rm = TRUE)) return(as.Date(x - 719529, origin = "1970-01-01"))  # MATLAB datenum
    x <- as.character(x)                                                               # yyyymmdd
  }
  x <- as.character(x)
  d <- suppressWarnings(as.Date(lubridate::parse_date_time(x, c("Ymd", "Ymd HMS", "mdY", "dmY"))))
  if (any(is.na(d) & !is.na(x))) stop("Unparseable dates, e.g. '", x[is.na(d)][1], "'")
  d
}

read_one <- function(path, site) {
  df <- read_csv(path, show_col_types = FALSE, guess_max = 10000)
  names(df) <- str_trim(names(df))
  date_col <- names(df)[tolower(names(df)) %in% c("date", "time", "day", "datenum")][1]
  if (is.na(date_col)) date_col <- names(df)[1]
  val_col <- names(df)[tolower(names(df)) %in% VALUE_COLS][1]
  if (is.na(val_col)) {
    num <- setdiff(names(df)[map_lgl(df, is.numeric)], date_col)
    if (length(num) != 1) {
      stop(basename(path), ": cannot tell which column is sea-ice concentration (columns: ",
           paste(names(df), collapse = ", "), ")")
    }
    val_col <- num
  }
  tibble(site = site, date = parse_dates(df[[date_col]]),
         value = as.numeric(df[[val_col]]),
         source_file = file.path(site, basename(path)))
}

raw <- map_dfr(SITE_ORDER, function(s) {
  files <- list.files(file.path(BREMEN_ICE_DIR, s), pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No CSV files for site ", s, " in ", file.path(BREMEN_ICE_DIR, s))
  message(sprintf("%s: %d file(s): %s", s, length(files), paste(basename(files), collapse = ", ")))
  map_dfr(files, read_one, site = s)
})

# -----------------------------------------------------------------------------
# 2. Units, duplicates, recording period
# -----------------------------------------------------------------------------
raw <- raw %>%
  group_by(site) %>%
  mutate(percent = max(value, na.rm = TRUE) > 1.5,
         ice_conc = if_else(percent, value / 100, value)) %>%
  ungroup()
units <- distinct(raw, site, percent)
message("Units detected (TRUE = percent, divided by 100):")
print(units)
bad <- filter(raw, !is.na(ice_conc), ice_conc < 0 | ice_conc > 1)
if (nrow(bad) > 0) stop(nrow(bad), " values outside 0-1 after unit conversion, e.g. ",
                        bad$site[1], " ", bad$date[1], " = ", bad$value[1])

dups <- raw %>%
  filter(!is.na(ice_conc)) %>%
  group_by(site, date) %>%
  filter(n() > 1, diff(range(ice_conc)) > 1e-6) %>%
  ungroup()
if (nrow(dups) > 0) {
  print(dups)
  stop("Conflicting values for the same site and date in different files (see above)")
}
ice <- raw %>%
  filter(!is.na(ice_conc)) %>%
  group_by(site, date) %>%
  summarise(ice_conc = first(ice_conc), source_file = first(source_file), .groups = "drop")

periods <- SITES %>%
  transmute(site, start = as.Date(rec_start), end = as.Date(rec_end))

daily <- periods %>%
  rowwise() %>%
  reframe(site, date = seq(start, end, by = "day")) %>%
  left_join(ice, by = c("site", "date")) %>%
  arrange(factor(site, levels = SITE_ORDER), date)

# -----------------------------------------------------------------------------
# 3. Gaps inside the recording periods
# -----------------------------------------------------------------------------
gaps <- daily %>%
  group_by(site) %>%
  mutate(run = cumsum(!is.na(ice_conc) != lag(!is.na(ice_conc), default = TRUE))) %>%
  filter(is.na(ice_conc)) %>%
  group_by(site, run) %>%
  summarise(gap_start = min(date), gap_end = max(date), n_days = n(), .groups = "drop") %>%
  select(-run) %>%
  left_join(periods, by = "site") %>%
  # A gap at the very start or end of a period cannot be interpolated.
  mutate(at_edge = gap_start == start | gap_end == end,
         filled = n_days <= BREMEN_MAX_GAP_DAYS & !at_edge) %>%
  select(site, gap_start, gap_end, n_days, filled)
write_csv(gaps, F_GAPS)

if (nrow(gaps) > 0) {
  message("Days without Bremen data inside the recording periods:")
  print(gaps, n = Inf)
}
if (any(!gaps$filled)) {
  stop(sum(!gaps$filled), " gap(s) longer than BREMEN_MAX_GAP_DAYS = ", BREMEN_MAX_GAP_DAYS,
       " days (or at the start/end of a recording period); see ", F_GAPS,
       ". Add the missing days to the Bremen files or raise BREMEN_MAX_GAP_DAYS.")
}

daily <- daily %>%
  group_by(site) %>%
  mutate(filled = is.na(ice_conc),
         ice_conc = approx(as.numeric(date[!filled]), ice_conc[!filled],
                           xout = as.numeric(date))$y,
         source_file = if_else(filled, "interpolated", source_file)) %>%
  ungroup() %>%
  select(date, site, ice_conc, filled, source_file)

write_csv(daily, F_BREMEN_ICE)
message(sprintf("Wrote %s: %d site-days (%d interpolated)", F_BREMEN_ICE,
                nrow(daily), sum(daily$filled)))
print(daily %>% group_by(site) %>%
        summarise(days = n(), interpolated = sum(filled),
                  mean_ice = mean(ice_conc), max_ice = max(ice_conc),
                  days_ge_15pct = sum(ice_conc >= 0.15)))

# -----------------------------------------------------------------------------
# 4. Check against Copernicus siconc (same days)
# -----------------------------------------------------------------------------
if (file.exists(F_COP_SURFACE)) {
  cop <- read_csv(F_COP_SURFACE, show_col_types = FALSE) %>%
    transmute(date = as.Date(date), site, copernicus = ice_conc)
  both <- daily %>%
    transmute(date, site, bremen = ice_conc) %>%
    inner_join(cop, by = c("date", "site"))
  cmp <- both %>%
    group_by(site) %>%
    summarise(n_days = n(),
              bias = mean(bremen - copernicus),
              rmse = sqrt(mean((bremen - copernicus)^2)),
              r = if (sd(bremen) > 0 && sd(copernicus) > 0) cor(bremen, copernicus) else NA_real_,
              mean_bremen = mean(bremen), mean_copernicus = mean(copernicus)) %>%
    arrange(factor(site, levels = SITE_ORDER))
  write_csv(cmp, F_CMP_CSV)
  message("Bremen vs Copernicus sea-ice concentration (fraction):")
  print(cmp)

  p <- both %>%
    pivot_longer(c(bremen, copernicus), names_to = "product", values_to = "ice_conc") %>%
    mutate(site = factor(site, levels = SITE_ORDER),
           product = recode(product, bremen = "AMSR2 ASI (Bremen)",
                            copernicus = "GLORYS12 (Copernicus)")) %>%
    ggplot(aes(date, 100 * ice_conc, colour = product)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ site, ncol = 1, scales = "free_x") +
    scale_colour_manual(values = c("#1E88E5", "#D81B60")) +
    labs(x = NULL, y = "Sea-ice concentration (%)", colour = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top")
  ggsave(F_CMP_PNG, p, width = 7, height = 6, dpi = 300)
  message("Wrote ", F_CMP_PNG)
} else {
  message("Skipping the Copernicus comparison: ", F_COP_SURFACE, " not found")
}
