# =============================================================================
# Script  : _values.R
# Purpose : Read the analysis outputs and define every value, phrase and table
#           that manuscript.qmd prints, so the manuscript updates itself when
#           the analysis scripts are re-run.
# Inputs  : CSVs in final_scripts_odontocetes/outputs/ (paths relative to
#           manuscript/, where Quarto runs this file)
# Outputs : objects used inline in manuscript.qmd (`r ...`) and table functions
# Author  : Natalie Posdaljian
# Date    : 2026-11-19
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(knitr); library(kableExtra)
})

OUT <- "../outputs"
rd <- function(...) read_csv(file.path(OUT, ...), show_col_types = FALSE)

SPECIES <- c("Pm", "Gm", "Oo", "BW29", "BW37")
SITES <- c("EI", "KGI", "CI")
SP_TEXT <- c(Pm = "sperm whales", Gm = "long-finned pilot whales", Oo = "killer whales",
             BW29 = "southern bottlenose whales", BW37 = "BW37")
SP_TABLE <- c(Pm = "Sperm whale", Gm = "Long-finned pilot whale", Oo = "Killer whale",
              BW29 = "Southern bottlenose whale", BW37 = "BW37 beaked whale")
TERM_LABEL <- c(FSLE = "FSLE", fsle_orient = "FSLV orientation",
                EKE_mid = "EKE, mid-water layer", EKE_deep = "EKE, deep layer",
                mixed_layer_anom = "Mixed layer depth anomaly",
                temperature_surface_anom = "Temperature anomaly, surface layer",
                temperature_mid_anom = "Temperature anomaly, mid-water layer",
                temperature_deep_anom = "Temperature anomaly, deep layer",
                chla_lag3_anom = "Chl-a standardized anomaly, 3-month lag")

# ---- Formatting helpers ------------------------------------------------------
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
fp <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
and_list <- function(x) {
  if (length(x) <= 1) return(paste(x))
  if (length(x) == 2) return(paste(x, collapse = " and "))
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}
order_rows <- function(d) {
  if (!"site" %in% names(d)) return(d[order(match(d$species, SPECIES)), ])
  d[order(match(d$species, SPECIES), match(d$site, SITES)), ]
}
# "sperm whales at KGI, killer whales at EI and KGI, and BW37 at CI"
combo_phrase <- function(d) {
  d <- order_rows(d)
  parts <- vapply(unique(d$species), function(s)
    paste(SP_TEXT[[s]], "at", and_list(d$site[d$species == s])), "")
  and_list(parts)
}

# ---- Sea-ice cutoffs and bins removed -----------------------------------------
ice <- rd("06_statistical_analyses", "stage1", "tables", "ice_cutoffs.csv")
cutoff <- with(distinct(ice, species, ice_cutoff), setNames(100 * ice_cutoff, species))
removed <- function(sp, st) with(filter(ice, species == sp, site == st),
                                 paste(n_bins_removed, "of", n_bins_before))

# ---- Bin lengths -------------------------------------------------------------
bins <- rd("06_statistical_analyses", "stage1", "tables", "bin_lengths.csv") |> order_rows()
capped <- filter(bins, capped)
capped_phrase <- and_list(SP_TEXT[capped$species])
capped_lags <- and_list(capped$max_site_lag)

# ---- Model inclusion -----------------------------------------------------------
hab_models <- rd("08_results_gams", "main_stage3_habitat_gam_table.csv")
not_modeled <- filter(hab_models, status == "descriptive only")
low_sample <- filter(hab_models, low_sample)
n_models <- sum(startsWith(hab_models$status, "modeled"))

# Seasonal share of EKE, FSLE and FSLV orientation (<= 8%) is not in the
# outputs (exploratory check, 2026-11-19); it is typed in manuscript.qmd.

# ---- Habitat model terms and diagnostics -------------------------------------
hab_terms <- rd("08_results_gams", "main_stage3_habitat_gam_table_terms.csv")
edf_median <- f1(median(hab_terms$edf)); edf_min <- f1(min(hab_terms$edf)); edf_max <- f1(max(hab_terms$edf))
hab_diag <- rd("08_results_gams", "supplementary_stage3_habitat_gam_diagnostics.csv") |> order_rows()
n_lb_fail <- sum(hab_diag$ljung_box_p <= 0.05)
lb_fail_phrase <- combo_phrase(filter(hab_diag, ljung_box_p <= 0.05))

ar1_models <- rd("06_statistical_analyses", "stage3", "tables", "ar1_sensitivity_models.csv")
ar1_terms <- rd("06_statistical_analyses", "stage3", "tables", "ar1_sensitivity_terms.csv")
n_ar1_supported <- sum(ar1_models$ar1_supported, na.rm = TRUE)
ar1_supported_phrase <- combo_phrase(filter(ar1_models, ar1_supported))
n_terms_testable <- sum(!is.na(ar1_terms$holds))
n_terms_hold <- sum(ar1_terms$holds, na.rm = TRUE)
n_terms_untestable <- sum(is.na(ar1_terms$holds))

# ---- Model validation ----------------------------------------------------------
lay <- rd("05b_hycom_comparison", "insitu_metrics_by_layer.csv")
lv <- function(m, v, l, col) lay[[col]][lay$model == m & lay$variable == v & lay$layer == l]
r_T <- sapply(c("surface", "mid", "deep"), function(l) lv("Copernicus", "temperature", l, "r"))
r_T_range <- paste0(f2(min(r_T)), "--", f2(max(r_T)))
r_S_mid_cop <- f2(lv("Copernicus", "salinity", "mid", "r"))
r_S_mid_hyc <- f2(lv("HYCOM", "salinity", "mid", "r"))
r_O2_max <- f2(max(lay$r[lay$variable == "oxygen"]))
bias_T <- function(m, l) f2(lv(m, "temperature", l, "bias"))

cmp <- rd("05b_hycom_comparison", "comparison_metrics.csv") |> filter(values == "raw")
modis_bias <- function(prefix) {
  b <- cmp$bias[startsWith(cmp$comparison, prefix) & grepl("MODIS", cmp$comparison)]
  paste0(f2(min(b)), " to ", f2(max(b)))
}

# ---- Occurrence across sea-ice classes -----------------------------------------
ice_occ <- rd("08_results_gams", "main_sea_ice_occurrence_table.csv")
heavy_ice_species <- ice_occ |> filter(ice_class == ">50", n_detection_bins > 0) |>
  distinct(species) |> pull(species)

# ---- Tables ------------------------------------------------------------------
# `hlines`: rows after which a horizontal rule is drawn (e.g. between species).
# "%" in captions is escaped (in LaTeX it starts a comment and cuts the caption).
latex_table <- function(d, caption, widths = NULL, size = 8, long = FALSE, align = NULL,
                        hlines = NULL) {
  k <- kbl(d, format = "latex", booktabs = TRUE, longtable = long,
           caption = gsub("%", "\\%", caption, fixed = TRUE),
           align = align, linesep = "")
  k <- kable_styling(k, font_size = size,
                     latex_options = c(if (long) "repeat_header" else "hold_position"))
  if (!is.null(widths)) for (i in seq_along(widths)) if (!is.na(widths[i]))
    k <- column_spec(k, i, width = widths[i])
  if (length(hlines) > 0) k <- row_spec(k, hlines, hline_after = TRUE)
  k
}
# Species shown once per species; site and model columns once per model.
# The rows ending each species are kept as attribute "species_end".
blank_repeats <- function(d) {
  species_end <- which(d$Species != dplyr::lead(d$Species))
  first_model <- !duplicated(paste(d$Species, d$Site))
  first_species <- !duplicated(d$Species)
  d$Species[!first_species] <- ""
  for (col in intersect(c("Site", "Bins (det.)", "Dev. (%)"), names(d))) d[[col]][!first_model] <- ""
  attr(d, "species_end") <- species_end
  d
}

habitat_table <- function() {
  modeled <- hab_terms |>
    mutate(Term = paste0(TERM_LABEL[variable], ifelse(robust_ar1 %in% FALSE, "†", "")))
  rows <- bind_rows(
    modeled |> transmute(species, site, n_bins, n_detection_bins, low_sample,
                         dev = f1(dev_expl_model_pct), Term, edf = f2(edf),
                         p = p_value, tdev = f1(dev_expl_term_pct), Shape = shape),
    not_modeled |> transmute(species, site, n_bins, n_detection_bins, low_sample = FALSE,
                             dev = "--", Term = "not modeled (<10 bins with detections)",
                             edf = "", p = "", tdev = "", Shape = "")
  )
  rows |>
    mutate(ord = row_number()) |>
    order_rows() |>
    transmute(Species = SP_TABLE[species], Site = site,
              `Bins (det.)` = paste0(n_bins, " (", n_detection_bins, ")", ifelse(low_sample, "‡", "")),
              `Dev. (%)` = dev, Term, edf, p, `Term dev. (%)` = tdev, Shape) |>
    blank_repeats()
}

sea_ice_table <- function() {
  rd("08_results_gams", "main_sea_ice_gam_table.csv") |> order_rows() |>
    transmute(Species = SP_TABLE[species], Site = site,
              `Bins (det.)` = paste0(n_bins, " (", n_detection_bins, ")", ifelse(low_sample, "‡", "")),
              edf = ifelse(is.na(edf), "--", f2(edf)), p = fp(p_value),
              `Dev. (%)` = ifelse(is.na(dev_expl_pct), "--", f1(dev_expl_pct)),
              Shape = case_when(status == "descriptive only" ~ "not modeled",
                                shape == "not significant" ~ paste(fitted_shape, "(n.s.)"),
                                TRUE ~ shape),
              `Robust to AR(1)` = ifelse(is.na(robust_ar1), "", ifelse(robust_ar1, "yes", "no")))
}

diagnostics_table <- function() {
  hab_diag |>
    transmute(Species = SP_TABLE[species], Site = site, `Bin (d)` = bin_days,
              `Tweedie p` = f2(tweedie_p), `Min k-index p` = formatC(k_pvalue_min, format = "f", digits = 3),
              `Max concurvity` = f2(concurvity_max), `Lag-1 ACF` = f2(resid_acf_lag1),
              `Ljung–Box p` = fp(ljung_box_p),
              `ΔAIC with AR(1)` = ifelse(is.na(delta_aic_ar1), "--", f1(delta_aic_ar1)))
}

ice_cutoff_table <- function() {
  ice |>
    mutate(cell = ifelse(species == "Pm" & site == "KGI", "-- (not detected)",
                         paste0(n_bins_removed, "/", n_bins_before))) |>
    select(species, site, cell) |>
    tidyr::pivot_wider(names_from = site, values_from = cell) |>
    order_rows() |>
    transmute(Species = SP_TABLE[species], `Cutoff (% ice)` = f1(cutoff[species]),
              `EI removed/total` = EI, `KGI removed/total` = KGI, `CI removed/total` = CI)
}

duty_table <- function() {
  simplify <- function(x) sub("^100\\.0 \\(100\\.0–100\\.0\\)$", "100.0", x)
  read_csv(file.path(OUT, "04_acoustic_detection", "duty_cycle_simulation_table.csv"),
           show_col_types = FALSE) |>
    transmute(Species, `Duty cycle`,
              Events = simplify(`% events retained`),
              `5-min bins` = simplify(`% 5-min bins retained`),
              `Presence days` = simplify(`% presence days retained`))
}

cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))

# ---- AR(1) check: terms that did not hold ------------------------------------
strip_s <- function(x) gsub("^s[(]|[)]$", "", x)   # "s(FSLE)" -> "FSLE"
term_text <- function(v) {   # lower-case the first word unless it is an acronym
  x <- sub("Chl-a", "chlorophyll-a", TERM_LABEL[v])
  unname(sub("^([A-Z])([a-z])", "\\L\\1\\2", x, perl = TRUE))
}

# ---- Site/year comparison (05_environmental_data/07_site_year_comparison.R) ----
site_year <- rd("05_environmental_data", "site_year_variance.csv")
site_year_pct <- function(v, col) f1(site_year[[col]][site_year$variable == v])
site_year_table <- function() {
  labels <- c(ice_conc = "Sea-ice concentration", temperature_surface = "Temperature, surface layer",
              temperature_mid = "Temperature, mid-water layer", temperature_deep = "Temperature, deep layer",
              mixed_layer = "Mixed layer depth", FSLE = "FSLE", EKE_mid = "EKE, mid-water layer",
              EKE_deep = "EKE, deep layer", chla_lag3 = "log10 chlorophyll-a, 3-month lag",
              fsle_orient_cos = "FSLV orientation (cos 2θ)", fsle_orient_sin = "FSLV orientation (sin 2θ)")
  site_year |>
    transmute(Variable = labels[variable], `Site (%)` = f1(pct_site), `Year (%)` = f1(pct_year),
              `Site × year (%)` = f1(pct_site_x_year), `Residual (%)` = f1(pct_residual),
              `Seasonal share of total (%)` = f1(pct_total_seasonal))
}
not_holding <- ar1_terms |>
  filter(holds %in% FALSE) |>
  mutate(variable = strip_s(term),
         reason = ifelse(!is.na(p_tmb) & p_tmb > 0.05, "ml", "ar1")) |>
  order_rows()
holding_phrase <- function(why) {
  d <- filter(not_holding, reason == why)
  if (nrow(d) == 0) return("none")
  and_list(paste0(term_text(d$variable), " (", SP_TEXT[d$species], " at ", d$site, ")"))
}

# ---- Pairwise concurvity of the final habitat models (worst case) -------------
CONCURVITY_FLAG <- 0.5
concurvity_pairs <- local({
  files <- list.files(file.path(OUT, "06_statistical_analyses", "stage3", "models"),
                      pattern = "[.]rds$", full.names = TRUE)
  purrr::map_dfr(files, function(f) {
    obj <- readRDS(f)
    m <- obj$model
    if (length(m$smooth) < 2) return(NULL)
    w <- mgcv::concurvity(m, full = FALSE)$worst
    w <- w[rownames(w) != "para", colnames(w) != "para", drop = FALSE]
    idx <- which(upper.tri(w) & w > CONCURVITY_FLAG, arr.ind = TRUE)
    if (nrow(idx) == 0) return(NULL)
    tibble(species = obj$species, site = obj$site,
           term1 = strip_s(rownames(w)[idx[, 1]]),
           term2 = strip_s(colnames(w)[idx[, 2]]),
           worst = round(w[idx], 2))
  })
})
concurvity_phrase <- if (nrow(concurvity_pairs) == 0) "none" else
  with(order_rows(concurvity_pairs),
       and_list(paste0(term_text(term1), " and ", term_text(term2), " (",
                       SP_TEXT[species], " at ", site, "; ", f2(worst), ")")))

# ---- Acoustic occurrence (07_results_occurrence/occurrence_summary.R) ----------
occ <- rd("07_results_occurrence", "occurrence_by_site.csv")
occ_row <- function(sp, st) occ[occ$species == sp & occ$site == st, ]
occ_days <- function(sp, st) occ_row(sp, st)$detection_days
occ_pct <- function(sp, st) round(occ_row(sp, st)$pct_days_with_detections)
rec_days <- function(st) occ$recording_days[occ$site == st][1]
occ_txt <- function(sp, st) paste0(occ_days(sp, st), " days (", occ_pct(sp, st), "%)")
SP_OCC <- c(SP_TABLE, BW58 = "BW58 beaked whale")
occurrence_table <- function() {
  occ |>
    mutate(cell = ifelse(detection_days == 0, "0",
                         paste0(detection_days, " (", round(pct_days_with_detections), "%); ",
                                f1(presence_hours), " h"))) |>
    select(species, site, cell) |>
    tidyr::pivot_wider(names_from = site, values_from = cell) |>
    mutate(species = factor(species, names(SP_OCC))) |> arrange(species) |>
    transmute(Species = SP_OCC[as.character(species)],
              `EI (2014)` = EI, `KGI (2015–16)` = KGI, `CI (2016)` = CI)
}
# Sea-ice classes: bins with detections / bins in the class
ice_occ_txt <- function(sp, st, cls) with(ice_occ[ice_occ$species == sp & ice_occ$site == st &
                                                   ice_occ$ice_class == cls, ],
                                          paste0(n_detection_bins, " of ", n_bins))
ice_bins_class <- function(st, cls) ice_occ$n_bins[ice_occ$species == "Gm" & ice_occ$site == st &
                                                     ice_occ$ice_class == cls]

# ---- Sea-ice exploration (03_sea_ice/sea_ice_exploration.R; 15% threshold) ----
ice_ph <- rd("03_sea_ice", "ice_phenology_by_year.csv") |> filter(threshold_pct == 15)
ice_rank <- rd("03_sea_ice", "ice_recording_year_rank.csv") |> filter(threshold_pct == 15)
ice_anova <- rd("03_sea_ice", "ice_site_year_anova.csv")
ice_sync <- rd("03_sea_ice", "ice_site_synchrony.csv")
ice_phase <- rd("03_sea_ice", "ice_phase_occurrence.csv")
rec_ph <- function(st, col) ice_ph[[col]][ice_ph$site == st & ice_ph$recording]
ice_pctl <- function(st, m) ice_rank$percentile[ice_rank$site == st & ice_rank$metric == m]
ice_med <- function(st, m) ice_rank$median_other[ice_rank$site == st & ice_rank$metric == m]
ice_n_years <- function(st) ice_rank$n_years[ice_rank$site == st & ice_rank$metric == "ice_days"]
ice_anova_pct <- function(m, col) f1(ice_anova[[col]][ice_anova$metric == m])
ice_sync_r <- function(pair) f2(ice_sync$spearman_r[ice_sync$metric == "cumulative_ice_days" &
                                                      ice_sync$pair == pair])
phase_pct <- function(st, sp, ph) ice_phase$pct_days_detected[ice_phase$site == st &
                                                                ice_phase$species == sp &
                                                                ice_phase$phase == ph]
