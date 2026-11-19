# =============================================================================
# Script  : 05_ar1_sensitivity.R
# Section : 2.6 Statistical analyses (sensitivity analysis; supplement)
# Purpose : Check whether the significant terms of the final habitat models
#           stay significant when residual autocorrelation is modeled
#           explicitly. Each final model is refit with glmmTMB, with and without
#           an AR(1) process; AR(1) is used where it lowers AIC by > 2. Same
#           smooth terms (s(x, k = GAM_K), fitted as penalized random-effect
#           smooths), a Tweedie family with log link, and an AR(1) process
#           on the bin index. Each term is tested with a likelihood-ratio test
#           (model with vs without the term, both with AR(1)), and the result
#           is compared with the term's p-value in the original GAM.
# Inputs  : <OUTPUT_ROOT>/06_statistical_analyses/<stage>/models/*.rds
#             (final models from 01, 03 or 04; each holds the fitted GAM, the
#             data it was fit to, and its diagnostics)
#           <stage> = first command-line argument; default "stage3" if it
#           exists, otherwise "stage1"
# Outputs : <OUTPUT_ROOT>/06_statistical_analyses/<stage>/tables/
#             ar1_sensitivity_terms.csv  (species, site, term, edf and p in the
#               GAM, likelihood-ratio p with AR(1), significant in each)
#             ar1_sensitivity_models.csv (species, site, Ljung-Box p and
#               lag-1 residual ACF of the GAM, AR(1) correlation, number of
#               terms significant in both fits, terms lost, fit status)
# Author  : Natalie Posdaljian
# Date    : 2026-11-19
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : - glmmTMB (>= 1.1.10) fits s() smooths by maximum likelihood
#             (Laplace approximation), which is more stable for these sample
#             sizes (34-151 bins) than mgcv::gamm, whose penalized quasi-
#             likelihood fits of Tweedie models failed or changed with the
#             optimizer settings (tested 2026-11-19).
#           - AR(1) term: ar1(bin + 0 | series), with `bin` a factor over ALL
#             bin indices of the series, so bins removed by the ice filter are
#             gaps in the latent AR(1) process, not neighbours.
#           - Likelihood-ratio df = difference in the number of estimated
#             parameters (for a smooth: its unpenalized linear part plus its
#             smoothing variance). Tests on a variance at its boundary are
#             conservative.
#           - Cyclic smooths (fsle_orient, bs = "cc") are not supported by
#             glmmTMB's s(); orientation enters as sin(2 theta) + cos(2 theta)
#             (periodic over 180 degrees) and is tested as a pair.
#           - Run from the repository root:
#               Rscript final_scripts_odontocetes/06_statistical_analyses/05_ar1_sensitivity.R [stage]
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))

suppressPackageStartupMessages(library(glmmTMB))

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else
  if (dir.exists(file.path(STATS_DIR, "stage3", "models"))) "stage3" else "stage1"
message("AR(1) sensitivity for ", stage)

models    <- read_stage_models(stage)
table_dir <- stage_dir(stage, "tables")

# -----------------------------------------------------------------------------
# 1. Refit each final model with glmmTMB, without and with AR(1)
# -----------------------------------------------------------------------------
# AR(1) is judged supported when it lowers AIC by more than AR1_DELTA_AIC.
# Where it is not supported, its correlation is poorly identified (estimates
# near 1 with no AIC gain were seen in testing) and the model without AR(1)
# is the reference.
AR1_DELTA_AIC <- 2

# Cyclic (axial) terms such as fsle_orient cannot be cyclic smooths in
# glmmTMB; they enter as the pair sin(2 theta) + cos(2 theta) (columns
# <var>_sin2 / <var>_cos2 added in refit_ar1()), which is periodic over 180
# degrees, and are tested as a pair.
term_rhs <- function(tm) {
  if (tm %in% names(CYCLIC_VARS)) paste0(tm, "_sin2 + ", tm, "_cos2")
  else paste0("s(", tm, ", k = ", GAM_K, ")")
}

fit_tmb <- function(response, terms, data, ar1) {
  rhs <- if (length(terms) > 0) vapply(terms, term_rhs, character(1)) else "1"
  if (ar1) rhs <- c(rhs, "ar1(bin + 0 | series)")
  fml <- as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
  tryCatch(
    suppressWarnings(glmmTMB(fml, data = data, family = tweedie(link = "log"))),
    error = function(e) e)
}

# Likelihood-ratio p-value of each term: full model vs model without it
lrt_p <- function(full, response, terms, data, ar1) {
  map_dbl(terms, function(tm) {
    red <- fit_tmb(response, setdiff(terms, tm), data, ar1)
    if (inherits(red, "error") || !is.finite(as.numeric(logLik(red)))) return(NA_real_)
    df <- attr(logLik(full), "df") - attr(logLik(red), "df")
    pchisq(2 * (as.numeric(logLik(full)) - as.numeric(logLik(red))), df = df,
           lower.tail = FALSE)
  })
}

# AR(1) correlation of the fitted latent process
ar1_phi <- function(fit) {
  vc <- VarCorr(fit)$cond
  attr(vc[[grep("series", names(vc))[1]]], "correlation")[1, 2]
}

refit_ar1 <- function(m) {
  gam_fit <- m$model
  terms   <- term_variable(rownames(summary(gam_fit)$s.table))
  if (length(terms) == 0) return(NULL)        # intercept-only model
  response <- all.vars(formula(gam_fit))[1]

  gam_tbl <- as_tibble(summary(gam_fit)$s.table, rownames = "term") |>
    transmute(term, edf_gam = edf, p_gam = `p-value`)
  failed <- function(msg) list(terms = mutate(gam_tbl, p_tmb = NA_real_, p_ar1 = NA_real_),
                               phi = NA_real_, delta_aic = NA_real_, status = msg)
  data <- m$data |>
    mutate(bin = factor(bin_index, levels = 0:max(bin_index)), series = factor(1))
  for (cv in intersect(terms, names(CYCLIC_VARS))) {
    data[[paste0(cv, "_sin2")]] <- sin(2 * data[[cv]] * pi / 180)
    data[[paste0(cv, "_cos2")]] <- cos(2 * data[[cv]] * pi / 180)
  }
  no_ar <- fit_tmb(response, terms, data, ar1 = FALSE)
  ar    <- fit_tmb(response, terms, data, ar1 = TRUE)
  if (inherits(no_ar, "error")) return(failed(paste("failed:", conditionMessage(no_ar))))
  ar_ok <- !inherits(ar, "error") && is.finite(AIC(ar))

  list(terms = mutate(gam_tbl,
                      p_tmb = lrt_p(no_ar, response, terms, data, ar1 = FALSE),
                      p_ar1 = if (ar_ok) lrt_p(ar, response, terms, data, ar1 = TRUE)
                              else NA_real_),
       phi = if (ar_ok) ar1_phi(ar) else NA_real_,
       delta_aic = if (ar_ok) AIC(ar) - AIC(no_ar) else NA_real_,
       status = if (ar_ok) "ok" else "AR(1) fit failed")
}

results <- imap(models, function(m, id) {
  message("  ", id)
  r <- refit_ar1(m)
  if (is.null(r)) return(NULL)
  list(
    terms = mutate(r$terms, species = m$species, site = m$site, .before = 1),
    model = tibble(species = m$species, site = m$site,
                   ljung_box_p_gam = m$summary$ljung_box_p %||% NA_real_,
                   resid_acf_lag1_gam = m$summary$resid_acf_lag1,
                   delta_aic_ar1 = r$delta_aic, ar1_phi = r$phi, status = r$status)
  )
}) |>
  compact()

# -----------------------------------------------------------------------------
# 2. Summaries
# -----------------------------------------------------------------------------
# Reference fit per model: with AR(1) where it is supported, otherwise
# without. A term "holds" when it is significant in the reference fit.
model_tbl <- map(results, "model") |>
  list_rbind() |>
  mutate(ar1_supported = !is.na(delta_aic_ar1) & delta_aic_ar1 < -AR1_DELTA_AIC)

term_tbl <- map(results, "terms") |>
  list_rbind() |>
  left_join(select(model_tbl, species, site, ar1_supported), by = c("species", "site")) |>
  mutate(p_reference = if_else(ar1_supported, p_ar1, p_tmb),
         holds = p_reference <= P_THRESH,
         edf_gam = round(edf_gam, 2),
         across(c(p_gam, p_tmb, p_ar1, p_reference), ~ signif(.x, 3)))

model_tbl <- model_tbl |>
  left_join(term_tbl |>
              group_by(species, site) |>
              summarise(n_terms = n(),
                        n_hold = sum(holds, na.rm = TRUE),
                        terms_not_holding = paste(term[!holds & !is.na(holds)],
                                                  collapse = "; "),
                        .groups = "drop"),
            by = c("species", "site")) |>
  mutate(across(c(ljung_box_p_gam, resid_acf_lag1_gam, delta_aic_ar1, ar1_phi),
                ~ round(.x, 3)))

write_csv(term_tbl,  file.path(table_dir, "ar1_sensitivity_terms.csv"))
write_csv(model_tbl, file.path(table_dir, "ar1_sensitivity_models.csv"))
print(model_tbl, n = Inf, width = Inf)
