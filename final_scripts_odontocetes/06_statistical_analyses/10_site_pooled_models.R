# =============================================================================
# Script  : 10_site_pooled_models.R
# Section : 2.6 Statistical analyses - shared vs site-specific responses
#           (PI feedback 2026-11-25; where it is reported is decided after the
#           results, decision N. Posdaljian)
# Purpose : Test whether a species' habitat response is shared among sites
#           (site/years) or differs between them, by pooling the binned data of
#           two or three sites into one hierarchical GAM (Pedersen et al. 2019,
#           PeerJ 7:e6876).
#           A. Term selection on the pooled data with a shared response:
#              presence ~ site + s(x1) + s(x2) + ...  (site = parametric
#              intercept per site), from the species' hypothesis set, same
#              backward selection, term cap (pooled bins) and VIF check as the
#              main models.
#           B. Each selected term refit three ways (other terms shared):
#                G  : s(x)                              shared curve
#                GI : s(x) + s(x, by = site, m = 1)     shared + site deviations
#                I  : s(x, by = site)                   separate curve per site
#              Verdict: "shared" if G is within 2 AIC of the best and the GI
#              deviations are not significant; "site-specific" if GI or I is
#              better than G by > 2 AIC and at least one deviation has
#              p <= 0.05; otherwise "unclear".
#           C. AR(1) check (glmmTMB, Tweedie, ar1(bin + 0 | site): one AR(1)
#              series per site; AR(1) used where it lowers AIC by > 2): LRT of
#              each shared term, and LRT of site-specific (I) vs shared (G)
#              curves for each term. glmmTMB cannot fit by-factor smooths, so
#              the site-specific curve for site s uses x_s = x at site s and a
#              constant elsewhere (the constant is absorbed by the site
#              intercepts).
#           D. Concurvity of each shared term with the rest of the model,
#              including the site term (raw FSLE/EKE differ in level mostly
#              between sites).
# Inputs  : stage1 binned_data/<species>_binned.csv, tables/model_inclusion.csv
#           (stage3), HYPOTHESIS_SETS (00_gam_parameters.R)
# Outputs : outputs/06_statistical_analyses/site_pooled/
#             pooled_selection.csv (selected terms per species x site set)
#             pooled_term_tests.csv (G / GI / I per term, AR(1) checks)
#             pooled_curves.csv, pooled_curves_<species>.png
# Author  : Natalie Posdaljian
# Date    : 2026-11-25
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Sites below the 10-detection-bin rule are not pooled (BW37 CI: 7
#           bins, too few for a site curve; decision 2026-11-25). Killer whales
#           are modeled at CI only, so no pooled model. "Site" = site/year.
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_parameters.R"))
source(file.path(FINAL_DIR, "06_statistical_analyses", "00_gam_functions.R"))
suppressPackageStartupMessages({ library(glmmTMB); library(patchwork) })

out_dir <- file.path(STATS_DIR, "site_pooled")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
DELTA_AIC <- 2

inclusion <- read_csv(file.path(STATS_DIR, "stage3", "tables", "model_inclusion.csv"),
                      show_col_types = FALSE) |>
  filter(status != "descriptive only")
modeled_sites <- split(inclusion$site, inclusion$species)
# Site sets: all modeled sites of the species and every pair of them
site_sets <- imap(modeled_sites, function(st, sp) {
  st <- SITE_ORDER[SITE_ORDER %in% st]
  if (length(st) < 2) return(NULL)
  sets <- combn(st, 2, simplify = FALSE)
  if (length(st) == 3) sets <- c(list(st), sets)
  map(sets, \(s) list(species = sp, sites = s))
}) |> compact() |> flatten()

# -----------------------------------------------------------------------------
# Model formulas
# -----------------------------------------------------------------------------
smooth_txt <- function(x, type) {
  cyc <- x %in% names(CYCLIC_VARS)
  bs <- if (cyc) ', bs = "cc"' else ""
  switch(type,
         G  = paste0("s(", x, bs, ", k = ", GAM_K, ")"),
         I  = paste0("s(", x, bs, ", by = site, k = ", GAM_K, ")"),
         GI = paste0("s(", x, bs, ", k = ", GAM_K, ") + s(", x, bs, ", by = site, k = ", GAM_K,
                     if (cyc) "" else ", m = 1", ")"))
}
fit_pooled <- function(dat, sp, terms, types = setNames(rep("G", length(terms)), terms)) {
  rhs <- c("site", map2_chr(terms, types[terms], smooth_txt))
  knots <- if (any(terms %in% names(CYCLIC_VARS)))
    CYCLIC_VARS[intersect(terms, names(CYCLIC_VARS))] else NULL
  gam(as.formula(paste(sp, "~", paste(rhs, collapse = " + ")), env = globalenv()),
      data = dat, family = gam_family(), method = GAM_METHOD, select = TRUE, knots = knots)
}
term_p <- function(m) setNames(summary(m)$s.pv, rownames(summary(m)$s.table))

# A. backward selection with shared curves (site intercept always kept)
backward_pooled <- function(dat, sp, preds, cap) {
  repeat {
    m <- fit_pooled(dat, sp, preds)
    pv <- setNames(summary(m)$s.pv, term_variable(rownames(summary(m)$s.table)))
    worst <- names(pv)[which.max(pv)]
    if (length(preds) == 1 || (max(pv) <= P_THRESH && length(preds) <= cap)) break
    preds <- setdiff(preds, worst)
  }
  list(model = m, terms = preds, pv = pv)
}
# VIF of each term with the site factor in the model (GVIF; for 1-df terms
# GVIF = VIF, comparable with VIF_THRESH)
vif_site <- function(dat, sp, terms) {
  if (length(terms) == 0) return(NULL)
  v <- car::vif(glm(reformulate(c("site", terms), sp), data = dat))
  if (is.matrix(v)) v <- setNames(v[, "GVIF"], rownames(v))   # 2-site sets give a vector
  v[terms]
}
# Backward selection, then the VIF check on the final model as in the main
# models: if a term has VIF > 5 (with site), the least significant of those
# terms is removed from the candidates and selection is repeated.
select_pooled <- function(dat, sp, preds, cap) {
  preds <- remove_zero_var(dat, preds)
  dropped <- character(0)
  repeat {
    be <- backward_pooled(dat, sp, preds, cap)
    vif <- vif_site(dat, sp, be$terms)
    if (length(be$terms) < 1 || max(vif) <= VIF_THRESH) break
    high <- names(vif)[vif > VIF_THRESH]
    drop <- high[which.max(be$pv[high])]
    dropped <- c(dropped, paste0(drop, " (VIF ", round(vif[[drop]], 1), ")"))
    preds <- setdiff(preds, drop)
  }
  list(model = be$model, terms = be$terms, no_sig = max(be$pv) > P_THRESH, vif = vif,
       vif_dropped = paste(dropped, collapse = "; "))
}

# -----------------------------------------------------------------------------
# C. glmmTMB AR(1) helpers
# -----------------------------------------------------------------------------
tmb_data <- function(dat, terms) {
  d <- dat |> mutate(bin = factor(bin_index, levels = 0:max(bin_index)))
  for (cv in intersect(terms, names(CYCLIC_VARS))) {
    d[[paste0(cv, "_sin2")]] <- sin(2 * d[[cv]] * pi / 180)
    d[[paste0(cv, "_cos2")]] <- cos(2 * d[[cv]] * pi / 180)
  }
  # site-specific copies: x at site s, constant (median at s) elsewhere
  for (tm in setdiff(terms, names(CYCLIC_VARS))) for (s in levels(d$site)) {
    med <- median(d[[tm]][d$site == s])
    d[[paste0(tm, "__", s)]] <- if_else(d$site == s, d[[tm]], med)
  }
  d
}
tmb_rhs <- function(tm, type, sites) {
  if (tm %in% names(CYCLIC_VARS)) {
    if (type == "G") return(paste0(tm, "_sin2 + ", tm, "_cos2"))
    return(paste0("site:", tm, "_sin2 + site:", tm, "_cos2"))
  }
  if (type == "G") paste0("s(", tm, ", k = ", GAM_K, ")")
  else paste(paste0("s(", tm, "__", sites, ", k = ", GAM_K, ")"), collapse = " + ")
}
fit_tmb <- function(d, sp, terms, types, ar1) {
  sites <- levels(d$site)
  rhs <- c("site", if (length(terms)) map2_chr(terms, types[terms], tmb_rhs, sites = sites))
  if (ar1) rhs <- c(rhs, "ar1(bin + 0 | site)")
  tryCatch(suppressWarnings(glmmTMB(as.formula(paste(sp, "~", paste(rhs, collapse = " + "))),
                                    data = d, family = tweedie(link = "log"))),
           error = function(e) e)
}
ok <- function(f) !inherits(f, "error") && is.finite(suppressWarnings(AIC(f)))
lrt <- function(full, red) {
  if (!ok(full) || !ok(red)) return(NA_real_)
  df <- attr(logLik(full), "df") - attr(logLik(red), "df")
  pchisq(2 * (as.numeric(logLik(full)) - as.numeric(logLik(red))), df = max(df, 1),
         lower.tail = FALSE)
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
sel_rows <- list(); test_rows <- list(); curve_rows <- list()
for (ss in site_sets) {
  sp <- ss$species; sites <- ss$sites; set_id <- paste(sites, collapse = "+")
  dat <- read_csv(file.path(STATS_DIR, "stage1", "binned_data", paste0(sp, "_binned.csv")),
                  show_col_types = FALSE) |>
    filter(site %in% sites, !ice_excluded) |>
    mutate(site = factor(site, levels = sites))
  n_det <- sum(dat[[sp]] > 0)
  cap <- term_cap(nrow(dat), n_det <= LOW_SAMPLE_MAX_BINS)
  sel <- select_pooled(dat, sp, HYPOTHESIS_SETS[[sp]], cap)
  terms <- if (sel$no_sig) character(0) else sel$terms
  conc <- if (length(terms)) {
    cc <- mgcv::concurvity(sel$model, full = TRUE)["worst", ]
    cc[names(cc) != "para"]
  } else NULL
  sel_rows[[paste(sp, set_id)]] <- tibble(
    species = sp, sites = set_id, n_bins = nrow(dat), n_detection_bins = n_det,
    term_cap = cap, terms = if (length(terms)) paste(terms, collapse = " + ") else "(none)",
    dev_expl_pct = round(100 * summary(sel$model)$dev.expl, 1),
    max_vif_with_site = if (is.null(sel$vif)) NA_real_ else round(max(sel$vif), 2),
    vif_dropped = sel$vif_dropped,
    max_concurvity = if (is.null(conc)) NA_real_ else round(max(conc), 2))
  message("\n", sp, " [", set_id, "]: ", sel_rows[[paste(sp, set_id)]]$terms,
          " (", sel_rows[[paste(sp, set_id)]]$dev_expl_pct, "%)")
  if (!length(terms)) next

  base_types <- setNames(rep("G", length(terms)), terms)
  # AR(1): is it supported for the shared model?
  d_tmb <- tmb_data(dat, terms)
  g0 <- fit_tmb(d_tmb, sp, terms, base_types, FALSE)
  g1 <- fit_tmb(d_tmb, sp, terms, base_types, TRUE)
  use_ar <- ok(g0) && ok(g1) && AIC(g1) - AIC(g0) < -DELTA_AIC
  g_ref <- if (use_ar) g1 else g0

  for (tm in terms) {
    fits <- map(c(G = "G", GI = "GI", I = "I"), function(tp) {
      ty <- base_types; ty[tm] <- tp
      tryCatch(fit_pooled(dat, sp, terms, ty), error = function(e) NULL)
    })
    aic <- map_dbl(fits, \(f) if (is.null(f)) NA_real_ else AIC(f))
    p_g <- term_p(fits$G)[[paste0("s(", tm, ")")]]
    p_dev <- if (!is.null(fits$GI)) {
      pp <- term_p(fits$GI); pp[grepl(":site", names(pp)) & grepl(tm, names(pp), fixed = TRUE)]
    } else NA
    best <- min(aic, na.rm = TRUE)
    verdict <- case_when(
      aic[["G"]] - best <= DELTA_AIC & all(p_dev > P_THRESH, na.rm = TRUE) ~ "shared",
      aic[["G"]] - best > DELTA_AIC & any(p_dev <= P_THRESH, na.rm = TRUE) ~ "site-specific",
      TRUE ~ "unclear")
    # AR(1) checks: shared term (drop it) and site-specific vs shared curve
    ty_i <- base_types; ty_i[tm] <- "I"
    red <- fit_tmb(d_tmb, sp, setdiff(terms, tm), base_types, use_ar)
    i_ref <- fit_tmb(d_tmb, sp, terms, ty_i, use_ar)
    test_rows[[paste(sp, set_id, tm)]] <- tibble(
      species = sp, sites = set_id, term = tm,
      aic_G = round(aic[["G"]], 1), aic_GI = round(aic[["GI"]], 1), aic_I = round(aic[["I"]], 1),
      dAIC_GI_vs_G = round(aic[["GI"]] - aic[["G"]], 1),
      dAIC_I_vs_G = round(aic[["I"]] - aic[["G"]], 1),
      p_shared = signif(p_g, 3),
      p_site_deviations = paste(signif(p_dev, 2), collapse = "; "),
      verdict = verdict,
      ar1_used = use_ar,
      p_shared_ar1 = signif(lrt(g_ref, red), 3),
      p_site_specific_vs_shared_ar1 = signif(lrt(i_ref, g_ref), 3),
      concurvity = if (!is.null(conc)) round(conc[[paste0("s(", tm, ")")]], 2) else NA_real_)
    message("  ", tm, ": dAIC GI ", round(aic[["GI"]] - aic[["G"]], 1), ", I ",
            round(aic[["I"]] - aic[["G"]], 1), " -> ", verdict)

    # Curves for the all-sites (or only) set: shared (G) and per-site (I)
    if (length(sites) == length(modeled_sites[[sp]])) {
      nd_base <- dat |> summarise(across(all_of(terms), median))
      # shared curve once over the pooled range; site curves over each site's range
      combos <- c(list(c("G", sites[1])), map(sites, \(s) c("I", s)))
      for (cb in combos) {
        tp <- cb[1]; s <- cb[2]
        rng <- if (tm %in% names(CYCLIC_VARS)) CYCLIC_VARS[[tm]] else
          if (tp == "G") range(dat[[tm]]) else range(dat[[tm]][dat$site == s])
        nd <- nd_base[rep(1, 100), ] |> mutate(site = factor(s, levels = sites))
        nd[[tm]] <- seq(rng[1], rng[2], length.out = 100)
        pr <- predict(fits[[tp]], newdata = nd, type = "terms", se.fit = TRUE)
        col <- if (tp == "G") paste0("s(", tm, ")") else paste0("s(", tm, "):site", s)
        curve_rows[[length(curve_rows) + 1]] <- tibble(
          species = sp, sites = set_id, term = tm, model = tp, site = s, x = nd[[tm]],
          fit = pr$fit[, col], se = pr$se.fit[, col], verdict = verdict)
      }
    }
  }
}
write_csv(bind_rows(sel_rows), file.path(out_dir, "pooled_selection.csv"))
tests <- bind_rows(test_rows)
write_csv(tests, file.path(out_dir, "pooled_term_tests.csv"))
curves <- bind_rows(curve_rows)
write_csv(curves, file.path(out_dir, "pooled_curves.csv"))

# -----------------------------------------------------------------------------
# Figures: per species, one panel per term; black = shared curve (G) over the
# pooled range, colors = site-specific curves (I) over each site's range
# -----------------------------------------------------------------------------
for (sp in unique(curves$species)) {
  cs <- filter(curves, species == sp)
  p <- ggplot(cs, aes(x, fit)) +
    geom_ribbon(data = filter(cs, model == "I"),
                aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se, fill = site), alpha = 0.15) +
    geom_line(data = filter(cs, model == "I"), aes(colour = site), linewidth = 0.7) +
    geom_line(data = filter(cs, model == "G"),
              colour = "black", linewidth = 0.8, linetype = "dashed") +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey50") +
    scale_colour_manual(values = SITE_COLORS) + scale_fill_manual(values = SITE_COLORS) +
    facet_wrap(~ paste0(variable_label(term), "\n(", verdict, ")"), scales = "free",
               labeller = label_wrap_gen(35)) +
    labs(x = NULL, y = "Partial effect (log scale)", colour = "Site", fill = "Site",
         title = paste0(species_label(sp), " (", cs$sites[1], "): dashed = shared curve; ",
                        "colors = site-specific curves")) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, paste0("pooled_curves_", sp, ".png")), p,
         width = 7.5, height = 3.2 + 2.2 * ceiling(n_distinct(cs$term) / 3), dpi = 300)
}

print(bind_rows(sel_rows), n = Inf, width = Inf)
print(select(tests, species, sites, term, dAIC_GI_vs_G, dAIC_I_vs_G, p_site_deviations,
             verdict, ar1_used, p_shared_ar1, p_site_specific_vs_shared_ar1, concurvity),
      n = Inf, width = Inf)
message("Site-pooled models complete. Outputs in ", out_dir)
