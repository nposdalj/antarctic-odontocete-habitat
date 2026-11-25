# Handoff — odontocete habitat manuscript (status 2026-11-24)

Read this first when resuming work (e.g., in Claude Code on another computer).
It summarizes every decision made in the cleanup/analysis sessions, how to run
the pipeline, the current results, and what is still open.

---

## 1. How to run things

- **Working directory:** run every R script from the **repository root**
  (`antarctic-odontocete-habitat/`), e.g.
  `Rscript final_scripts_odontocetes/06_statistical_analyses/03_stage3_harmonized_models.R`.
  Scripts use `BASE_DIR <- "."`; running from inside `final_scripts_odontocetes/`
  creates a stray nested `final_scripts_odontocetes/outputs/` folder.
- **R:** 4.5.1. Packages used: tidyverse, mgcv, gratia, car, glmmTMB, terra,
  ncdf4, patchwork, ragg, jpeg, png, knitr, kableExtra.
- **Credentials** (Copernicus/AVISO downloads only): `Documents/.Renviron`; run
  with `R_ENVIRON_USER=".../.Renviron"`. Never print or commit credentials.
- **Manuscript:** Quarto, `final_scripts_odontocetes/manuscript/manuscript.qmd`
  + `_values.R` (all numbers/tables read from `outputs/`). Render with the Quarto
  bundled with RStudio, with `QUARTO_R` pointing at R 4.5.1:
  `quarto render manuscript.qmd` (PDF via xelatex/MiKTeX; MiKTeX auto-install on).
  `main.tex`, `main.pdf`, `figures/`, `overleaf_upload.zip` in `manuscript/` are
  an obsolete LaTeX version (safe to delete after checking).
- **Rule:** any change to figures/text must also be made in `manuscript/`.

### Pipeline order (no downloads)
1. `03_sea_ice/bremen_01_download.R` (AMSR2 2013–2021; only if data missing)
2. `03_sea_ice/amsr2_site_series.R` → `outputs/03_sea_ice/amsr2_daily_site.csv`
3. `05_environmental_data/04b_eke_per_cell.R` (per-cell EKE; only if velocities change)
4. `05_environmental_data/05_compile_dataset.R`
5. `04_acoustic_detection/build_presence.R`, `duty_cycle_simulation.R` (only if detections change)
6. `06_statistical_analyses/01_stage1_automated_selection.R` (binning, ice cutoffs, inclusion)
7. `06_statistical_analyses/03_stage3_harmonized_models.R` (final habitat models)
8. `06_statistical_analyses/04_sea_ice_models.R`
9. `06_statistical_analyses/05_ar1_sensitivity.R stage3` and `... sea_ice`
10. `07_results_occurrence/occurrence_summary.R`, `occurrence_figures.R`, `rare_detections.R`
11. `08_results_gams/results_tables.R`, `response_curves.R`, `sea_ice_occurrence.R`, `ice_edge_detections.R`
12. `02_deep_oceanography/deep_oceanography.R`
13. `03_sea_ice/sea_ice_exploration.R`, `bremen_02_compare_maps.R` (needs `modis_01_fetch.R` images)
14. `05_environmental_data/07_site_year_comparison.R`, `08_site_year_tests.R`
15. Optional checks: `06_statistical_analyses/07_permutation_hypotheses.R C` (~40 min),
    `08_daily_ar1_check.R`
16. Render the manuscript.

---

## 2. Decisions (all by N. Posdaljian)

### Data and predictors
- Sites: EI (2014), KGI (2015–16), CI (2016); ±20 km box (justified by detection ranges + BGC grid).
- **Site and year are confounded**: always write "site/year differences".
- Water-column variables averaged in 3 layers: surface 0–100, mid 100–500, deep 500–1000 m.
- Water mass = **temperature only** (validated against WOD profiles; salinity/O2 poorly reproduced).
- Predictors are anomalies from a 10-year climatology except EKE, FSLE, FSLV orientation.
- **EKE = per-cell definition** (u′, v′ per grid cell vs that cell's 10-yr mean, then box/layer mean; Ducet et al. 2000 [CONFIRM]). EKE comparison NOT reported.
- **Chl-a REMOVED completely (2026-11-25, after PI feedback)** and not mentioned anywhere in the paper or supplement. Reason: no mechanism for a specific lag; the 3-month lag had been chosen from the same data. Unreported check `06_statistical_analyses/09_chla_sensitivity.R`: adding lag-0 or lag-3 chl-a changed 3 / 5 of 11 models, never the same ones, rarely robust to AR(1); prewhitened CCF (lags 0–180 d) had no lag outside the 95% bounds.
- **Pilot whales = mid-water forager** (2026-11-25, PI feedback): config depths 0/16/375 m; candidates surface + mid temperature, mid EKE (tag studies: Baird 2002, Heide-Jørgensen 2002, Visser 2017).
- Species-depth and predictor rationale with literature: Methods 2.6 + Table 3; references in `manuscript/references_pending.md`.
- **Sea ice = AMSR2 ASI (Univ. Bremen, 6.25 km)** everywhere. Copernicus ice eliminated, with **no comparison reported**. Context record 2013–2021; gaps ≤5 days interpolated in the exploration.
- Daily AMSR2 metrics: box-mean concentration, open-water fraction, signed distance to ice edge (15%; coast buffer 12.5 km, isolated pixels ignored).

### Habitat models
- Tweedie GAMs, REML, k = 4, `select = TRUE` (shrinkage; edf ≈ 1 = effectively linear, Marra & Wood 2011), 4-day bins (killer whales 2-day).
- Inclusion: ≥10 detection bins; 10–20 = low sample (≤2 terms); term cap floor(n_bins/10).
- Ice filter: bins with ice ≤ the species' max ice at detection (AMSR2: Pm 94.7%, Oo 96.6%, Gm 17.8%, BW29 0%, BW37 0.1%).
- **Candidate sets = a priori hypothesis sets** (`HYPOTHESIS_SETS` in `00_gam_parameters.R`): all species FSLE, FSLV, MLD anomaly, chl-a; + temperature/EKE in the species' layers; + surface temperature for pilot and killer whales. (Stage 1 union sets abandoned.)
- **VIF used as a check on the final model, not a filter** (threshold 5; `select_vif_checked()`); highest final VIF 1.6.
- Selection: backward elimination (p ≤ 0.05) with term cap.
- **Reporting decision (agreed 2026-11-23):** report the models as exploratory, with AR(1) robustness markers; permutation test at most one supplementary sentence. No pooling of sites (responses differ among sites).
- **Site as a variable (2026-11-25):** `06_statistical_analyses/10_site_pooled_models.R` (shared smooth + site intercept vs site-specific smooths; all sites and pairs). Mostly site-specific or unclear; only BW29 ~ MLD shared in every combination, and not AR(1)-robust. Decision: one Discussion sentence, not in the supplement.
- **Site/year tests (2026-11-25, PI request):** `05_environmental_data/08_site_year_tests.R`. Annual means 2011–2019 with years as replicates (two-way ANOVA, year block + site; paired t-tests, Holm); recording-window ranks; EI 2014 temperature figure. Sites differ in every variable except surface temperature; EI ≈ KGI except EKE/FSLE; CI distinct. EI 2014 cold mid/deep anomaly from late Feb 2014 (similar in 2009–2010 at EI); not confirmed in situ, HYCOM does not reproduce the deep part. Supplement S7 rewritten (old variance-partition table dropped; monthly figure kept).

### Sea-ice models (agreed 2026-11-23)
- Concentration only (fit better than distance to ice edge by AIC everywhere).
- Where all detections were at 0% ice (Gm KGI, BW29 KGI & CI, BW37 KGI): **no fitted curve**; report "never detected in a bin with sea ice".
- For those, open-water bins only: distance-to-ice-edge models (none significant).
- EI: no ice variation during recording (recorded only before ice advance).

### Manuscript decisions (2026-11-25, after first full read)
- Species-depth bullets deleted (Table 3 covers them); predictor-rationale bullets condensed to one paragraph.
- Ice-class occurrence figure moved to Supplement S4; ice-edge detections figure stays in the main text.
- Figure 3 = weekly mean presence minutes per recorded day (the model response units), not days per week.
- Response stays minutes/day (hourly/daily 0/1 tested 2026-11-19, rejected).
- FSLV and mid-water variables stay as candidates (a priori sets; rarely selected is a result).
- Target: npj Acoustics collection https://www.nature.com/collections/gfdagjjhge (Natalie 2026-11-25). Guidelines pages require sign-in for automated fetches; word/figure limits not yet confirmed.
- Permutation test C rerun with final hypothesis sets; supplement figure from `06_statistical_analyses/11_permutation_figure.R`.

### Figures
- Response curves: kernel densities at panel bottom (all bins grey, detection bins site colour), clipped to data range. Main figure = models with ≥2 terms; single-term models in supplement.
- BW58 merged into the main weekly-presence figure.
- Supplement: AMSR2 + MODIS maps (4 clear CI dates), combined ice climatology + detections figure, year×day heatmap, phenology metrics.

---

## 3. Current results (AMSR2, hypothesis sets, no chl-a, Gm mid-water; 2026-11-25)

| Model | Final terms | Dev. % |
|---|---|---|
| Pm EI | FSLE + MLD + deep temperature | 75.6 |
| Pm CI | FSLE + deep temperature + mid-water EKE | 48.9 |
| Gm EI | surface temperature + mid-water temperature | 77.2 |
| Gm KGI | MLD + surface temperature | 55.7 |
| Gm CI | MLD + surface temperature | 53.4 |
| Oo CI | surface temperature + mid-water EKE | 26.7 |
| BW29 EI | deep EKE | 29.9 |
| BW29 KGI | MLD + deep temperature | 21.3 |
| BW29 CI | MLD + deep temperature | 42.7 |
| BW37 EI | FSLE + FSLV | 47.8 |
| BW37 KGI | MLD + deep temperature | 43.3 |

- Highest final VIF 1.26. AR(1)-robust terms (binned): Pm EI FSLE, MLD; Pm CI deep T; Gm all except CI MLD; BW37 EI both, KGI deep T; none for Oo CI or BW29.
- Daily AR(1) check (`08_daily_ar1_check.R`) agrees with binned check, slightly stricter → daily GAMM NOT adopted.
- MLD + deep temperature recurs in BW29 KGI/CI and BW37 KGI (candidate shared response for the site-as-variable step).
- Permutation test (time-shift null, test C): only Gm EI p = 0.04 unadjusted; none after BH correction. Flaws: low power, conservative for stable predictors, overlapping null samples.
- Sea-ice story (descriptive): CI has a long heavy ice season, EI/KGI short and patchy; KGI 2015 = heaviest/longest ice year in its record, CI 2016 one of the lightest, EI 2014 essentially no ice season (recording ended before ice). Beaked whales never detected with any ice; sperm whales at CI detected up to 77 km inside the ice edge; killer whales mostly at advance/retreat (CI) and once 35 km inside 96% pack (KGI). MODIS shows leads/broken pack on heavy-ice days.

---

## 4. Open tasks (in order)

1. **Next:** sea-ice detections vs distance to the ice edge (Results/supplement) + leads in the Discussion; then rerun the permutation test once the models are final. (Site as a variable and site/year tests: done.)
2. **Manuscript pass** (`manuscript.qmd`):
   - Replace all Copernicus ice text with AMSR2 (Spreen et al. 2008 [CONFIRM]); add ice-edge distance metric and coast buffer to Methods 2.3.
   - Methods 2.6: hypothesis sets (Table 2 with references), VIF-check rule, shrinkage/edf sentence, sea-ice model changes, exploratory framing; optional one-sentence permutation sensitivity.
   - Results: occurrence-led; sea-ice section with phenology/recording-year context and "never detected with ice"; habitat section exploratory with AR(1) markers.
   - Update figures (combined ice figure, maps, curves) and supplement numbering; check stale references (e.g., Supplement S3 EKE comparison was removed).
   - Discussion notes: site/year confound, sea-ice interpretation, concurvity, detection vs ambient noise, raw-vs-anomaly validation.
3. **Cleanup:** delete exploratory outputs (sensitivity_* folders, `_previous_run_*`, box40 comparison, `06_permutation_test.R` Stage-1 version, tests A/B files) — ask first. Update `INVENTORY.md`/READMEs.
4. **Placeholders still open:** detection ranges (killer/pilot whale; sperm whale 20 km ref), killer whale LTSA method, site characterization, species-depth literature, product DOIs, BW29 species attribution (southern bottlenose whale vs *Mesoplodon* sp.? — user's old table used *Mesoplodon* sp.), hypothesis-set citations.
5. References → BibTeX at the end.
6. **After the analysis is final:** update Natalie's existing summary artifact (made in another Claude Code session): https://claude.ai/artifact/G3oyxeQrGcrpDWAv3o9hA9 — read it first, then republish to the same URL.
