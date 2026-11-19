# 06_statistical_analyses — Tweedie GAMs (Section 2.6)

Two model types per species × site (species in `MODELED_SPECIES`, `config.R`):

* **Sea-ice model**: presence ~ s(ice_conc, k = 4), fitted to all bins (no ice
  filter). Main-text results.
* **Habitat model**: fitted to ice-filtered bins, without sea ice as a candidate.
  Built in three stages: Stage 1 automated (supplement) → Stage 2 decision table →
  Stage 3 harmonized (main text).

Both use the same species bin length (residual ACF, 4-day cap), Tweedie family
(log link, p in [1.1, 1.9]), REML, `select = TRUE`, k = 4, and inclusion rule
(≥ 10 bins with detections; 10–20 = low sample). FSLE orientation (axial angle,
−90° ≡ 90°) always enters as a cyclic smooth `s(fsle_orient, bs = "cc", k = 4)`
with knots at ±90°. Results are reported by
`08_results_gams/` (Section 3.2).

## Files

| File | Role |
|---|---|
| `00_gam_parameters.R` | Every analysis setting (family, k, VIF and p thresholds, bin-length cap, inclusion rule, term caps, climatology requirements, raw/anomaly rule, candidate variables, themes, diagnostic and shape-class rules). Sourced by all scripts below. |
| `00_gam_functions.R` | Shared functions: deseasoning, raw/anomaly choice, residual-ACF bin length, binning, ice cutoff, VIF screening, theme selection, backward elimination with term caps, diagnostics, shape classes, and the common result format. |
| `01_stage1_automated_selection.R` | **Stage 1** (habitat, supplement): 10-yr climatology; per site × variable raw OR anomaly by seasonal variance share → ACF bin length → bins → ice filter → inclusion rule → zero variance → VIF ≤ 5 → theme winners by single-term AIC → backward elimination (p ≤ 0.05; at most floor(n_bins / 10) terms, 2 in low-sample models; a last remaining term with p > 0.05 = "no significant terms"). Saves all intermediate tables, binned data (used by Stages 3 and 4), final models (.rds) and diagnostics. |
| `02_stage2_decision_table.R` | **Stage 2**: decision table per species × site × candidate variable (raw/anomaly version, screening, single-term p / deviance explained / AIC rank, theme winner, final p / edf / per-term deviance, direction, shape class), theme and variable summaries, selection heatmap, single-term curves of all candidates. |
| `03_stage3_harmonized_models.R` | **Stage 3 (template)**, habitat, main text: fits the manually harmonized term sets — `[INSERT harmonized formulas after reviewing Stage 2]` — on the Stage 1 binned data, then backward-eliminates within them with the Stage 1 rules (no non-significant terms in final models). Writes the Stage 1 result format. |
| `04_sea_ice_models.R` | **Sea-ice models**, main text: presence ~ s(ice_conc) on all bins; replaces `GAMs/seasonalModels/` (concentration only). |

## Run order (from the repository root)

```
Rscript final_scripts_odontocetes/06_statistical_analyses/01_stage1_automated_selection.R
Rscript final_scripts_odontocetes/06_statistical_analyses/02_stage2_decision_table.R
#   -> review stage2/stage2_decision_table.csv and figures, fill in harmonized_terms
Rscript final_scripts_odontocetes/06_statistical_analyses/03_stage3_harmonized_models.R
Rscript final_scripts_odontocetes/06_statistical_analyses/04_sea_ice_models.R   # needs 01 only
Rscript final_scripts_odontocetes/08_results_gams/results_tables.R
Rscript final_scripts_odontocetes/08_results_gams/response_curves.R
```

The 08 scripts always report the sea-ice models and the habitat models of
`ODONTOCETE_RESULTS_STAGE` (default `stage3`, file names `main_stage3_*`; set
`stage1` for the supplementary tables, file names `supplementary_stage1_*`).

## Inputs

* `F_MODEL_INPUT` (`05_environmental_data/model_input_daily.csv`): recording days,
  daily environmental variables, presence minutes per species, `rec_minutes`.
* `F_ENV_DAILY` (`05_environmental_data/env_daily_wide.csv`): the full
  climatology window (deployment years ± `CLIM_HALF_WINDOW_YEARS`). Stage 1 stops
  if any deseasoning-eligible variable has fewer than `CLIM_MIN_YEARS` years or
  misses a calendar month in that window.

## Outputs

`outputs/06_statistical_analyses/stage1|stage2|stage3|sea_ice/` and
`outputs/08_results_gams/`; file lists are in each script header.

## Testing overrides

Environment variables `ODONTOCETE_OUTPUT_ROOT`, `ODONTOCETE_ENV_DAILY`,
`ODONTOCETE_MODEL_INPUT` redirect inputs/outputs; `ODONTOCETE_ALLOW_SHORT_CLIMATOLOGY=TRUE`
lets a short climatology run with a loud warning. These exist only to test the
code before the 10-yr downloads exist; production runs use the defaults.
