# 03_sea_ice

> **Update 2026-11-23:** the models now take sea-ice concentration from the
> AMSR2 ASI product (University of Bremen, 3.125 km), set by
> `SEA_ICE_SOURCE <- "bremen"` in `config.R`. The box means were extracted
> outside this pipeline and are imported by
> `05_environmental_data/04c_import_bremen_ice.R` from
> `Environmental Data/Sea Ice Bremen/extracted/<site>/*.csv` into
> `outputs/05_environmental_data/bremen_ice_daily.csv`. `05_compile_dataset.R`
> then swaps them in for the Copernicus `ice_conc`. The import also writes a
> Bremen vs Copernicus comparison (`bremen_vs_copernicus_ice.csv/.png`).
> `ice_thickness` is still Copernicus and is not used in any model. Set
> `SEA_ICE_SOURCE <- "copernicus"` to restore the description below. The
> manuscript text (Section 2.3, Table tab:pred) still says Copernicus.

No scripts live in this folder. Sea ice is extracted together with the other
Copernicus model fields in
`05_environmental_data/04_extract_copernicus.R`.

- **Source:** Copernicus Marine GLOBAL_MULTIYEAR_PHY_001_030 (GLORYS12 physics
  reanalysis, 1/12°, daily), variables `siconc` (sea ice area fraction) and
  `sithick` (sea ice thickness, m), over the same ±20 km box around each site
  as all other predictors (±40 km comparison in 05_environmental_data).
- **Not a satellite product:** sea ice is NOT from GCOM-W1/AMSR2 or NSIDC (an
  earlier briefing assumed so; the AMSR2 code in the repository is unused).
- **Missing values:** in the source files `siconc`/`sithick` are NaN (their
  `_FillValue`) over open-water ocean cells and never defined over land. NaN
  over ocean cells (cells with sea surface height) is set to 0 before the box
  mean, so `ice_conc` is the ice-covered fraction of the ocean part of the box.
- **Outputs:** `ice_conc`, `ice_conc_sd`, `ice_thickness` in
  `outputs/05_environmental_data/copernicus_surface.csv`, carried into
  `env_daily_wide.csv` and `model_input_daily.csv` by
  `05_environmental_data/05_compile_dataset.R`. `ice_diff` and `ice_regime`
  from the original analysis are not produced.
- **Use:** (1) the only predictor in the sea-ice models
  (`06_statistical_analyses/04_sea_ice_models.R`, all bins); (2) the
  per-species ice filter for the habitat models: bins with mean ice
  concentration above the highest value at which the species was ever
  detected (pooled across sites) are removed (`ice_cutoff()` in
  `06_statistical_analyses/00_gam_functions.R`; values in
  `outputs/06_statistical_analyses/stage1/tables/ice_cutoffs.csv`). Sea ice
  is not a habitat-model candidate. Raw values are used (not deseasoned).
