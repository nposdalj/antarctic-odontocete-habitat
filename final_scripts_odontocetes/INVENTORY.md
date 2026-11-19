# Script Inventory — Antarctic Odontocete Habitat

Generated 2026-11-16. Covers all 83 `.R` files in the repo (no `.Rmd` files exist). Originals are untouched.
"Fit" column maps each script to manuscript sections: 2.1 Study area · 2.2 Deep oceanography · 2.3 Sea ice · 2.5 Env. data · 2.5.1 HYCOM vs Copernicus · 2.6 Statistics · 3.1 Occurrence results · 3.2 GAM results · — = not for odontocete manuscript.

## Core pipeline (top level, Autocorrelation, Environmental Data)

| Path | Summary | Inputs | Outputs | Fit |
|---|---|---|---|---|
| FinalData.R | Builds the daily site × date modeling table: merges daily odontocete presence with AAO, Copernicus PHY/BGC (single depth levels 0.5–1665 m), STL/anomaly deseasoned surface series, AVISO FSLE, and 1–6 month lags; pivots depth vars to wide (`var_<depth>`); adds `ice_diff`, `ice_regime`, `bathymetry`, `julian_day`. Sets **all** NAs to 0 (L85). Then plots time series. | `data/Antarc_Odontocetes_daily_5min_bins.csv`; `Daily_AAO.csv`; `Copernicus/copernicus_40km.csv`; `Copernicus/copernicusLagged.csv`; `Copernicus/Deseasoning/{site}_deseasoned.csv`; `AVISO/{site}_fsle_40km` | `data/allData.csv` (absolute path, L193); on-screen plots | 2.5, 2.3, 3.1 |
| FinalData_Updated.R | **Plotting only — does not read or write data.** Requires `allData` in memory; filters to surface (depth 0.5), converts presence to minutes/day by guessing units from max value; saves 3-panel env + presence figures per site × species. Lines 1–137 are duplicated at 139–273. | in-memory `allData` | `plots/timeseries_3panel/{site}_{sp}_{set}_surface.png`; `saved_plot_files_surface.txt` | 3.1, 2.3 |
| Autocorrelation/ACF_SiteSpecies.R | Builds hourly/daily presence tables per site (Call ∈ Clicks/regular-click); fits binomial GLM `sp ~ bs(Day, k=4)` per site × species and records the first lag where residual ACF falls within the 95% CI; counts detection days. Note `bs()` has no `k` arg (partial-matches `knots`). | `Data/Antarc_Odontocetes.csv` | `Data/dailyDetections.csv`; `Autocorrelation/acf_table.csv`; `binned_detections.csv` | 2.6, 3.1 |
| Environmental Data/GetAVISO.R | Reads daily AVISO FSLE NetCDFs, crops to hardcoded ±40 km box per site (lat ±0.361°, lon ±0.74°), computes daily spatial mean/SD of `fsle_max` and `theta_max`. Bug L80: orientation SD computed from FSLE. | `D:/FSLE/AVISO_Global_FSLE/dt_global_allsat_madt_fsle_YYYYMMDD_20180704.nc` | `AVISO/{EI,KGI,CI}_fsle` (no `_40km` suffix, no extension); plots | 2.5 |
| Environmental Data/GetCopernicus.R | Reads pre-downloaded 40 km-box PHY (`cmems_mod_glo_phy_my_0.083deg_P1D-m`; thetao, so, uo, vo, zos, mlotst, siconc, sithick) and BGC (`cmems_mod_glo_bgc_my_0.25deg_P1D-m`; chl, o2, nppv) NetCDFs; daily spatial mean/SD (MAD for velocities) per depth; matches PHY/BGC depths; computes EKE and EKE_mad. | `{site}40km_cmems_mod_glo_phy_my_*.nc`, `{site}40km_cmems_mod_glo_bgc_my_*.nc` | `Copernicus/copernicus.csv` (10 depths); `Copernicus/copernicus_depths.csv` (all matched depths incl. 54, 130, 183, 222 m); plots | 2.5, 2.3, 2.2 |
| Environmental Data/TimeLaggedCopernicus.R | 1–6 "month" lags (30-day steps; single-day value at date − 30×lag) of surface chla, productivity, temp, salinity, EKE per deployment day. | `Copernicus/lagged/{site}-lag_*.nc` | `Copernicus/lagged/copernicusLagged.csv` (FinalData reads a different path) | 2.5, 2.6 |
| Environmental Data/DeseasoningMethods.R | Compares additive/multiplicative decomposition, STL, and day-of-year anomaly on ~10 yr (±5 yr) surface chla, npp, sst, sss, mlayer; saves deseasoned tables; tests raw vs STL vs anomaly in binomial GAMs (BW29, CI). | `Copernicus/Deseasoning/{site}-deseasoning_*.nc`; `Data/allData_40km.csv` | `Copernicus/Deseasoning/{site}_deseasoned.csv`; plots | 2.6 |
| Environmental Data/GetHYCOM.R | Reads pre-downloaded daily HYCOM `.mat` files (~1° box per site), daily spatial mean temp/salt per depth. Writes to drive root (bug). | HYCOM `W303E304S61N59_*.mat` etc. (Trisha's Mac path) | `/HYCOM.csv`; plots | 2.5.1 |
| Environmental Data/Cop_vs_HYCOM.R | Compares surface temp/salinity from HYCOM (depth 0) and Copernicus PHY against Aqua MODIS SST: overlays; Shapiro-Wilk, paired t, Wilcoxon, Diebold-Mariano (AE loss). Notes SMOS SSS "seemed unreasonable". | HYCOM `.mat`; `{site}_ALLsurf_cmems_mod_glo_phy_my_*.nc`; `SST_ERDDAP/{site}_SST.csv` | Plots + printed tests only | 2.5.1 |
| Environmental Data/GetERDDAP_SST.R | Downloads Aqua MODIS daily SST (`erdMH1sstd1dayR20190SQ_Lon0360`) for ±0.5° box; daily mean. | ERDDAP | `{site}_SST_{yr}.nc`; `SST_ERDDAP/{site}_SST.csv` | 2.5.1 |
| Environmental Data/GetERDDAP_SSS.R | Downloads SMOS daily SSS (`coastwatchSMOSv662SSS1day`) ±0.5° box; daily mean. Not used downstream. | ERDDAP | `Salinity_SMOS/*` | 2.5.1 (validation) |
| Environmental Data/GetERDDAP_Chla.R | Downloads Aqua MODIS daily chl-a (`erdMH1chla1day_Lon0360`) ±0.5° box. Not used downstream. | ERDDAP | `{site}_ChlA_{yr}.nc`; plots | 2.5 (validation) |
| Environmental Data/Chla-OCNET.R | Unfinished OCNET chl-a crop/average; undefined variable, placeholder output path. | `I:/Chlorophyll - OCNET/*.nc` | none functional | — |
| Environmental Data/HabtiatAssociations.R | Exploratory: presence vs absence histograms, lag heatmaps, mean-vs-SD phase plots. Several saves reference undefined objects. | `data/allData.csv` | `F:/Antarctica/hist_*.png` | 2.6 (exploratory) |
| Environmental Data/Copernicus/Extended_.../Extending_Env_TimeSeries.R | KGI only: circular 40/60/100 km haversine masks on extended PHY/BGC; daily means + EKE_mad; deseasons by cyclic GAM on day-of-year (`bs="cc"`, k=20) and STL. Built for krill side project. | extended PHY/BGC `.nc` | `KGI_{40,60,100}km_mask_with_deseasoned.csv`; PNGs | 2.6 (sensitivity) / — |
| Environmental Data/Copernicus/Extended_.../Visualizaing_MAD_EKE.R | KGI 60 km, 2016: maps EKE field on lowest/highest MAD days; illustrates EKE_mad metric. | PHY `.nc` (u, v) | `EKE_MAD_visualizations/*.png` | 2.5 (supp. figure) |

## GAMs (current lineage)

| Path | Summary | Inputs | Outputs | Fit |
|---|---|---|---|---|
| GAMs/Gm_surface_depth_GAM_tweedie.R | Original hand-built Gm Tweedie workflow: bins from `acf_table.csv`, manual VIF and backward elimination, raw predictors. Bug: CI models fit on `KGI_binned`. | `data/allData.csv`; `acf_table.csv` | printed plots only | 2.6 (history) |
| GAMs/Gm_surface_depth_GAM_tweedie_variability.R | Adds in-script block-model ACF (EI, CI), ice filter (ice_conc < 0.15), variability predictors. Same CI/KGI data bug. | same | printed plots | 2.6 (history) |
| GAMs/Gm_surface_depth_GAM_tweedie_variability_deSeasoned.R | Predictors converted to site × month-mean anomalies; manual VIF/elimination. Fixes CI bug. | same | printed plots | 2.6 (history) |
| GAMs/Gm_..._deSeasoned_moreAutomated.R | First automated version (per-site copy-paste): block-model ACF bins, 95th-pct ice filter, monthly-mean deseasoning, `build_preds` → VIF ≤ 5 → best-AIC per theme → `auto_gam` (backward, p ≤ 0.05, `select=TRUE`). Depths 0/16/635. Interactive save prompt. | `allData.csv` | `F:/Antarctica/GAMs/GAM_<site>_<date>.png` (on prompt) | 2.6, 3.2 |
| GAMs/Oo_..._deSeasoned_moreAutomated.R | Same as Gm with Oo, depths 0/455; KGI block model `s(SSH)` only. | same | same | 2.6, 3.2 |
| GAMs/Pm_..._deSeasoned_moreAutomated.R | Gm template for Pm, depths 0/375/1665, KGI dropped (no detections); non-interactive save. **Bug: anomaly lists use 1655 instead of 1665**, so deepest anomalies silently skipped. | same | `F:/Antarctica/GAMs/GAM_Pm_<site>.png` | 2.6, 3.2 |
| GAMs/Automated_GAM_tweedie_variability_deSeasoned_allSpecies.R | Config-driven single-species pipeline looping over sites: block-model ACF bins → ice filter (odontocetes) → monthly-mean deseasoning → `build_preds` → VIF ≤ 5 → theme selection → `auto_gam` → plots → diagnostics; saves everything. Configured for Mn. | `F:/Antarctica/GAMs/allData_humpback.csv` | `pre_model_checks/`, `final_gams/` (PNG, summary, term CSVs, metrics), `diagnostics/` | 2.6, 3.2 |
| GAMs/DEPTH_BINNED_Automated_..._allSpecies.R | Newest (2026-07-02). Same as above plus: baleen depth bins 0–54 m (`111054`) and 130–222 m (`130222`); site-specific ACF block-model terms (incl. julian_day); **periodic cubic spline monthly climatology deseasoning** (`splinefun(method="periodic")`); stacked time-series step. Configured for Mn on L: drive. | `L:/.../mysticetes/allData_humpback.csv` | same folders + `stacked_timeseries/` | 2.6, 3.2 |

## GAMs (superseded)

| Path | Summary | Inputs | Outputs | Fit |
|---|---|---|---|---|
| GAMs/seasonalModels/{BW29,BW37,Gm,Oo,Pm}_surface_depth_GAM_tweedie_seasonalModel.R (5) | Tweedie `tw(a=1.1,b=1.9)`, REML (no select) GAMs with **sea-ice predictors only** (ice_conc, ice_diff, ice_thickness); bins from in-script residual ACF; manual term removal. Finals mostly `s(ice_conc)`; CI often `s(ice_thickness)`. Mar 2026. | `data/allData.csv` | none (plots in session) | 2.6 / 3.2 — the **sea-ice model** of the two-model design (presence vs. sea ice on all bins); rewritten in 06 with `s(ice_conc, k = 4)` |
| GAMs/Trishas Old Models/BW29_depthsGAM.R | BW29 surface + 768 m, ACF bins, chla/NPP lag 0–6 screening, VIF, stepwise binomial GAMs. | `allData_40km.csv`; `acf_table.csv` | none | — (superseded) |
| GAMs/Trishas Old Models/BW29_surfaceGAM.R | BW29 surface-only binomial GAMs with presence upweighting. | same | none | — |
| GAMs/Trishas Old Models/BW37_depthsGAM.R | BW37 (67/920 m) binomial workflow. | same | none | — |
| GAMs/Trishas Old Models/BW37_surfaceGAM.R | BW37 surface binomial GAMs. | same | none | — |
| GAMs/Trishas Old Models/Gm_depthsGAM.R | Gm (16/635 m) binomial workflow. | same | none | — |
| GAMs/Trishas Old Models/Gm_surfaceGAM.R | Longest Gm binomial workflow (2147 lines). | same | none | — |
| GAMs/Trishas Old Models/Gm_surfaceGAM_tweedie.R | Bridge: first Tweedie version of Gm surface model. | `allData.csv` | none | — |
| GAMs/Trishas Old Models/OLDGm_GAM.R | Earliest Gm script (±0.5° env box); header documents problems and move to 40 km box. | `allData.csv` (hfdeg era) | none | — |
| GAMs/Trishas Old Models/Oo_depthsGAM.R | Oo at CI (11/455 m) binomial. | `allData_40km.csv` | none | — |
| GAMs/Trishas Old Models/Oo_surfaceGAM.R | Oo at CI surface binomial; final `s(ice_conc)`. | same | none | — |
| GAMs/Trishas Old Models/Pm_depthsGAM.R | Pm at EI, CI (375/1665 m) binomial. | same | none | — |
| GAMs/Trishas Old Models/Pm_surfaceGAM.R | Pm at EI, CI surface binomial. | same | none | — |
| GAMs/Trishas Old Models/allSite_surfaceGAM.R | Pooled-site BW29, Gm, BW37 binomial GAMs incl. julian_day, AAO. | `allData_40km.csv`; `acf_table.csv` | none | — |

## Data Visualizations

| Path | Summary | Inputs | Outputs | Fit |
|---|---|---|---|---|
| Data Visualizations/siteMap.R | Map of EI/KGI/CI on GMRT bathymetry, mako fill, isobaths every 500 m, magenta site points. Theme, scale bar, and 40 km boxes commented out; polar projection marked "IGNORE". No save. | `Environmental Data/GMRT/GMRTbathymetry_lowerres_wideView.grd` | none (existing `SiteMap*.png` likely manual exports) | 2.1 |
| Data Visualizations/TimeseriesPlots.R | Merges three per-site detection logs; daily hours by call type; weekly days detected per site/species with effort shading; builds 5-min presence grid → daily minutes. Assumes continuous recording (ignores duty cycle). Bug: `_daily_5min_bins.csv` written from the minutes table (files are identical). | `data/Antarc_{EIE,SSI,EI}_01_Odontocetes.csv`; `data/dailyDetections.csv` | `data/Antarc_Odontocetes.csv`; `data/Antarc_Odontocetes_daily_minutes.csv`; `data/Antarc_Odontocetes_daily_5min_bins.csv` | 3.1, 2.6 |
| Data Visualizations/DepthProfiles-Copernicus.R | Date × depth profiles of temp, salinity, O2, density (`oce::swRho`) per deployment window. Bug: salinity profile filters temperature. No save. | `Copernicus/copernicus_depths.csv` | ggplot objects | 2.2 |
| Data Visualizations/antarctic_3D_TS_pipeline_sitecolor_sigma_icefilter.R | 3D O2–T–S scatters, LDA/PCA by site, T–S diagrams with σθ contours (`oce::swSigmaTheta`), depth profiles, time–depth sections; keeps only ice_conc ≤ 0.15. No water-mass thresholds. | `data/allData.csv` | PNG/HTML in `F:/3D_T_S_0/10302025/` | 2.2 |
| Data Visualizations/antarctic_depth_timeseries_presence_minutes.R | Per site × species × species depth: stacked env time series + presence minutes bars. First half dead code. | in-memory long-format `allData` | `plots/timeseries_3panel/*_depth_*.png` | 3.1, 2.5 (supp.) |
| Data Visualizations/antarctic_depth_timeseries_presence_minutes_baleen.R | Same for mysticetes (Bp). | `L:/.../allData_fin.csv` | PNGs | — |
| Data Visualizations/FSLE_GIF.R | Daily FSLE maps + GIF (July 2016). Uses `lon_min` before definition. | AVISO FSLE `.nc` | PNGs | 2.5 (presentation) |
| Data Visualizations/Nicole_FSLE_Magnitude_Orientation_GIFs/FSLE_Mag_Orient_GIF.R | Cleaner two-panel FSLE magnitude/orientation maps with GMRT isobaths and sites; GIF for EI 2014. | AVISO `.nc`; GMRT `.grd` | daily PNGs; `FSLE_2014.gif` | 2.5 (supp./presentation) |
| Data Visualizations/SeaIceGIF.R | Stitches Uni-Bremen sea-ice screenshots into GIF. | `Presentation Plots/Ice Images/*.png` | GIF object | 2.3 (presentation) |

## Mysticetes (companion baleen paper)

| Path | Summary | Inputs | Outputs | Fit |
|---|---|---|---|---|
| Mysticetes/FinalDataBlueZ.R | Copy of old FinalData.R for Bm Z-calls (hours/day). | baleen data paths; `copernicus_depths.csv`; FSLE 40 km | `mysticetes/allData.csv` | — |
| Mysticetes/FinalDataFin20.R | Same for fin 20 Hz (Bp) at 130 m. | same | plots | — |
| Mysticetes/FinalDataHumpback.R | Same for Mn (minutes/day) at 130 m. | same | plots | — |
| Mysticetes/Mn_surface_depth_GAM_tweedie_seasonalModel.R | Copy of BW37 seasonal template; GAMs still fit `BW37 ~` (bug). | `I:/allData_humpback.csv` | none | — |

## Proposal, Krill, Sophia (side projects)

| Path | Summary | Inputs | Outputs | Fit |
|---|---|---|---|---|
| Proposal/SiteMap.R | Map function (GMRT highest-res, mako, isobaths, ggrepel labels, scale bar) for proposed sites; KML includes EI/KGI/CI points. | GMRT highest-res `.grd`; `Antarctica.kml` | `Site_Locations.csv`; region maps PNG/PDF | 2.1 (styling template) |
| Proposal/SiteMap_SixOnly.R | Single map of 8 proposed sites (no EI/KGI/CI). Missing `+` bug. | GMRT `.grd` | `Antarctic_Proposed_Sites_Map.png/.pdf` | — |
| Proposal/RegionalSeaIce.R | Sea-ice conc/thickness within 25 km of proposed sites, last 3 yr. | `Proposal/*.nc`; KML | PNG/PDF | — |
| Proposal/Krill_EKE_MAD.R | Krill catch near KGI vs EKE MAD, 30-day rolling. | KGI 60 km CSV; CCAMLR `.Rds` | plot | — |
| Krill/CCAMLR.R | CCAMLR Statistical Bulletin krill/fishery catch bar plots. | `AggregatedKrillCatch.csv` etc. | PNGs | — |
| Krill/KRILLBASE.R | KRILLBASE hauls within 600 km of sites: bubble maps and time series. | `krillbase_data.csv` | PNGs | — (Discussion at most) |
| Krill/SeaIce_and_Krill.R | GLORYS physics + sea ice box means joined to CCAMLR catch (APEI, APDPE). | PHY `.nc`; catch CSV | `{APEI,APDPE}_krill_env_*.csv` | — |
| Krill/SIC vs Krill Catch Series APDPE.R | Catch vs SIC scatter + lm, lags 0–3 yr. | APDPE CSVs | PNGs | — |
| Krill/SIC vs Krill Catch Series APEI.R | Identical to APDPE version for APEI. | APEI CSVs | PNGs | — |
| Sophia/Plot_OBIS.R | Downloads OBIS fin whale records for WAP box (no plot). | OBIS API | CSV/GPKG (placeholder path) | — |

## Sofia (krill-fishery side project)

| Path | Summary | Inputs | Outputs | Fit |
|---|---|---|---|---|
| Sofia/CCAMLR R code/Habitat Modeling/All Whale Detections with Zero Recording effort.R | Duplicate of Time Series version (comment differences only). | — | — | — |
| Sofia/CCAMLR R code/Habitat Modeling/Bathymetry working code.R | CCAMLR 756 krill-haul maps on NOAA bathymetry (`marmap`); approximate island points, not HARPs. | CCAMLR `.Rds`; NOAA bathy | plots | — (2.1 code idea only) |
| Sofia/CCAMLR R code/Habitat Modeling/Krill Catch for {Clarence,Elephant,King George} Island.R (3) | Bathymetry + krill catch within ~70 km of each island. | CCAMLR `.Rds`; NOAA bathy | plots | — |
| Sofia/CCAMLR R code/Sofie's Models/(USE THIS) GAMs with historical Mesoscale data.R | Gaussian GAMs of sqrt(daily krill catch) within 60 km of KGI vs deseasoned Copernicus anomalies; ACF, VIF, backward selection. | `KGI_60km_mask_with_deseasoned.csv`; CCAMLR `.Rds` | console/plots | — |
| Sofia/CCAMLR R code/Sofie's Models/GAMs with historical Mesoscale data3.R | Byte-identical to "(USE THIS)". | — | — | — |
| Sofia/CCAMLR R code/Sofie's Models/GAMS with historical mesoscale data (DRAFT ONE).R | Earlier non-deseasoned draft. | `KGI_60km_mask.csv`; `.Rds` | console/plots | — |
| Sofia/CCAMLR R code/Sofie's Models/GAMs with autocorrelation for krill hauls (DRAFT 2) .R | Single-variable krill GAMs + residual ACF. | `KGI_deseasoned.csv`; `.Rds` | console/plots | — |
| Sofia/CCAMLR R code/Sofie's Models/Export haul catch table (756, no buffer).R | Haul × taxon catch table for subareas 756/759. | CCAMLR `.Rds` | `CCAMLR_756_759_haul_catch_table.xlsx` | — |
| Sofia/CCAMLR R code/Time Series/All Whale Detections with Zero Recording effort.R | **Mysticetes only**: daily faceted plot of detections, krill catch within 100 km, env lines, effort (row count/day) with zero-effort shading. | `D:/Mystecedes Time Series/allData_*.csv`; `.Rds` | plot | — (3.1 plotting idea) |
| Sofia/CCAMLR R code/Time Series/All Whale Detections with Zero Recording effort_datingto2011.R | Despite name: lm of KGI deseasoned env var vs krill catch from 2011. | `KGI_deseasoned.csv`; `756_USA.Rds` | plot | — |
| Sofia/CCAMLR R code/Time Series/All whale detections not seperated by island.R | Same as zero-effort plot without site buffer on krill. | same | plot | — |
| Sofia/CCAMLR R code/Time Series/Krill Catch vs ENV_VAR panels.R | Krill catch vs 5 deseasoned env vars (loads CI file but buffers KGI — mismatch). | `CI_deseasoned.csv`; `.Rds` | plot | — |
| Sofia/DailyDetections/PlotData.R | Base-R daily presence spike plots (Pm, Oo at CI). | `data/dailyDetections.csv` | plots | 3.1 (minimal) |
| Sofia/Krill prev Semester/SIC vs Krill Catch Series APDPE.R | SIC vs krill catch loess. | `APDPE_krill_env_monthly.csv` | plot | — |
| Sofia/Krill prev Semester/SIC vs Krill Catch Series APEI.R | Same for APEI. | `APEI_krill_env_monthly.csv` | plot | — |
| Sofia/Krill prev Semester/SIT scaling.R | Weekly Pm presence at CI with sea-ice thickness overlay (column-position selection). | `dailyDetections.csv` | plot | 3.1 / 2.3 (marginal) |

## Key-parameter verification (briefing table vs. code)

| Parameter | Expected | Found in code | Status |
|---|---|---|---|
| Spatial box | 40 × 40 km | ±40 km half-width (~80 × 80 km), hardcoded corners in GetAVISO.R L18–39; Copernicus clipped at download (`{site}40km_*`) | **DIFFERENT** |
| Surface depth bin | 0–54 m | `111054` columns only in DEPTH_BINNED GAM and only for baleen species; odontocetes use point depths (0/375/1665 in allSpecies; Gm 0/16/635, Oo 0/455, Pm 0/375/1665 in species scripts). No script in repo produces the binned columns. | **NOT FOUND for odontocetes** |
| Sub-surface depth bin | 130–222 m | `130222`, same as above | **NOT FOUND for odontocetes** |
| Deseasoning | 10-yr monthly climatology, periodic spline | Periodic spline only in DEPTH_BINNED (L528), climatology computed from the modeling data per site (deployment months), not 10 yr. Other GAMs: site × month mean subtraction. DeseasoningMethods.R uses ~10 yr but STL/day-of-year methods. | **DIFFERENT** |
| Lag range | 1–6 mo per predictor | Lags 1–6 mo exist in data (single-day values at 30-day steps), but GAMs keep only chla/productivity at 3 and 6 mo (`keep_lags <- c(3, 6)`), competing within "productivity" theme by AIC | **DIFFERENT** |
| ACF bins | EI 3, KGI 4, CI 4 d | Computed at runtime per site × species (first lag where block-model residual ACF < 95% CI). Saved `acf_table.csv`: no column matches 3/4/4. | **NOT FOUND** |
| VIF | ≤ 5 sequential | `vif_stepwise_drop`, `car::vif` on gaussian GLM, GVIF^(1/(2Df)), threshold 5 | FOUND |
| GAM family | `tw(link="log")`, REML, select=TRUE | `tw(link="log", a=1.1, b=1.9)`, REML, `select=TRUE` in final `auto_gam`; `s(x, k=4)` TPRS; one model per site; no offset/weights | FOUND |
| Backward selection | p ≤ 0.05 | `auto_gam`, drop max p until all ≤ 0.05 (last term always kept) | FOUND |
| FSLE | AVISO+, 0.04°, backward | `dt_global_allsat_madt_fsle_*_20180704.nc`, `fsle_max`/`theta_max`; resolution & direction not stated in code | PARTIAL |
| PHY product | PHY_001_030, 1/12° | `cmems_mod_glo_phy_my_0.083deg_P1D-m` | FOUND |
| BGC product | BGC_001_029, 1/4° | `cmems_mod_glo_bgc_my_0.25deg_P1D-m` | FOUND |
| Sea ice source | GCOM-W1/NSIDC? | Copernicus PHY `siconc`/`sithick`; AMSR2/Bremen code exists but README says it does not work | Copernicus |
| Response variable | — | Daily presence **minutes** (5-min bins × 5), averaged within ACF bin; no duty-cycle correction | note |
| Pilot whale detection | — | No method in any script; all Gm log entries are Call = "Clicks" (262 events) | not in code |
