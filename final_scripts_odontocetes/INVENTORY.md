# Script Inventory — Antarctic Odontocete Habitat

Generated 2026-11-16. Covers all 83 `.R` files in the repo (no `.Rmd` files exist). Originals are untouched.
"Fit" column maps each script to manuscript sections: 2.1 Study area · 2.2 Deep oceanography · 2.3 Sea ice · 2.5 Env. data · 2.5.1 HYCOM vs Copernicus · 2.6 Statistics · 3.1 Occurrence results · 3.2 GAM results · — = not for odontocete manuscript.

## Core pipeline (top level, Autocorrelation, Environmental Data)

| Path | Summary | Inputs | Outputs | Fit |
|---|---|---|---|---|
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
