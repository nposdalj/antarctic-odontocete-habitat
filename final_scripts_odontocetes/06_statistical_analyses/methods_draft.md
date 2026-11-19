# 2.6 Statistical Analyses — methods draft

*Drafted from `06_statistical_analyses/00_gam_parameters.R`, `00_gam_functions.R`, `01_stage1_automated_selection.R`, `02_stage2_decision_table.R`, `03_stage3_harmonized_models.R`, `04_sea_ice_models.R`, and `config.R`. Flags: `[INSERT]` = information not in the code or result-dependent; `[CONFIRM]` = value to verify.*

---

Relationships between odontocete acoustic presence and environmental conditions were modeled with generalized additive models (GAMs), consistent with the methods described for baleen whale targets in [COMPANION PAPER CITATION]. Separate models were fit for each species at each site. To reduce temporal autocorrelation, daily presence was aggregated into consecutive multi-day bins: for each species and site, a GAM of daily presence minutes on a smooth of days since the start of recording was fit, and the bin length was the first lag at which the autocorrelation of its residuals fell within the 95% confidence bounds. For comparability across sites, each species was assigned the longest bin length among its sites, capped at 4 days to retain at least 30 bins at the shortest deployment (EI, 135 days). Bins were 4 days long for all species except killer whales (2 days); the cap applied to sperm whales, southern bottlenose whales, and BW37, whose residual autocorrelation lasted 6, 11, and 17 days, respectively [values from the preliminary run]. The response was the mean daily presence minutes within each bin. Species–site combinations with fewer than 10 bins containing detections were not modeled and are described only in Section 3.1 (sperm whales at KGI, killer whales at EI and KGI, and BW37 at CI); models with 10–20 bins with detections (long-finned pilot whales at all sites, killer whales at CI, southern bottlenose whales at KGI, and BW37 at EI) were limited to two terms and identified as low-sample models [values from the preliminary run].

Two types of model were fit. Sea-ice models related presence to sea-ice concentration alone, using all bins. Habitat models related presence to the remaining environmental predictors, using only bins with sea-ice concentration no greater than the highest concentration at which the species was detected (Section 2.3). Candidate habitat predictors, grouped into seven themes, are listed in Table [INSERT]. Water-column variables were evaluated at the surface and at each species' typical foraging depths. Strongly seasonal predictors were deseasoned: for sea surface height, mixed layer depth, and temperature, salinity, and oxygen at each depth, a 10-year monthly climatology (deployment years ± 5 years) was interpolated to daily values with a periodic cubic spline, and anomalies from this climatology were used whenever the climatology explained more than 50% of the predictor's daily variance at a site; otherwise raw values were used. Sea-ice concentration, chlorophyll-*a*, productivity and their lags, EKE, FSLE, and FSLV orientation were always used as raw values.

All models used a Tweedie distribution with a log link (power parameter estimated between 1.1 and 1.9), smoothing parameters estimated by restricted maximum likelihood (REML), and thin-plate regression splines with a basis dimension of 4; FSLV orientation was fit with a cyclic cubic spline joined at ±90°, because orientations differing by 180° are equivalent (Haller 2015). Predictors were screened for collinearity by sequentially removing the predictor with the highest variance inflation factor until all were ≤ 5. Remaining predictors were grouped into themes (temperature, salinity, oxygen, mesoscale activity, front orientation, stratification, and productivity), and the predictor with the lowest Akaike information criterion in a single-term model was retained from each theme. Starting from these predictors, terms were removed by backward selection, using additional shrinkage of smooth terms (`select = TRUE`), until all remaining terms were significant (*p* ≤ 0.05) and the number of terms did not exceed one per 10 bins. A common set of predictors was then chosen for all models based on ecological relevance and consistency across species and sites [CONFIRM wording once harmonization is done], and the final models were refit with these predictors followed by the same backward selection. Model fit was assessed with residual diagnostics, basis-dimension checks, concurvity, and the autocorrelation of model residuals. Analyses were conducted in R 4.5.1 (R Core Team 2025) with the packages mgcv 1.9.3 (Wood 2017), car 3.1.3 (Fox & Weisberg 2019), and gratia 0.11.1 (Simpson 2024) [CONFIRM versions after the final run].

---

**Table [INSERT] (draft).** Candidate predictors for the habitat models. Depth-resolved variables were evaluated at the surface (0.5 m) and at each species' foraging depths: sperm whale 375 and 1665 m; long-finned pilot whale 16 and 635 m; killer whale 11 and 455 m; southern bottlenose whale 768 m; BW37 67 and 920 m (nearest model levels; Section 2.5). Raw/anomaly: "rule" = anomaly used when the 10-year seasonal cycle explains > 50% of daily variance at that site, otherwise raw.

| Theme | Predictor | Depths | Units | Raw / anomaly | Source |
|---|---|---|---|---|---|
| Temperature | Potential temperature | Surface + species depths | °C | Rule | Copernicus PHY |
| Salinity | Salinity | Surface + species depths | — | Rule | Copernicus PHY |
| Oxygen | Dissolved oxygen | Surface + species depths | mmol m⁻³ | Rule | Copernicus BGC |
| Mesoscale activity | Sea surface height | Surface | m | Rule | Copernicus PHY |
| | Eddy kinetic energy, ½(u′² + v′²), relative to the 10-year mean flow | Surface + species depths | cm² s⁻² | Raw | Copernicus PHY |
| | Within-box velocity variability, ½(MAD_u² + MAD_v²) | Surface + species depths | cm² s⁻² | Raw | Copernicus PHY |
| | FSLE (backward-in-time) | Surface | day⁻¹ | Raw | AVISO+ |
| Front orientation | FSLV orientation (cyclic) | Surface | ° | Raw | AVISO+ |
| Stratification | Mixed layer depth | Surface | m | Rule | Copernicus PHY |
| Productivity | Chlorophyll-*a* | ≤ 50 m | mg m⁻³ | Raw | Copernicus BGC |
| | Net primary productivity | ≤ 50 m | mg C m⁻³ d⁻¹ | Raw | Copernicus BGC |
| | Chlorophyll-*a*, 1–6-month lags | Surface | mg m⁻³ | Raw | Copernicus BGC |
| | Net primary productivity, 1–6-month lags | Surface | mg C m⁻³ d⁻¹ | Raw | Copernicus BGC |

Where the biogeochemistry grid had no water within ±20 km of the recorder at a species' depth, oxygen was taken from the deepest available level (KGI: 773 m in place of 920 and 1665 m; CI: 1046 m in place of 1665 m). Sea-ice concentration (Copernicus PHY) is used only in the sea-ice models and for the ice cutoff, not as a habitat-model candidate. *This table lists candidates; after Stage 3 harmonization, add a column or a second table with the final predictor set.*

---

### Notes for Natalie

**Things that will be filled in from the real run**
- Bin lengths, not-modeled and low-sample combinations, residual-autocorrelation counts and package versions are now filled from the preliminary run (19 Nov 2026). Re-check them after the final run.

**Wording to check**
- **Harmonization sentence** (paragraph 3): written to be accurate about the process — the automated screening informed the common predictor set. Adjust once you've done Stage 3, but please keep a mention of the screening (see our discussion).
- **Residual autocorrelation:** you chose to report it as a limitation. Suggested sentence for the Discussion or end of this section: *"Residual autocorrelation remained in 7 of 11 habitat models and 10 of 11 sea-ice models, so p-values for those models may be optimistic."* (preliminary run)
- **Length:** this is longer than the 1–3 paragraph target in the briefing; the statistics carry most of the decisions. Cut or move detail (e.g., the predictor list) to a table if needed.

**Differences from the original briefing (for your records)**
- Box: ±20 km (briefing: 40 × 40 km — consistent); original code used ±40 km.
- Depth bins 0–54 and 130–222 m are baleen-only; odontocetes use species-specific depths.
- ACF bins are computed per species (max across sites, cap 4 d), not the fixed EI 3 / KGI 4 / CI 4 d from the briefing.
- Lags: 1–6 months for chlorophyll-*a* and productivity only, as 30-day window means (described in 2.5).
- Sea ice: separate sea-ice models + habitat models without ice (two-model design).
- No effort offset (justified by the duty-cycle simulation in 2.4).

**References to complete:** companion paper; Haller 2015 (Annu. Rev. Fluid Mech. 47:137–162); Wood 2017 (GAMs book, 2nd ed.) [CONFIRM which Wood reference you prefer]; any VIF reference (e.g., Zuur et al. 2010).
