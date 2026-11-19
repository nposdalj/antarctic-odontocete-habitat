# 2.6 Statistical Analyses — methods draft

*Drafted from `06_statistical_analyses/00_gam_parameters.R`, `00_gam_functions.R`, `01_stage1_automated_selection.R`, `02_stage2_decision_table.R`, `03_stage3_harmonized_models.R`, `04_sea_ice_models.R`, `05_ar1_sensitivity.R`, and `config.R`. Flags: `[INSERT]` = information not in the code or result-dependent; `[CONFIRM]` = value to verify. Rewritten 2026-11-19 (depth layers, anomalies, temperature-only water mass, revised selection order, Ljung–Box and AR(1) checks). Values are from the preliminary run; re-check after the final Stage 3 run.*

---

Relationships between odontocete acoustic presence and environmental conditions were modeled with generalized additive models (GAMs), consistent with the methods described for baleen whales in [COMPANION PAPER CITATION]. Separate models were fit for each species at each site. To reduce temporal autocorrelation, daily presence was aggregated into consecutive multi-day bins. For each species and site, a GAM of daily presence minutes on a smooth of days since the start of recording was fit, and the bin length was the first lag at which the autocorrelation of its residuals fell within the 95% confidence bounds. For comparability across sites, each species was assigned the longest bin length among its sites, capped at 4 days to retain at least 30 bins at the shortest deployment (EI, 135 days). Bins were 4 days long for all species except killer whales (2 days); the cap applied to sperm whales, southern bottlenose whales and BW37, whose residual autocorrelation lasted 6, 11 and 17 days, respectively. The response was the mean daily presence minutes within each bin. Species–site combinations with fewer than 10 bins containing detections were not modeled and are described only in Section 3.1 (sperm whales at KGI, killer whales at EI and KGI, and BW37 at CI). Models with 10–20 bins with detections (long-finned pilot whales at all sites, killer whales at CI, southern bottlenose whales at KGI, and BW37 at EI) were identified as low-sample models.

Two types of model were fit. Sea-ice models related presence to sea-ice concentration alone, using all bins, and describe the seasonal pattern of occurrence. Habitat models related presence to the remaining environmental predictors, using only bins with sea-ice concentration no greater than the highest concentration at which the species was detected (Section 2.3). Because the seasonal pattern is described by the sea-ice models, the habitat models were designed to test responses to conditions that were unusual for the time of year. All habitat predictors with a seasonal cycle were therefore expressed as anomalies from a 10-year climatology (deployment years ± 5 years): the 12 monthly means were interpolated to daily values with a periodic cubic spline and subtracted from the daily values. EKE, FSLE and FSLV orientation have almost no seasonal cycle (the seasonal cycle explained ≤ 8% of their daily variance at every site) and were used as raw values.

Each candidate habitat predictor represents a different process (Table [INSERT]). Water-column conditions were represented by the temperature anomaly of each depth layer the species uses: the surface layer (0–100 m) for all species, plus the mid-water layer (100–500 m) for sperm and killer whales and the deep layer (500–1000 m) for sperm whales, long-finned pilot whales, southern bottlenose whales and BW37, according to their foraging depths (Table [INSERT]; Section 2.2). Temperature was the only water-mass variable because, of the variables available, it was the one the reanalysis reproduced well in all three layers when compared with in-situ profiles (Section 2.5.1; r = 0.75–0.90 for layer means, versus 0.19 for mid-water salinity and ≤ 0.15 for oxygen). Eddy activity was represented by EKE in the mid-water and deep layers used by each species, and surface fronts by FSLE and FSLV orientation. Stratification was represented by the mixed layer depth anomaly. Productivity was represented by the anomaly of log10-transformed surface chlorophyll-*a* averaged over the 30-day window ending 3 months before each day. A single lag was used for all species: among the candidate lags, the 3-month lag was the most consistently associated with presence in exploratory single-term models (significant in 9 of 11 species–site models, versus 6 of 11 for a 5-month lag).

All models used a Tweedie distribution with a log link (power parameter estimated between 1.1 and 1.9), smoothing parameters estimated by restricted maximum likelihood (REML), and thin-plate regression splines with a basis dimension of 4. FSLV orientation was fit with a cyclic cubic spline joined at ±90°, because orientations differing by 180° are equivalent (Haller 2015). For each species and site, the candidates were screened for collinearity by sequentially removing the predictor with the highest variance inflation factor (VIF) until all were ≤ 5, and terms were then removed by backward selection, using additional shrinkage of smooth terms (`select = TRUE`), until all remaining terms were significant (*p* ≤ 0.05). To limit overfitting, the number of smooth terms was restricted following the guideline of at least 10 observations per estimated term (Peduzzi et al. 1996, Harrell 2015) [CONFIRM references]. Each model was allowed at most one term per 10 bins; in low-sample models (10–20 bins with detections), where the bins with detections rather than all bins carry most of the information, at most two terms were allowed (i.e., about one term per 10 bins with detections). The limit was applied during backward selection: once all remaining terms were significant, the least significant term continued to be removed until the model met the limit. A common set of predictors was then chosen for each species (the predictors selected at any of its sites), and the final model at each site was refit from this set with the same backward selection [CONFIRM wording once Stage 3 is final].

The VIF threshold of 5 follows Zuur et al. (2010). Model fit was assessed with residual diagnostics, basis-dimension checks and concurvity. Residual autocorrelation was tested with a Ljung–Box test over lags 1–5 bins, computed on the bin series with bins removed by the ice filter treated as gaps. As a sensitivity analysis, each final model was refit with an additional first-order autoregressive (AR(1)) process in glmmTMB (Brooks et al. 2017), keeping the same smooth terms (FSLV orientation as sin 2θ and cos 2θ). AR(1) was retained where it lowered AIC by more than 2, and each term was then tested by a likelihood-ratio test (Supplement S[INSERT]). Analyses were conducted in R 4.5.1 (R Core Team 2025) with the packages mgcv 1.9.3 (Wood 2017), car 3.1.3 (Fox & Weisberg 2019), gratia 0.11.1 (Simpson 2024) and glmmTMB 1.1.15.2 (Brooks et al. 2017) [CONFIRM versions after the final run].

---

**Table [INSERT] (draft).** Candidate predictors for the habitat models. Layers: surface 0–100 m (all species); mid-water 100–500 m (sperm whale, killer whale); deep 500–1000 m (sperm whale, long-finned pilot whale, southern bottlenose whale, BW37).

| Predictor | Represents | Layer / depth | Units | Raw / anomaly | Source |
|---|---|---|---|---|---|
| Potential temperature | Water mass (e.g., Antarctic Surface Water, Winter Water, Circumpolar Deep Water) | Mean of each layer used | °C | Anomaly | Copernicus PHY |
| FSLE (backward-in-time) | Surface fronts and stirring | Surface | day⁻¹ | Raw | AVISO+ |
| FSLV orientation (cyclic) | Orientation of surface fronts | Surface | ° | Raw | AVISO+ |
| EKE, ½(u′² + v′²), relative to the 10-year mean flow | Eddy activity at foraging depth | Mean of mid-water and/or deep layer | cm² s⁻² | Raw | Copernicus PHY |
| Mixed layer depth | Stratification | Surface | m | Anomaly | Copernicus PHY |
| log10 chlorophyll-*a*, 30-day mean ending 3 months earlier | Delayed productivity (prey) | Surface | log10 mg m⁻³ | Anomaly | Copernicus BGC |

Sea-ice concentration (Copernicus PHY) is used only in the sea-ice models and for the ice cutoff, not as a habitat-model candidate.

---

**Supplement S[INSERT] (draft): sensitivity to residual autocorrelation.** Refitting the final habitat models with an AR(1) process, where supported by AIC, retained [INSERT x] of [INSERT y] significant terms (Table S[INSERT]). [INSERT one sentence on which effects were not robust, e.g., lagged chlorophyll-*a* in models where AR(1) was supported.] Effects that are not robust are interpreted cautiously in the Discussion.

---

### Notes for Natalie

**Preliminary values to re-check after Stage 3:** bin lengths, not-modeled and low-sample combinations, Ljung–Box and AR(1) results.

**Decisions reflected here (all 2026-11-19):**
- Habitat predictors are anomalies, except EKE, FSLE and FSLV orientation.
- Depth layers replace species-specific depths.
- Temperature is the only water-mass variable (validation-based).
- SSH, salinity, oxygen, EKE MAD, surface EKE, productivity and same-day chl-a were removed as candidates.
- One common 3-month chl-a lag, log10-transformed.
- VIF screening runs after the theme winners are chosen.
- Ljung–Box diagnostic, plus the AR(1) check in the supplement.

**Stage 1/2 wording:** you planned not to describe Stage 1 and 2 in detail. The selection paragraph describes the automated screening because it is what produced the predictor sets; it can be shortened to "predictors were screened by single-term AIC within themes, collinearity (VIF ≤ 5) and backward selection", but keep a mention so the harmonized sets are reproducible.

**Chl-a lag sentence:** it says honestly that 3 months was chosen by comparing it with 5 months. If you'd rather give an ecological rationale first (prey response time), add a reference.

**Rare species–site combinations (for Results/Discussion, qualitative only):** killer whales at EI and KGI (3 bins with detections each) and BW37 at CI (7) are not modeled. Describe qualitatively when and under what conditions detections occurred (e.g., time of year, sea ice, temperature and EKE in the relevant layers on detection days). An exploratory model run (outputs/sensitivity_selection/C_model_rare) gave no significant terms for killer whales at EI, mid-water EKE and temperature for killer whales at KGI, and deep-layer EKE and surface-layer temperature for BW37 at CI; with 3–7 detection bins these fits are not reliable and should not be reported as models, but they can point to which conditions to describe.

**Decisions confirmed 2026-11-19:** keep the term limits and the 10-bin inclusion rule; keep lagged chl-a (its robustness is reported via the AR(1) check; a run without chl-a is in outputs/sensitivity_selection/D_no_chla); keep VIF ≤ 5 (as in the original analysis; cite Zuur et al. 2010).

**Residual autocorrelation limitation** (for the Discussion), once Stage 3 is final: "Residual autocorrelation (Ljung–Box, p ≤ 0.05) remained in [n] of [N] models; refitting with an AR(1) process retained [x] of [y] terms."

**References to complete:** companion paper; Zuur et al. 2010 (Methods Ecol. Evol. 1:3–14); Peduzzi et al. 1996 (J. Clin. Epidemiol. 49:1373–1379); Harrell 2015 (Regression Modeling Strategies, 2nd ed., Springer); Haller 2015; Wood 2017; Fox & Weisberg 2019 (car); Simpson 2024 (gratia); Brooks et al. 2017 (glmmTMB, *R Journal* 9:378–400); R Core Team 2025; Grimm 1987 (CONISS; cited in 2.2).
