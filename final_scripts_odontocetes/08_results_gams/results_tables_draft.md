# Results tables — draft

*Generated from `08_results_gams/results_tables.R` outputs (final Stage 3 habitat models and sea-ice models, run 2026-11-19). Table numbers are placeholders. The CSV versions are in `outputs/08_results_gams/`.*

---

**Table [INSERT] (main text).** Final habitat models (Tweedie GAMs) for each species and site. Predictors are anomalies from the 10-year climatology except FSLE, FSLV orientation and EKE. Term deviance explained is the decrease in deviance explained when the term is removed and the model refit, so terms need not sum to the model total. Shape describes the fitted smooth. ‡ Low-sample model (10–20 bins with detections; at most two terms). † Term not significant when residual autocorrelation was modeled with an AR(1) process, or when refit by maximum likelihood (Supplement S[INSERT]). Bins were 4 days long (killer whales, 2 days).

| Species | Site | Bins (with detections) | Dev. expl. (%) | Term | edf | p | Term dev. expl. (%) | Shape |
|---|---|---|---|---|---|---|---|---|
| Sperm whale | EI | 34 (23) | 85.2 | FSLE | 2.30 | <0.001 | 13.1 | peaked |
|  |  |  |  | Mixed layer depth anomaly† | 1.45 | 0.004 | 8.4 | decreasing |
|  |  |  |  | Chl-a anomaly, 3-month lag† | 1.88 | <0.001 | 23.4 | U-shaped |
| Sperm whale | KGI | 89 (0) | — | not modeled (< 10 bins with detections) | | | | |
| Sperm whale | CI | 76 (65) | 53.1 | FSLV orientation | 1.61 | 0.007 | 5.5 | cyclic, peak ° |
|  |  |  |  | EKE, mid-water layer | 0.87 | 0.004 | 2.6 | decreasing |
|  |  |  |  | EKE, deep layer† | 1.58 | 0.021 | 2.5 | peaked |
|  |  |  |  | Mixed layer depth anomaly† | 1.52 | 0.029 | 4.8 | peaked |
|  |  |  |  | Temperature anomaly, surface layer | 1.79 | <0.001 | 37.2 | decreasing |
|  |  |  |  | Chl-a anomaly, 3-month lag | 1.76 | <0.001 | 10.6 | U-shaped |
| Long-finned pilot whale | EI | 34 (19)‡ | 84.8 | Temperature anomaly, surface layer | 2.78 | <0.001 | 21.1 | increasing |
|  |  |  |  | Chl-a anomaly, 3-month lag | 1.56 | <0.001 | 26.3 | increasing |
| Long-finned pilot whale | KGI | 64 (18)‡ | 31.9 | Temperature anomaly, surface layer | 1.64 | <0.001 | 31.9 | increasing |
| Long-finned pilot whale | CI | 64 (11)‡ | 48 | Temperature anomaly, surface layer† | 0.82 | 0.027 | 6.8 | decreasing |
|  |  |  |  | EKE, deep layer† | 1.50 | 0.010 | 15.7 | peaked |
| Killer whale | EI | 68 (3) | — | not modeled (< 10 bins with detections) | | | | |
| Killer whale | KGI | 177 (3) | — | not modeled (< 10 bins with detections) | | | | |
| Killer whale | CI | 151 (11)‡ | 21.1 | Chl-a anomaly, 3-month lag† | 1.57 | 0.007 | 21.1 | U-shaped |
| Southern bottlenose whale | EI | 34 (26) | 59 | FSLE† | 1.41 | 0.031 | 9.8 | peaked |
|  |  |  |  | Chl-a anomaly, 3-month lag† | 2.76 | <0.001 | 45.2 | complex |
| Southern bottlenose whale | KGI | 56 (20)‡ | 38.7 | Chl-a anomaly, 3-month lag | 1.93 | 0.003 | 38.7 | increasing |
| Southern bottlenose whale | CI | 52 (24) | 66.7 | FSLE | 1.44 | 0.003 | 11.8 | peaked |
|  |  |  |  | Chl-a anomaly, 3-month lag | 0.96 | <0.001 | 30.2 | decreasing |
|  |  |  |  | Mixed layer depth anomaly | 0.90 | 0.003 | 6.9 | decreasing |
|  |  |  |  | Temperature anomaly, surface layer | 1.87 | <0.001 | 27 | peaked |
| BW37 beaked whale | EI | 34 (17)‡ | 47.8 | FSLE | 1.84 | <0.001 | 34.7 | U-shaped |
|  |  |  |  | FSLV orientation | 1.35 | 0.029 | 11.5 | cyclic, peak ° |
| BW37 beaked whale | KGI | 56 (27) | 42.5 | Mixed layer depth anomaly† | 1.40 | 0.018 | 8.4 | U-shaped |
|  |  |  |  | Temperature anomaly, deep layer | 1.80 | <0.001 | 34.1 | peaked |
| BW37 beaked whale | CI | 52 (7) | — | not modeled (< 10 bins with detections) | | | | |

---

**Table S[INSERT] (supplement).** Sea-ice models (presence ~ s(sea-ice concentration)) for each species and site, using all bins. "Robust to AR(1)": whether the effect remained significant when residual autocorrelation was modeled with an AR(1) process (shown for significant effects only). ‡ Low-sample model (10–20 bins with detections).

| Species | Site | Bins (with detections) | edf | p | Dev. expl. (%) | Shape | Robust to AR(1) |
|---|---|---|---|---|---|---|---|
| Sperm whale | EI | 34 (23) | 0.73 | 0.064 | 6.5 | not significant |  |
| Sperm whale | KGI | 89 (0) | — | — | — | not modeled | |
| Sperm whale | CI | 76 (65) | 1.60 | 0.002 | 11.7 | peaked | no |
| Long-finned pilot whale | EI | 34 (19)‡ | 2.21 | <0.001 | 34.7 | peaked | no |
| Long-finned pilot whale | KGI | 89 (18)‡ | 1.67 | 0.183 | 27.9 | not significant |  |
| Long-finned pilot whale | CI | 76 (11)‡ | 0.60 | 0.208 | 3.8 | not significant |  |
| Killer whale | EI | 68 (3) | — | — | — | not modeled | |
| Killer whale | KGI | 177 (3) | — | — | — | not modeled | |
| Killer whale | CI | 152 (11)‡ | 1.66 | <0.001 | 28.5 | peaked | no |
| Southern bottlenose whale | EI | 34 (26) | 0.09 | 0.292 | 0.5 | not significant |  |
| Southern bottlenose whale | KGI | 89 (20)‡ | 0.85 | 0.041 | 32.9 | decreasing | yes |
| Southern bottlenose whale | CI | 76 (24) | 0.92 | <0.001 | 56.1 | decreasing | yes |
| BW37 beaked whale | EI | 34 (17)‡ | 0.83 | 0.201 | 5.3 | not significant |  |
| BW37 beaked whale | KGI | 89 (27) | 0.90 | 0.009 | 35.5 | decreasing | yes |
| BW37 beaked whale | CI | 76 (7) | — | — | — | not modeled | |

---

**Table S[INSERT] (supplement).** Diagnostics of the final habitat models. Tweedie p: estimated power parameter. Min k-index p: smallest p-value of the basis-dimension check across terms (low values suggest too few basis functions). Max concurvity: largest "worst-case" concurvity among terms. Lag-1 ACF and Ljung–Box p (lags 1–5 bins): residual autocorrelation; p ≤ 0.05 indicates remaining autocorrelation. ΔAIC with AR(1): change in AIC when an AR(1) process was added in glmmTMB (values < −2 indicate support for AR(1); — = AR(1) model did not converge).

| Species | Site | Bin (d) | Tweedie p | min k-index p | Max concurvity | Lag-1 ACF | Ljung–Box p | ΔAIC with AR(1) |
|---|---|---|---|---|---|---|---|---|
| Sperm whale | EI | 4 | 1.47 | 0.220 | 0.74 | 0.37 | 0.066 | -25.6 |
| Sperm whale | CI | 4 | 1.57 | 0.122 | 0.88 | 0.04 | 0.17 | -23.1 |
| Long-finned pilot whale | EI | 4 | 1.31 | 0.148 | 0.49 | -0.16 | 0.28 | — |
| Long-finned pilot whale | KGI | 4 | 1.32 | 0.948 | 0.00 | 0.44 | 4.4e-05 | -5.7 |
| Long-finned pilot whale | CI | 4 | 1.44 | 0.400 | 0.25 | 0.11 | 0.16 | 2.9 |
| Killer whale | CI | 2 | 1.35 | 0.660 | 0.00 | 0.12 | 0.34 | -19.1 |
| Southern bottlenose whale | EI | 4 | 1.50 | 0.295 | 0.71 | -0.21 | 0.35 | -3.6 |
| Southern bottlenose whale | KGI | 4 | 1.18 | 0.050 | 0.00 | 0.27 | 0.14 | — |
| Southern bottlenose whale | CI | 4 | 1.39 | 0.265 | 0.71 | 0.35 | 0.034 | 3.2 |
| BW37 beaked whale | EI | 4 | 1.10 | 0.540 | 0.26 | 0.01 | 0.47 | 2.0 |
| BW37 beaked whale | KGI | 4 | 1.17 | 0.530 | 0.06 | -0.10 | 0.58 | 2.2 |

---

### Notes for Natalie
- **Low-sample marker (‡) and the AR(1) marker (†)** are the two footnote symbols; change them to the journal's style later.
- **Pilot whales at CI:** both terms are marked † because they were not significant when refit by maximum likelihood even without AR(1) (p = 0.08 and 0.10); that is the "fitting method" category from the AR(1) check. Consider describing this model as weak.
- **Sea-ice diagnostics** (residual autocorrelation strong in 9 of 11 models) are in `supplementary_sea_ice_gam_diagnostics.csv`; add a table like the habitat one if you want them in the supplement.
- **Basis-dimension check:** no habitat model is flagged at p < 0.05 (southern bottlenose whale at KGI is at 0.050).
