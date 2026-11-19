# 2.5.1 Model comparison and validation — methods draft (with Results text and supplementary legends)

*Drafted from `05b_hycom_comparison/01_download_hycom.R`, `02_download_modis_sst.R`, `03_compare_hycom_copernicus.R`, `04_download_insitu_profiles.R`, `05_download_models_at_profiles.R`, `06_compare_models_insitu.R`, `config.R`, and their outputs (run 2026-11-19). Flags: `[INSERT]` = information not in the code; `[CONFIRM]` = value to verify. References: placeholders only.*

---

## Methods

To assess how much the water-column predictors depended on the choice of ocean model, and how well they reproduced observations, the Copernicus reanalysis was compared with an independent ocean model and with satellite and in-situ observations. HYCOM GOFS 3.1 (GLBv0.08; 1/12° × 1/25°, 40 depth levels) [REF] daily temperature and salinity were obtained for each site's recording period, using the 3-hourly time step closest to 12:00 UTC. No single HYCOM experiment spans 2014–2016, so the reanalysis (experiment 53.X) was used for 2014 and the analysis experiments 56.3 and 57.2 thereafter. Daily sea surface temperature from the Aqua MODIS sensor (4 km, daytime, quality-masked; NOAA CoastWatch ERDDAP dataset erdMH1sstd1dayR20190SQ) [REF] was used where at least 20% of the box's ocean pixels were valid. All products were averaged over the same ±20 km box around each site and matched by date. HYCOM was interpolated to the Copernicus model levels, and agreement was summarized as the mean difference (bias), root-mean-square error (RMSE) and Pearson correlation (r), for raw values and for anomalies from each product's monthly means over the matched days (Fig. S3).

Both models were also compared with in-situ profiles from the World Ocean Database [REF] in a wider region around the sites (60–63°S, 52–60°W) from 2014 to 2016: 915 profiles from instrumented seals, Argo floats, ship CTDs, bottle casts and expendable bathythermographs (Fig. S4a). Only observations and profiles that passed the World Ocean Database quality control were used. For each profile, both models were sampled at the nearest grid cell with water on the profile's day (Copernicus daily mean; HYCOM time step nearest the profile time). Observed profiles were interpolated to the Copernicus model levels within their measured depth range, and oxygen was converted from µmol kg⁻¹ to mmol m⁻³ using in-situ density. Bias, RMSE and r were computed at each level and for thickness-weighted means of the three depth layers (Section 2.2), using only profiles covering at least 80% of a layer and only profile–level combinations available in both models. Because both models assimilate some in-situ profiles (e.g., Argo), this comparison is not fully independent.

---

## Results (draft text for Section 3.x)

The Copernicus reanalysis agreed with the independent observations at least as well as HYCOM did, and it was used for all predictors. Against MODIS surface temperature, Copernicus had smaller warm biases than HYCOM at all sites (−0.01 to 0.18 °C vs 0.24 to 0.44 °C; Fig. S3a). At depth, the two models differed by a near-constant offset, with HYCOM up to ~1 °C warmer than Copernicus at EI (Fig. S3b). Against in-situ profiles, Copernicus reproduced layer-mean temperature well in all three layers (r = 0.89, 0.75 and 0.87 for the surface, mid-water and deep layers), with smaller biases than HYCOM in the mid-water (0.20 vs 0.38 °C) and deep layers (0.02 vs 0.20 °C) and similar skill at the surface (Fig. S4b). In contrast, Copernicus mid-water salinity was poorly reproduced (r = 0.19; HYCOM r = 0.52), and Copernicus oxygen showed no positive correlation with observations below the surface layer. Copernicus also provides sea-ice concentration, mixed layer depth and a matching biogeochemistry hindcast from one consistent multi-year system, whereas no single HYCOM experiment covered the study period. Temperature was therefore used as the only water-mass variable in the habitat models (Section 2.6).

---

**Figure S3 (draft legend; `FigS3_hycom_copernicus.png`).** Comparison of the Copernicus reanalysis with HYCOM and MODIS over each site's recording period (±20 km box). (a) Daily surface temperature from Copernicus (blue), HYCOM (orange) and Aqua MODIS (points; days with ≥ 20% valid pixels) at EI, KGI and CI. (b) Bias (HYCOM − Copernicus) and root-mean-square error of daily temperature and salinity at every Copernicus model level to 1000 m, by site; dashed lines mark the layer boundaries (100 and 500 m).

**Figure S4 (draft legend; `FigS4_insitu_validation.png`).** Validation of the Copernicus and HYCOM models against in-situ profiles from the World Ocean Database, 2014–2016. (a) Locations of the 915 profiles by instrument type, with the ±20 km boxes around the recording sites; bathymetry from GMRT. (b) Bias (model − observed, top) and root-mean-square error (bottom) of temperature, salinity and dissolved oxygen (Copernicus only) at each model level, computed on the same profiles and levels for both models (levels with ≥ 10 matched profiles). Dashed lines mark the layer boundaries (100 and 500 m).

---

### Notes for Natalie
- **Layer-mean metrics** (`outputs/05b_hycom_comparison/insitu_metrics_by_layer.csv`) could become a short supplementary table; the numbers in the Results paragraph come from it.
- **MODIS match-up days are few** in the ±20 km box (EI 10, KGI 35, CI 34), because of cloud and ice; the surface comparison is therefore indicative.
- **Profile counts by type:** seal tags 778, Argo 73, ship CTD 32, bottle 23, XBT 9; most seal-tag profiles are shallower than ~120 m, so the mid-water and deep layers rely mainly on Argo and ship CTDs (97 and 82 profiles).
- **[CONFIRM]** whether to mention that the Copernicus biogeochemistry hindcast does not assimilate oxygen (I believe it does not; check the product manual).
- **References to add:** HYCOM GOFS 3.1; MODIS Aqua SST (NASA OBPG) / ERDDAP; World Ocean Database 2023 (Boyer et al. / Mishonov et al.) [CONFIRM edition].
