# 2.3 Sea Ice — methods draft

*Drafted from `05_environmental_data/01_download_copernicus.R`, `04_extract_copernicus.R`, `03_sea_ice/README.md`, `06_statistical_analyses/00_gam_functions.R` (`ice_cutoff()`), `04_sea_ice_models.R`, `config.R`, and the Stage 1 output `outputs/06_statistical_analyses/stage1/tables/ice_cutoffs.csv`. Flags: `[INSERT]` = information not in the code; `[CONFIRM]` = value to verify.*

---

Daily sea-ice concentration was obtained from the Copernicus Marine Service global ocean physics reanalysis (GLOBAL_MULTIYEAR_PHY_001_030, dataset version 202311; 1/12° horizontal resolution; accessed 17–18 November 2026) and averaged over a box extending ±20 km from each HARP (approximately 40 × 40 km). The product reports open-water ocean cells as missing rather than zero. These cells were therefore set to zero before averaging, so the daily box mean is the fraction of the whole box covered by ice, not the mean concentration of the ice-covered cells only.

Acoustic occurrence was first described across sea-ice conditions as the percentage of time bins with detections in four sea-ice concentration classes (0–5, 5–15, 15–50 and > 50%; Fig. [INSERT]). Sea ice was then used in two ways in the models. First, for each species and site, acoustic presence was modeled against sea-ice concentration using all recording bins (Section 2.6), to describe how occurrence varied with the seasonal advance and retreat of the ice. Second, the habitat models (Section 2.6) were meant to test the influence of other environmental conditions only when each species could plausibly be present. They were therefore restricted, for each species, to time bins whose mean sea-ice concentration was no greater than the highest bin-mean concentration at which that species was detected, pooled across sites. The cutoffs were 88.9% for killer whales, 84.5% for sperm whales, 33.4% for long-finned pilot whales, and 12.6% for both southern bottlenose whales and BW37 beaked whales (Table S[INSERT]). By construction, no bins with detections were excluded. The restriction removed no bins at EI. It removed 1 of 152 two-day bins for killer whales at CI, 12 of 76 (CI) and 25 of 89 (KGI) four-day bins for long-finned pilot whales, and 24 of 76 (CI) and 33 of 89 (KGI) four-day bins for each beaked whale. Sea-ice concentration was therefore not a candidate predictor in the habitat models.

---

**Table S[INSERT] (draft).** Sea-ice cutoffs for the habitat models and the number of time bins removed at each site. The cutoff is the highest bin-mean sea-ice concentration at which the species was detected at any site. Bin length was 4 days for all species except killer whales (2 days).

| Species | Cutoff (% ice) | EI removed / total | KGI removed / total | CI removed / total |
|---|---|---|---|---|
| Sperm whale (Pm) | 84.5 | 0 / 34 | — (not detected) | 0 / 76 |
| Killer whale (Oo) | 88.9 | 0 / 68 | 0 / 177 | 1 / 152 |
| Long-finned pilot whale (Gm) | 33.4 | 0 / 34 | 25 / 89 | 12 / 76 |
| Southern bottlenose whale (BW29) | 12.6 | 0 / 34 | 33 / 89 | 24 / 76 |
| BW37 beaked whale | 12.6 | 0 / 34 | 33 / 89 | 24 / 76 |

---

### Notes for Natalie
- **Values are from the preliminary Stage 1 run** (±20 km box, 18 Nov 2026). They will be regenerated if the models are re-run, so check them against `ice_cutoffs.csv` before submission.
- **The two beaked whales share a cutoff (12.6%)** because both species' highest-ice detection is in the same EI bin (starting 2014-07-11). That's correct, not a bug, but a reviewer may ask about it.
- **Interpretation:** the cutoff removes almost nothing for sperm and killer whales, which were detected under up to ~85–89% ice at CI. For pilot whales and the beaked whales it removes the winter ice-covered bins at KGI and CI (roughly 15–37% of bins).
- **Detection counts per model** (≥ 10 bins rule) belong in Section 2.6/Table [models]; not repeated here.
- **Source:** the original briefing listed GCOM-W1/NSIDC sea ice, but every script uses Copernicus `siconc` (PHY_001_030). An AMSR2/Bremen extraction exists in the repo, but its README says the code does not work, and it was never used. Switching to a satellite product would need new code and a re-run of the models.
- **Ice thickness** (`sithick`) is extracted but not used in any model (decision: concentration only).
- **Ice averaging fix:** the original pipeline averaged only ice-covered cells, which overstated concentration on partly ice-covered days (e.g., EI 2014-05-27: 0.053 vs. 0.003). The new extraction counts open water as 0.
- **Raw, not deseasoned:** ice concentration is used as raw values in both the ice filter and the sea-ice models (decision 2.6).
- **Box size:** ±20 km box, justified by detection ranges and the 1/4° biogeochemistry grid (Section 2.5).
- **Figure:** this section has no figure of its own; the ice time series appears as the top strip of the water-column figure (Section 2.2).
