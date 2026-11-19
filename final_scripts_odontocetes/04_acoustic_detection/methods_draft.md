# 2.4 Acoustic Detection — methods draft

*Drafted from `04_acoustic_detection/build_presence.R`, `duty_cycle_simulation.R`, `README.md`, `config.R`, and the simulation output. Flags: `[INSERT]` = information not in the code; `[CONFIRM]` = value to verify.*

---

For beaked whales, sperm whales (*Physeter macrocephalus*), and long-finned pilot whales (*Globicephala melas*), echolocation signals were first detected with an expert system using a two-step approach (Roch et al. 2011) and then verified by an expert analyst. For beaked whales at EI, a previously published dataset of echolocation signals was used (Baumann-Pickering et al. 2015); the same data processing methods were replicated for beaked whales at KGI and CI, and for sperm whale and long-finned pilot whale signals at all sites. Killer whale (*Orcinus orca*) signals were identified by an expert analyst through manual review of long-term spectral averages (LTSAs) [CONFIRM method and LTSA settings]. Beaked whale echolocation pulses were classified into southern bottlenose whale (*Hyperoodon planifrons*; BW29) and the BW37 and BW58 signal types. Previous work has hypothesized that the BW37 and BW58 echolocation pulses could be produced by Gray's (*Mesoplodon grayi*) and strap-toothed (*M. layardii*) beaked whales, respectively (Trickey et al. 2015, Baumann-Pickering et al. 2015, Barlow et al. 2021); Shepherd's beaked whale (*Tasmacetus shepherdi*) has also been proposed as a candidate for these unidentified signals (Barlow et al. 2021). Further research is still needed to confirm the species producing these echolocation pulses.

Only echolocation clicks were used to define presence; sperm whale slow clicks (n = 44 events) and killer whale high-frequency modulated signals (n = 3 events) were excluded. Verified detections were summarized into 5-min bins, and a bin was considered positive for a species if any detection overlapped it. Daily acoustic presence was quantified as the number of minutes in positive bins per day (positive bins × 5 min). Daily recording effort was calculated as the time within each deployment's recording period, accounting for partial first and last days and for the duty cycle at EI and KGI.

To assess whether duty cycling reduced detectability, the continuously recorded CI data were subsampled to the EI (5 min on every 6 min) and KGI (10 min on every 13.33 min) duty cycles at 100 random phase offsets, and the proportion of detection events, positive 5-min bins, and presence days retained was calculated for each species (Table S[INSERT]; Fig. S2). Subsampling retained ≥ 89.6% of positive 5-min bins and ≥ 89.6% of presence days on average for all species under either duty cycle (lowest: sperm whales under the KGI duty cycle; ≥ 96.5% of presence days for all other species), despite recording only 83.3% (EI) and 75.0% (KGI) of the time. Because losses were far smaller than proportional to recording time, presence was not scaled by recording effort.

---

**Table S[INSERT] (draft).** Duty-cycle subsampling of continuous CI recordings (mean % retained, range across 100 phase offsets).

| Species | Duty cycle | Events | 5-min bins | Presence days |
|---|---|---|---|---|
| Sperm whale | EI (5/6 min) | 87.1 (84.6–90.1) | 94.0 (92.7–95.4) | 93.5 (90.1–96.7) |
| Sperm whale | KGI (10/13.33 min) | 79.1 (75.9–84.4) | 89.6 (88.3–91.5) | 89.6 (84.1–94.0) |
| Long-finned pilot whale | EI | 100.0 | 99.1 (97.8–100) | 100.0 |
| Long-finned pilot whale | KGI | 95.0 (86.8–100) | 95.7 (90.5–99.3) | 97.3 (86.7–100) |
| Killer whale | EI | 98.2 (87.0–100) | 99.3 (97.8–100) | 100.0 |
| Killer whale | KGI | 94.3 (82.6–100) | 96.7 (94.0–98.5) | 100.0 |
| Southern bottlenose whale (BW29) | EI | 100.0 | 98.4 (97.5–100) | 100.0 |
| Southern bottlenose whale (BW29) | KGI | 90.0 (87.2–92.4) | 90.6 (88.3–92.6) | 98.5 (94.1–100) |
| BW37 beaked whale | EI | 100.0 | 99.3 (96.3–100) | 100.0 |
| BW37 beaked whale | KGI | 92.3 (83.3–100) | 92.5 (81.5–100) | 96.5 (90.9–100) |

BW58 was not detected at CI and could not be evaluated.

**Figure S2 (draft legend).** Effect of duty cycling on odontocete detections, estimated by subsampling the continuous recordings at Clarence Island (CI) to the duty cycles used at Elephant Island (EI; 5 min on every 6 min) and King George Island (KGI; 10 min on every 13.33 min) at 100 random phase offsets. Panels show the percentage of detection events, positive 5-min bins, and presence days retained relative to continuous recording. Points show the mean and bars the range across phase offsets. Dashed horizontal lines mark the proportion of time recorded under each duty cycle (83.3% for EI, 75.0% for KGI). Open symbols mark species not detected at the site whose duty cycle was simulated (sperm whale at KGI). BW58 is not shown because it was not detected at CI. Pm, sperm whale; Gm, long-finned pilot whale; Oo, killer whale; BW29, southern bottlenose whale; BW37, BW37 beaked whale.

---

### Notes for Natalie
- **Killer whale method** (manual LTSA review) is from the planning briefing, not your text — please confirm and add LTSA settings.
- **Species names:** I added scientific names (*P. macrocephalus*, *G. melas*, *O. orca*, *H. planifrons*, *M. grayi*, *M. layardii*, *T. shepherdi*) — please check against your usage.
- **References to complete:** Roch et al. 2011; Baumann-Pickering et al. 2015; Trickey et al. 2015; Barlow et al. 2021.
- **Duty-cycle table and figure** re-run 2026-11-18 (BW37 name updated; in-figure caption moved to the draft legend above).
- The ACF-based temporal binning used for modeling is described in 2.6, not here.
