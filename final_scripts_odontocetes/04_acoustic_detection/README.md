# 04_acoustic_detection

Odontocete signals were detected and verified **outside R**. This folder contains
only the processing of the verified detection logs into presence tables, plus a
duty-cycle simulation.

## Detection methods (performed outside R)

For beaked whales, sperm whales, and long-finned pilot whales, echolocation
signals were first detected with an expert system using a two-step approach
(Roch et al., 2011) and then verified by an expert analyst. For beaked whales
at site EI, a previously published dataset of echolocation signals was used
(Baumann-Pickering et al., 2015). These same data processing methods described
in Baumann-Pickering et al. (2015) were replicated for beaked whales at sites
KGI and CI, as well as sperm whales and long-finned pilot whale signals at all
sites.

Killer whale signals were identified by manual review of long-term spectral
averages (LTSAs) by an expert analyst. [CONFIRM killer whale method, LTSA
settings, and analyst]

## Detection logs (inputs, in `data/`)

| Site | Log file | Recording period (UTC) | Duty cycle |
|---|---|---|---|
| EI (Elephant Island) | `Antarc_EI_01_Odontocetes.csv` | 2014-03-05 00:00:00 to 2014-07-17 12:03:50 | 5 min on / 6 min |
| KGI (King George Island) | `Antarc_SSI_01_Odontocetes.csv` | 2015-02-10 00:00:00 to 2016-01-29 13:15:25 | 10 min on / 13.33 min |
| CI (Clarence Island) | `Antarc_EIE_01_Odontocetes.csv` | 2016-02-03 23:59:59 to 2016-12-02 05:15:05 | continuous |

The KGI and CI log file names keep their original internal deployment names;
the mapping is defined once in `config.R` (`SITES$detection_log`). Each log
has one row per detection event with columns `Species Code`, `Call`,
`Start time`, `End time` (`m/d/Y H:M:S`, UTC), and `Comments`.

**Presence = echolocation clicks only.** Only events with `Call` in
`PRESENCE_CALLS` (`"Clicks"`, `"regular-click"`) are used. Sperm whale slow
clicks (44 events) and killer whale HFM signals (3 events, EI) are excluded.

## Scripts

Run from the repository root, in this order:

1. **`build_presence.R`**: reads the three logs, keeps clicks only, gives
   zero-length events a 1-s duration, and builds a 5-min presence grid (a bin
   is present for a species if any detection overlaps it). Daily presence
   minutes = present bins × 5. Daily recording effort (`rec_minutes`) = clock
   minutes inside the recording period (partial first/last days) × duty-cycle
   on-fraction. Writes `outputs/04_acoustic_detection/presence_5min.csv` and
   `presence_daily.csv` (formats in `config.R`, section 5).
2. **`duty_cycle_simulation.R`** (supplementary analysis): subsamples the
   continuously recorded CI log with the EI and KGI duty cycles at 100 random
   phase offsets and reports, per species, the percentage of detection events,
   5-min presence bins, and presence days retained (mean and range). Writes
   `duty_cycle_simulation.csv` (per phase), `duty_cycle_simulation_table.csv`
   (supplementary table), and `duty_cycle_simulation.png/.pdf`.

## References

Baumann-Pickering, S., et al. (2015). [INSERT full reference]

Roch, M. A., et al. (2011). [INSERT full reference]
