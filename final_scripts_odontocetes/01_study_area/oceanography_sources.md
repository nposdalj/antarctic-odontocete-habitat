# Sources for schematic currents and water masses (study-area map)

Compiled 2026-11-18. All DOIs were resolved via Crossref/DataCite. Waypoints are **approximate eye-readings** (±0.2–0.5°) from published figures, suitable only for schematic arrows; for publication, digitize from the original figures (e.g., georeference in QGIS). None of these pathways exist as GIS files except the Orsi et al. (1995) fronts (`orsifronts` R package).

## Best single figure to adapt
**Damini et al. (2025), Fig. 1** (modified from Damini et al. 2022): CDW intrusions, Bransfield Current + recirculation, Antarctic Coastal Current (ACoC), Dense Shelf Water via Antarctic Sound, Antarctic Slope Front + bifurcation eddy east of Clarence Island, Peninsula Front, SACCF, SB.
- Damini B.Y., Brum A.L., Hall R.A., Dotto T.S., Azevedo J.L.L., Heywood K.J., Mata M.M., Garcia C.A.E., Kerr R. (2025) *Deep-Sea Res. I* 222, 104516. https://doi.org/10.1016/j.dsr.2025.104516 (accepted MS: nora.nerc.ac.uk/id/eprint/539533)
- Damini B.Y., Kerr R., Dotto T.S., Mata M.M. (2022) *Deep-Sea Res. I* 179, 103667. https://doi.org/10.1016/j.dsr.2021.103667

## Feature waypoints (lon, lat; approximate)

| Feature | Waypoints | Sources (figure) |
|---|---|---|
| mCDW on-shelf north of KGI | (−58.8,−61.55)→(−58.5,−61.85); (−57.6,−61.35)→(−57.4,−61.6) | Zhou et al. 2013 (Figs 1, 4–6, 12); Damini et al. 2025 (Fig 1) |
| mCDW into eastern basin (KGI–EI) | (−57.0,−61.4)→(−56.9,−61.8)→(−57.3,−62.1) | Damini et al. 2025 (Fig 1) |
| Shelf-water export north of EI | (−55.6,−61.0)→(−55.5,−60.5)→(−55.0,−60.0) | Zhou et al. 2013 |
| Bransfield Current (main) | (−61.0,−62.85)→(−59.8,−62.5)→(−58.7,−62.25)→(−57.8,−62.05)→(−56.8,−61.8)→(−55.8,−61.45)→(−54.7,−61.15)→(−54.8,−60.8) | Frey et al. 2023 (Fig 1); Veny et al. 2024 (Fig 1a); Damini et al. 2025; Sangrà et al. 2011, 2017; Veny et al. 2022 |
| BC recirculation (dashed) | (−56.5,−61.5)→(−58.0,−61.6)→(−59.5,−62.0)→(−61.0,−62.4) | Damini et al. 2025 (Fig 1); Dotto et al. 2021 (Fig 1b) |
| TWW / Antarctic Coastal Current into Bransfield | (−52.5,−63.5)→(−53.0,−62.9)→(−54.0,−62.5)→(−55.0,−62.4)→(−56.3,−62.65)→(−57.5,−62.95)→(−58.5,−63.2)→(−59.5,−63.45) | Thompson et al. 2009 (Figs 2–3); Heywood et al. 2004; Damini et al. 2025 (entry ~62.40°S, 55.00°W); Frey et al. 2023 |
| Antarctic Sound branch ("possible", dashed) | (−56.4,−63.5)→(−56.8,−63.35)→(−57.1,−63.2) | Frey et al. 2023 (Fig 1, after Sangrà et al. 2011) |
| Antarctic Slope Front | (−52.0,−63.8)→(−52.4,−62.9)→(−52.5,−62.1)→(−53.1,−61.4)→(−51.5,−60.9)→(−51.0,−60.0) | Thompson et al. 2009 (Fig 3); Damini et al. 2025; Dotto et al. 2021 |
| ASF westward branch | (−52.5,−62.4)→(−54.8,−62.5) | Thompson et al. 2009 (Fig 3) |
| ASF standing eddy | centre ≈ (−53.2,−61.7) | Thompson et al. 2009; Damini et al. 2025 |
| Peninsula Front | (−61.0,−63.2)→(−59.8,−62.95)→(−58.5,−62.7)→(−57.2,−62.45)→(−56.0,−62.3) | Veny et al. 2024 (Fig 1a,b); Dotto et al. 2021 (Fig 1b); Damini et al. 2025 |
| Bransfield Front (distinct from PF) | (−61.0,−62.75)→(−59.5,−62.4)→(−58.3,−62.15)→(−57.3,−61.95) | Sangrà et al. 2011; Moffat & Meredith 2018 |
| Polar Front | outside map: ~58.6°S at 61°W, 56.5°S at 51°W | Orsi et al. 1995 (`orsifronts`) |

Context: EI sits just downstream of the Shackleton Gap (Shackleton Transverse Ridge ≈ (−57.5,−60.0) to (−56.3,−60.8)); CI lies in the ASF bifurcation zone; KGI sits on the north South Shetland slope where Moffat & Meredith (2018) draw a southwestward Polar Slope Current (200–600 m).

## Water-mass acronyms
- **TBW**: "Transitional Zonal Water with Bellingshausen Sea influence" (Tokarczyk 1987; García et al. 1994; Sangrà et al. 2011, 2017). Θ > −0.4 °C, S < 34.45 (Veny et al. 2024 citing Sangrà et al. 2017); Moffat & Meredith 2018 use S < 34.35 for the equivalent "Shelf Slope Water".
- **TWW**: "Transitional Zonal Water with Weddell Sea influence". Θ < −0.4 °C, S > 34.45, vertically homogeneous. TBW/TWW separated by σθ ≈ 27.64 (Veny et al. 2024).
- **BSW**: no authoritative acronym definition found — avoid or define explicitly.
- **WSW**: "Weddell Sea shelf water" (Whitworth et al. 1994); others use HSSW/LSSW or Dense Shelf Water.
- **Peninsula Front ≠ Bransfield Front.** PF: shallow (0–100 m), ~20–30 km off the peninsula, TBW|TWW boundary. BF: subsurface (50–400 m), 10–30 km off the South Shetland coast.

## References (verified DOIs)
- Dotto T.S., Mata M.M., Kerr R., Garcia C.A.E. (2021) *ESSD* 13, 671–696. https://doi.org/10.5194/essd-13-671-2021 — NAPv1.0 gridded T/S/O₂ (data: https://doi.org/10.5281/zenodo.4420006); can contour PF and mCDW from data.
- Frey D. et al. (2023) *Front. Mar. Sci.* 10, 1111541. https://doi.org/10.3389/fmars.2023.1111541 (CC-BY)
- Freeman N.M., Lovenduski N.S. (2016) *ESSD* 8, 191–198. https://doi.org/10.5194/essd-8-191-2016
- Gordon A.L. et al. (2000) *JGR* 105, 11337–11346. https://doi.org/10.1029/2000JC900030
- Heywood K.J., Naveira Garabato A.C., Stevens D.P., Muench R.D. (2004) *JGR* 109, C06021. https://doi.org/10.1029/2003JC002053
- Huneke W.G.C., Huhn O., Schröeder M. (2016) *Polar Biol.* 39, 789–798. https://doi.org/10.1007/s00300-016-1936-8
- Moffat C., Meredith M.P. (2018) *Phil. Trans. R. Soc. A* 376, 20170164. https://doi.org/10.1098/rsta.2017.0164
- Orsi A.H., Whitworth T., Nowlin W.D. (1995) *Deep-Sea Res. I* 42, 641–673.
- Park Y.-H., Durand I. (2019) SEANOE fronts. https://doi.org/10.17882/59800
- Sangrà P. et al. (2011) *Deep-Sea Res. I* 58, 390–402. https://doi.org/10.1016/j.dsr.2011.01.011
- Sangrà P. et al. (2017) *Deep-Sea Res. I* 119, 1–15. https://doi.org/10.1016/j.dsr.2016.11.003
- Thompson A.F., Heywood K.J., Thorpe S.E., Renner A.H.H., Trasviña A. (2009) *J. Phys. Oceanogr.* 39, 3–26. https://doi.org/10.1175/2008JPO3995.1
- Thompson A.F., Heywood K.J. (2008) *Deep-Sea Res. I* 55, 1229–1251. https://doi.org/10.1016/j.dsr.2008.06.001
- van Caspel M., Hellmer H.H., Mata M.M. (2018) *Deep-Sea Res. II* 149, 25–30. https://doi.org/10.1016/j.dsr2.2017.09.006
- Veny M. et al. (2022) *Prog. Oceanogr.* 204, 102795. https://doi.org/10.1016/j.pocean.2022.102795
- Veny M. et al. (2024) *Ocean Sci.* 20, 389–415. https://doi.org/10.5194/os-20-389-2024 (CC-BY)
- Whitworth T. et al. (1994) *Deep-Sea Res. I* 41, 629–641. https://doi.org/10.1016/0967-0637(94)90046-9
- Zhou M. et al. (2013) *Deep-Sea Res. II* 90, 4–14. https://doi.org/10.1016/j.dsr2.2013.03.041
- Zhou M., Niiler P.P., Zhu Y., Dorland R.D. (2006) *Deep-Sea Res. I* 53, 1244–1252. https://doi.org/10.1016/j.dsr.2006.04.003

Not verified/read: Sangrà 2011/2017 figure numbers; Heywood 2004 Fig 6; Tokarczyk 1987; García et al. 1994.
