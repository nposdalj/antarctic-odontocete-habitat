# 2.1 Study Area — methods draft

*Drafted from `01_study_area/site_map.R`, `config.R`, and `data/DataSummary.xlsx`. Flags: `[INSERT]` = information not in the code; `[CONFIRM]` = value to verify.*

---

Passive acoustic data were collected at three sites along the northern shelf slope of the South Shetland Islands, Antarctica: Elephant Island (EI; 60°53.214′S, 55°57.238′W; 762 m depth), King George Island (KGI; 61°27.469′S, 57°56.515′W; 768 m), and Clarence Island (CI; 61°15.112′S, 53°29.006′W; 1033 m) (Fig. 1).

The South Shetland Islands lie at the confluence of water masses from the Drake Passage, the Bellingshausen Sea, and the Weddell Sea (Fig. 1). To the north, the Southern Antarctic Circumpolar Current Front and the Southern Boundary of the Antarctic Circumpolar Current (ACC) follow the continental slope (Orsi et al. 1995), and the Southern Boundary passes close to the EI and KGI sites. Upper Circumpolar Deep Water from the ACC intrudes onto the northern South Shetland shelf as modified Circumpolar Deep Water (mCDW), particularly west of the Shackleton Transverse Ridge, whereas shelf waters are exported offshore north of Elephant Island, just downstream of the EI site (Zhou et al. 2013, Damini et al. 2025). South of the islands, the Bransfield Current flows northeastward along the southern slope of the South Shetland Islands, carrying relatively warm and fresh Transitional Zonal Water with Bellingshausen Sea influence (TBW; Sangrà et al. 2011, 2017, Veny et al. 2022), with part of the flow recirculating southwestward along the northern side of the islands (Dotto et al. 2021, Damini et al. 2025). Colder and saltier Transitional Zonal Water with Weddell Sea influence (TWW) enters Bransfield Strait from the east via the Antarctic Coastal Current and the Antarctic Slope Front around the tip of the Antarctic Peninsula (Heywood et al. 2004, Thompson et al. 2009) and is separated from TBW by the Peninsula Front (Sangrà et al. 2011, Veny et al. 2024). East of Elephant Island, the Antarctic Slope Front bifurcates over the South Scotia Ridge near the CI site, forming a westward branch and a quasi-stationary eddy (Thompson et al. 2009). The three sites therefore sample contrasting oceanographic regimes: EI near the shelf-water export pathway at the Shackleton Gap, KGI on the northern South Shetland slope beneath the Southern Boundary of the ACC, and CI under stronger Weddell Sea influence [CONFIRM characterization of each site].

High-frequency Acoustic Recording Packages (HARPs; Wiggins & Hildebrand 2007) were deployed sequentially at the three sites, sampling at 200 kHz. Recordings spanned 5 March – 17 July 2014 at EI (134.5 days), 10 February 2015 – 29 January 2016 at KGI (353.6 days), and 3 February – 2 December 2016 at CI (302.2 days) (Table 1). The EI and KGI recorders were duty-cycled (5 min on every 6 min and 10 min on every 13.33 min, respectively), whereas the CI recorder recorded continuously.

Bathymetry for the study-area map (Fig. 1) was obtained from the Global Multi-Resolution Topography synthesis (GMRT v4.4.1; Ryan et al. 2009; accessed 11 May 2026 via GMRT MapTool) and contoured at 500-m intervals. Regional circulation (Fig. 1) was illustrated with long-term (2009–2021) mean surface currents from the Copernicus Marine Service global ocean physics reanalysis (GLOBAL_MULTIYEAR_PHY_001_030, dataset version 202311, monthly means at 0.5 m; accessed 17 November 2026) and the positions of the Southern Antarctic Circumpolar Current Front and the Southern Boundary of the Antarctic Circumpolar Current from Orsi et al. (1995).

---

**Figure 1 (draft legend).** Study area and regional circulation around the South Shetland Islands, Antarctica. Circles mark the High-frequency Acoustic Recording Package (HARP) sites at Elephant Island (EI; 2014), King George Island (KGI; 2015–2016), and Clarence Island (CI; 2016), shown on GMRT bathymetry (v4.4.1; Ryan et al. 2009) with 500-m isobaths. Black arrows show mean surface currents for 2009–2021 from the Copernicus Marine Service global ocean physics reanalysis (GLORYS12V1). Black dashed and dot-dashed lines mark the Southern Antarctic Circumpolar Current Front (SACCF) and the Southern Boundary of the Antarctic Circumpolar Current (SB), respectively (Orsi et al. 1995). Colored lines show schematic pathways digitized from Damini et al. (2025): intrusion of modified Circumpolar Deep Water (mCDW); the Bransfield Current (BC), carrying Transitional Zonal Water with Bellingshausen Sea influence (TBW), and its recirculation; the Antarctic Coastal Current (ACoC), carrying Transitional Zonal Water with Weddell Sea influence (TWW); the Antarctic Slope Front (ASF) and associated eddy; and the Peninsula Front (PenF). The inset shows the location of the study area (red box) relative to the Antarctic Peninsula and South America.

---

**Table 1 (draft).** HARP deployments.

| Site | Latitude | Longitude | Depth (m) | Recording start (UTC) | Recording end (UTC) | Duration (d) | Duty cycle |
|---|---|---|---|---|---|---|---|
| EI | 60°53.214′S | 55°57.238′W | 762 | 2014-03-05 00:00 | 2014-07-17 12:03 | 134.5 | 5 of 6 min |
| KGI | 61°27.469′S | 57°56.515′W | 768 | 2015-02-10 00:00 | 2016-01-29 13:15 | 353.6 | 10 of 13.33 min |
| CI | 61°15.112′S | 53°29.006′W | 1033 | 2016-02-03 23:59 | 2016-12-02 05:15 | 302.2 | Continuous |

---

### Notes for Natalie
- **Oceanographic-setting paragraph** (new, 2026-11-18): every statement is tied to a source verified in `oceanography_sources.md`. Water-mass names follow the source literature (TBW/TWW = Transitional Zonal Water with Bellingshausen/Weddell Sea influence). Points to check: (i) the per-site characterization in the last sentence is my synthesis from the sources and the map — confirm it matches your understanding; (ii) "relatively warm and fresh" (TBW) vs "colder and saltier" (TWW) follows Veny et al. (2024): TBW Θ > −0.4 °C, S < 34.45; TWW Θ < −0.4 °C, S > 34.45; (iii) Sangrà et al. 2011/2017 full texts were not read (citations verified; content from abstracts and from later papers citing them). Your Results (2.2/3.x) show warm deep water at KGI and cold water at CI — that fits this framing but belongs in Results, not here.
- **EI end date:** `DataSummary.xlsx` lists the EI data end date as 2014-07-14, but its own recording duration (134.5 d) and the detection log (detections until 2014-07-17 10:54) both indicate 2014-07-17 12:03:50. The draft uses 07-17; please correct the spreadsheet or confirm.
- **Depths** are from `DataSummary.xlsx` (762, 768, 1033 m); the earlier briefing rounded them to 760, 760, 1030 m.
- **Sample rate:** 200 kHz (confirmed by Natalie).
- **GMRT:** the grid file records only "Created with GMRT MapTool"; no version. The file was created 2026-05-11 (download date). GMRT lists its current version and how to cite it on gmrt.marine-geo.org; use the version that was current on 11 May 2026. Verify the Ryan et al. (2009) citation.
- **Currents:** mean surface currents and Orsi et al. (1995) fronts are drawn in Figure 1 (`study_area_currents_map.R`). The Southern Boundary of the ACC passes close to EI and KGI; CI lies south of it.
- **Duty cycles** stay in this section (study area and data collection), as in Table 1.
- **Figure 1** is `study_area_currents_map_schematic.png`, produced by `study_area_currents_map.R` (GMRT bathymetry, sites, mean surface currents, Orsi et al. 1995 fronts, and pathways digitized from Damini et al. 2025). Main text only; nothing in the supplement.
