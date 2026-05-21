# Sea Ice Concentration vs Krill Catch & Weekly Cetacean Presence

This file documents two visualization scripts: one producing scatter plots of sea ice concentration against krill catch for two management areas (APDPE and APEI), and one producing a weekly cetacean presence bar chart overlaid with sea ice thickness.

---

## Script 1: SIC vs Krill Catch Scatterplots

### Overview

Generates two separate scatter plots — one for the **APDPE** area and one for the **APEI** area — showing the relationship between monthly mean sea ice concentration (%) and krill catch weight, with a LOESS trend line overlaid.

### Inputs

| File | Area | Description |
|------|------|-------------|
| `APDPE_krill_env_monthly.csv` | APDPE | Monthly krill catch and environmental data |
| `APEI_krill_env_monthly.csv` | APEI | Monthly krill catch and environmental data |

**Required columns in each CSV:**

| Column | Description |
|--------|-------------|
| `sice_conc_mean_month_mean` | Monthly mean sea ice concentration (0–1 fraction) |
| `krill_weight` | Krill catch weight; rows with `NA` are excluded |

### What the Script Does

1. Loads the CSV and converts sea ice concentration to percentage (`× 100`).
2. Removes rows where `krill_weight` is `NA`.
3. Assigns a color to each point along a gradient keyed to sea ice concentration value:
   - **APDPE**: turquoise (`#40E0D0`) → coral (`#FF6F61`)
   - **APEI**: medium purple (`#9F79EE`) → pale turquoise (`#AFEEEE`)
4. Plots `ice_pct` (x) vs `krill_weight` (y) as filled circles (`pch = 19`, `cex = 1.5`).
5. Fits a LOESS smoother and draws it as a thick line:
   - **APDPE**: `royalblue3`
   - **APEI**: `royalblue4`

### Output

Two base-R scatter plots (one per area), rendered to the active graphics device. Each includes:
- Gradient-colored points scaled to SIC value
- LOESS trend line
- Axis labels and title

---

## Script 2: Weekly Cetacean Presence & Sea Ice Thickness

### Overview

Produces a dual-axis bar chart showing **weekly acoustic detection counts** for a selected cetacean species at a selected site, overlaid with a **sea ice thickness** time series.

### Input

| File | Description |
|------|-------------|
| `dailyDetections.csv` | Daily detection and environmental data per site |

**Required columns:**

| Column | Description |
|--------|-------------|
| `Day` | Date string in `%m/%d/%Y` format |
| `Site` | Site code (e.g., `CI`, `EI`, `KGI`) |
| `[species]` | Binary daily detection column (0 or 1); column name set by `species` variable |
| Column 8 | Sea ice thickness (m); assumed to be the 8th column of the data frame |

### Configuration

```r
species <- "Pm"    # Column name for the target species (e.g., Pm = sperm whale)
site    <- "CI"    # Site code to filter: CI, EI, or KGI
```

### What the Script Does

1. Filters `dailyDetections` to the selected site.
2. Assigns each day to its week-start (previous Sunday).
3. Aggregates weekly:
   - **Presence**: sum of daily detections per week (range 0–7)
   - **Sea ice thickness**: mean of column 8 per week
4. Plots a barplot of weekly presence (left y-axis, 0–7).
5. Scales sea ice thickness to the same 0–7 range and overlays it as a line (right y-axis shows original thickness in meters).

### Output

A single dual-axis base-R plot:
- **Bars** (light pink): weekly detection count (0–7 days present)
- **Line** (deep pink): weekly mean sea ice thickness, scaled to bar axis
- **Left y-axis**: weekly presence (days)
- **Right y-axis**: sea ice thickness (m)
- **X-axis**: week-start dates, labeled every other week at 45°

### Notes

- Sea ice thickness is read from **column 8 by index**, not by name. If the column order changes, update `seaice_column <- 8` accordingly.
- The LOESS smoother in Script 1 may behave poorly with sparse data at the tails of the ice concentration range; inspect the trend line visually before interpreting.
- The APEI block appears twice in the source — the second copy is a duplicate and can be removed.
