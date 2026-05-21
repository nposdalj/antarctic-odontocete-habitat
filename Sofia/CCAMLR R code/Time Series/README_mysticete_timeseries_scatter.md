# Mysticete Detections & Environmental Time Series / Krill–Environment Scatter Plots

This README covers four related scripts that visualize mysticete acoustic detections alongside environmental and krill catch data from CCAMLR site 756.

---

## Scripts

| Script | Description |
|--------|-------------|
| `timeseries_by_variable.R` | Multi-panel time series; bars colored by variable type |
| `timeseries_by_site.R` | Multi-panel time series; bars and lines colored by site |
| `scatter_single.R` | Scatter plot of krill catch vs one environmental variable with linear model |
| `scatter_panel.R` | Faceted scatter plots of krill catch vs all environmental variables |

---

## Requirements

### R Packages

```r
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(readxl)
library(sf)
library(units)
```

---

## Inputs

| File | Description |
|------|-------------|
| `allData_fin.csv` | Daily fin whale (`Bp`) detections and environmental data |
| `allData_blueZ (5).csv` | Daily blue whale (`Bm`) detections and environmental data |
| `allData_humpback (2).csv` | Daily humpback whale (`Mn`) detections and environmental data |
| `KGI_deseasoned.csv` / `CI_deseasoned.csv` | Deseasoned daily environmental data for a single site |
| `756_USA_2026-01-02.Rds` / `756_USA.Rds` | CCAMLR data list containing `C1` and `C1_CATCH` tables |

**CCAMLR tables used:**

| Table | Key columns |
|-------|-------------|
| `C1` | `c1_id`, `obs_haul_id`, `latitude_haul_start`, `longitude_haul_start`, `date_catchperiod_start` |
| `C1_CATCH` | `c1_id`, `taxon_code`, `greenweight_caught_kg` |

---

## Configuration

### Species (time series scripts)

```r
species <- "Mn"   # Options: "Bp" (fin), "Bm" (blue), "Mn" (humpback)
```

The correct CSV is loaded automatically based on this value.

### Environmental variable (scatter_single.R only)

```r
env_var <- "mlayer_0"
# Options: temperature_0, salinity_0, chla_0, productivity_0, mlayer_0
```

### Spatial buffer (timeseries_by_site.R and scatter scripts)

Krill hauls are spatially filtered to hauls within a buffer around each site:

```r
buffer_km     <- 100
antarctic_crs <- 3031   # Antarctic Polar Stereographic CRS (EPSG:3031)
```

Sites used:

| Site | Lon | Lat |
|------|-----|-----|
| EI (Elephant Island) | −55.954 | −60.887 |
| KGI (King George Island) | −57.942 | −61.458 |
| CI (Clarence Island) | −53.483 | −61.252 |

---

## Environmental Variables

| Column | Description |
|--------|-------------|
| `temperature_0` | Sea surface temperature |
| `salinity_0` | Sea surface salinity |
| `chla_0` | Chlorophyll-a |
| `productivity_0` | Net primary productivity |
| `mlayer_0` / `mixed_layer_anom` | Mixed layer depth |
| `EKE_0` | Eddy kinetic energy |
| `o2_0` | Dissolved oxygen |

In the deseasoned CSVs, raw column names differ (`sst_none`, `sss_none`, `chla_none`, `npp_none`, `mlayer_none`) and are renamed during summarisation.

---

## Script Details

### Script 1 — `timeseries_by_variable.R`

Multi-panel time series from 2014 onward. Data are **not** spatially filtered by buffer — krill is aggregated across all site 756 hauls. Bars and lines are colored by **variable type** (species bars in purple, krill bars in pink, env lines in fixed colors).

**Panels (top to bottom):** species detections → krill catch → productivity → chlorophyll → EKE → oxygen → salinity → temperature

**Grey shading:** periods when all sites have zero acoustic recording effort (acoustic panels) or zero fishing hauls (krill panel).

**Grey lines:** recording/fishing effort rescaled to each panel's y-axis range, shown per site with different linetypes (KGI solid, EI dashed, CI dotted).

---

### Script 2 — `timeseries_by_site.R`

Same structure as Script 1 but krill catch is **spatially filtered** to hauls within 100 km of EI, KGI, or CI. Bars and lines are colored by **site** (KGI blue, EI green, CI red). Effort linetype legend is hidden to reduce clutter.

---

### Script 3 — `scatter_single.R`

Scatter plot of daily krill catch (kg) vs a single deseasoned environmental variable at KGI. Fits a linear model and annotates the plot with R² and p-value.

**Output:** one `ggplot2` scatter plot with:
- Black points (daily observations)
- Linear regression line with grey 95% confidence ribbon
- R² and p-value label (top-right)

---

### Script 4 — `scatter_panel.R`

Same data preparation as Script 3 but plots **all five environmental variables** simultaneously as a faceted panel. Each facet has an independent y-axis. No R² annotation — intended for visual comparison across variables.

**Output:** one `ggplot2` faceted plot with five panels: Temperature, Salinity, Chlorophyll, Productivity, Mixed Layer.

> **Note:** Scripts 3 and 4 share identical data prep code. The site label in the `sites` data frame reads `"KGI"` but loads `CI_deseasoned.csv` — verify the site name and file path match your intended analysis area before running.

---

## Output

All four scripts produce `ggplot2` figures rendered to the active graphics device. To save:

```r
ggsave("output_plot.png", width = 12, height = 10, dpi = 300)
```
