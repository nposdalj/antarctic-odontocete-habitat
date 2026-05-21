# Krill Catch Spatial Distribution Maps — CCAMLR Site 756

This script produces six spatial visualizations of Antarctic krill catch from CCAMLR haul data near Elephant Island, King George Island, and Clarence Island (site 756).

---

## Requirements

### R Packages

```r
library(dplyr)
library(ggplot2)
library(sf)
library(marmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ggspatial)
library(units)
```

---

## Input

| File | Description |
|------|-------------|
| `756_USA_2026-01-02.Rds` | CCAMLR data list containing `C1`, `C1_CATCH`, and `OBS_HAUL` tables |

**Tables used:**

| Table | Key columns |
|-------|-------------|
| `C1` | `obs_haul_id`, `latitude_haul_start`, `longitude_haul_start` |
| `OBS_HAUL` | `obs_haul_id`, `greenweight_catch_kg` |

Update the file path in `readRDS()` to match your local directory.

---

## Data Preparation

1. **Krill catch**: `OBS_HAUL` rows with `NA` catch are dropped; haul start coordinates are joined from `C1`.
2. **Land masking**: Hauls that fall on land (via `rnaturalearth`) are removed.
3. **Island reference points**: Elephant Island, King George Island, and Clarence Island are defined manually and plotted as labeled markers where applicable.
4. **Grid averaging**: The study area is divided into a 0.5° × 0.5° grid; mean catch per cell is computed and represented at each cell centroid.
5. **Bathymetry**: 1-minute resolution seafloor depth is fetched from NOAA via `marmap::getNOAA.bathy()` and covers the full spatial extent of the haul data ± 1°.

---

## Plots

### Plot 1 — Raw Haul Locations
Scatter map of individual haul positions. Point size and color both encode `greenweight_catch_kg` using the `plasma` viridis palette. No island labels.

### Plot 2 — Raw Haul Locations with Island Labels
Same as Plot 1 with red island markers and bold labels for Elephant Island, King George Island, and Clarence Island.

### Plot 3 — Grid-Averaged Catch
Hauls are averaged within 0.5° grid cells and plotted at cell centroids. Point size and color encode mean catch. No island labels.

### Plot 4 — Grid-Averaged Catch with Island Labels
Same as Plot 3 with black island markers and labels.

### Plot 5 — Bathymetry + Grid-Averaged Catch
Adds a bathymetric raster (light to dark blue = shallow to deep) beneath the grid-averaged catch points and island labels. Legend sizes are reduced to avoid crowding.

### Plot 6 — Density Hotspots
2D kernel density estimation of haul locations, shown as filled contour polygons (plasma palette) with red contour lines overlaid. Highlights spatial clustering of fishing effort independent of catch weight.

---

## Output

Six `ggplot2` figures rendered to the active graphics device. To save any plot, wrap the relevant block in `ggsave()`:

```r
ggsave("plot5_bathymetry_catch.png", width = 10, height = 8, dpi = 300)
```

---

## Notes

- `dplyr::select` is explicitly assigned at the top of the script to resolve conflicts with other packages (e.g., `raster`).
- Bathymetry is fetched live from NOAA each run; cache it with `saveRDS(bathy, "bathy_756.Rds")` to avoid repeated downloads.
- Island label positions are manually defined. `nudge_y` controls vertical label offset: positive = label above point, negative = below. Clarence Island is nudged downward to avoid overlap with Elephant Island.
- Hauls on land are removed using a spatial intersection test, not a depth filter — this handles island coastlines more precisely.
