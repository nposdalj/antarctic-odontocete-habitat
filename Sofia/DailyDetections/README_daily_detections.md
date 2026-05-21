# Daily Cetacean Detection Presence Plots

This script generates simple time series plots of daily acoustic detection presence (0/1) for a selected species and site.

---

## Input

| File | Description |
|------|-------------|
| `dailyDetections.csv` | Daily binary detection records per species and site |

**Required columns:**

| Column | Description |
|--------|-------------|
| `Day` | Date string in `%m/%d/%Y` format |
| `Site` | Site code (e.g., `CI`, `KGI`, `EI`) |
| `[species]` | Binary detection column (0 = absent, 1 = present); one column per species |

---

## Configuration

Two variables control which data are plotted:

```r
species <- "Pm"   # Species column to plot — options: "Pm", "Oo", "BW29", "Gm", etc.
site    <- "CI"   # Site to filter  — options: CI, KGI, EI
```

Change these before each `plot()` call to switch species or site.

---

## What the Script Does

1. Reads `dailyDetections.csv` and parses the `Day` column as a date.
2. Filters rows to the selected `site`.
3. Plots daily presence as a vertical line (`type = "h"`) over time.

Repeat the configuration + plot block for each species or site combination you want to visualize. The current script produces two plots:

| Plot | Species | Site |
|------|---------|------|
| 1 | `Pm` (sperm whale) | Clarence Island (CI) |
| 2 | `Oo` (killer whale) | Clarence Island (CI) |

---

## Output

One base-R presence plot per block, rendered to the active graphics device:
- **X-axis**: date
- **Y-axis**: presence (0 or 1)
- **Plot type**: vertical line (`type = "h"`)

---

## Notes

- Update the file path in `read.csv()` to wherever `dailyDetections.csv` is saved locally.
- To add more species or sites, copy the configuration + `plot()` block and update `species` and `site`.
- Species codes used in this dataset include at minimum: `Pm`, `Oo`, `BW29`, `Gm` — check column headers in the CSV for the full list.
