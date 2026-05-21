# Antarctic Krill Habitat Analysis — King George Island

This pipeline models the relationship between krill catch (CCAMLR data) and environmental covariates (Copernicus reanalysis/hindcast) near King George Island using Generalized Additive Models (GAMs).

---

## Overview

The script proceeds through five major stages:

1. Load and prepare environmental and krill catch data
2. Check and correct for temporal autocorrelation
3. Screen for collinearity among environmental predictors (VIF)
4. Build and refine the GAM via stepwise variable selection
5. Visualize the final model

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
library(car)
library(mgcv)
library(gratia)
library(patchwork)   # for panel plots
library(ggfortify)   # for ACF confidence intervals
```

---

## Data Inputs

| Variable | Description |
|----------|-------------|
| `env_file` | Daily environmental data CSV (deseasoned/anomaly fields) from Copernicus |
| `catch_file` | CCAMLR haul-level catch `.Rds` file (contains `C1` and `C1_CATCH` tables) |

### Site Configuration

The focal site is **King George Island (KGI)**:

```r
sites <- data.frame(
  name = "KGI",
  lon  = -57.941917,
  lat  = -61.457817
)
buffer_km <- 60   # spatial buffer around site
```

Hauls within 60 km of KGI are retained for analysis.

---

## Pipeline Steps

### 1. Load & Prepare Data

- Environmental data are filtered from `start_date_env` onward and averaged to daily summaries.
- Krill catch (taxon code `"KRI"`) is spatially filtered to hauls within the site buffer, then aggregated to daily totals (`greenweight_caught_kg`).
- Environmental and krill data are merged by date.
- Krill catch is transformed to approximate normality:
  - **Square root**: `value_sqrt = sqrt(value)`
  - **Log**: `value_log = log1p(value)`

### 2. Temporal Autocorrelation

A simple GAM is fit on the raw daily data, and `acf()` is run on its residuals to estimate the autocorrelation length. Data are then **binned** into non-overlapping windows of that length to produce approximately independent observations.

```r
gam_model <- gam(value_sqrt ~ s(EKE_mad, k=4) + s(temp_anom, k=4) + s(salinity_anom, k=4), ...)
acf_res   <- acf(residuals(gam_model), lag.max = 100)
# ACFval = first lag where ACF drops below the confidence interval
model_data_binned <- model_data %>%
  group_by(bin = floor_date(date, unit = my_unit)) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
```

### 3. Collinearity Screening (VIF)

Variance Inflation Factors are computed via `car::vif()` on a GLM containing all candidate predictors. Variables with high VIF are removed iteratively.

**Removed due to collinearity or convergence issues:**
- `zos_mean`
- `sithick_mean`
- `EKE_mean`

**Final candidate set:**

```
mixed_layer_anom, salinity_anom, temp_anom, chla_anom, EKE_mad, siconc_mean, o2_anom
```

### 4. GAM Model Selection

Starting from the full candidate set, the least significant smoother (by p-value) is dropped at each step until only significant terms remain.

**Stepwise removal order:**
1. `mixed_layer_anom`
2. `salinity_anom`
3. `o2_anom`

**Final model (with sea ice):**

```r
gam_model <- gam(
  value_sqrt ~ s(chla_anom, k=4) + s(EKE_mad, k=4) + s(siconc_mean, k=4),
  data   = model_data_binned,
  family = gaussian()
)
```

An alternative model excluding `siconc_mean` is also evaluated for comparison.

### 5. Visualization

- **Single-variable plot**: Predicted partial effect of `EKE_mad` on `sqrt(krill catch)` with 95% confidence ribbon and rug plot. Annotated with p-value and deviance explained.
- **Panel plot**: Side-by-side partial effect plots for `chla_anom` and `siconc_mean`, with a shared title showing total deviance explained.

---

## Environmental Variables

| Variable | Description |
|----------|-------------|
| `chla_anom` | Chlorophyll-a anomaly (mg/m³) |
| `EKE_mad` | Mean absolute deviation of Eddy Kinetic Energy (cm²/s²) |
| `siconc_mean` | Mean sea ice concentration |
| `temp_anom` | Sea surface temperature anomaly |
| `salinity_anom` | Salinity anomaly |
| `o2_anom` | Dissolved oxygen anomaly |
| `mixed_layer_anom` | Mixed layer depth anomaly |
| `zos_mean` | Mean sea surface height |

---

## Output

- Model summary (`summary(gam_model)`) including deviance explained and smoother significance
- Diagnostic plots via `gratia::draw()` and base `plot.gam()`
- `ggplot2` figures of partial effects for key predictors

---

## Notes

- `siconc_mean` NAs are set to 0 prior to modeling (ice-free assumption for missing values).
- All smoothers use `k = 4` (3 effective degrees of freedom) to avoid overfitting given the limited post-binning sample size.
- The response variable is `value_sqrt` (square root of kg krill caught) throughout the GAM.

