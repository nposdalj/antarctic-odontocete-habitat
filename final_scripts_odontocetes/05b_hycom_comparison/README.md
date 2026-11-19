# 05b_hycom_comparison — Supplementary Figs. S3 and S4: HYCOM vs Copernicus

The habitat models use water-column temperature and salinity from the
Copernicus global reanalysis (GLOBAL_MULTIYEAR_PHY_001_030, 1/12°, daily). This
supplement checks how much those values depend on the choice of ocean model.
It compares Copernicus with:

- **Aqua MODIS satellite SST** at the surface, and
- **HYCOM GOFS 3.1**, an independent ocean model, at the surface and at each
  species-relevant depth.

Each product is averaged over the same ±20 km box (the main analysis box; downloads cover ±40 km and are cropped in `03`) around each recorder (EI,
KGI, CI) for each day of the recording period. The products are then joined
by date.

## Scripts (run in order from the repository root)

| Script | What it does | Output |
|---|---|---|
| `01_download_hycom.R` | Downloads daily HYCOM GLBv0.08 temperature and salinity from 0 to 2000 m over each box via OPeNDAP. You can re-run it after an interruption; it skips days that are already on disk. | `outputs/downloads/hycom/<site>/hycom_<site>_<date>.nc`, `hycom_download_log.csv` |
| `02_download_modis_sst.R` | Downloads Aqua MODIS daily SST (`sstMasked`) over each box from NOAA CoastWatch ERDDAP. | `outputs/downloads/modis/modis_sst_<site>.nc` |
| `03_compare_hycom_copernicus.R` | Computes daily box means, joins the products by date, calculates the agreement metrics and draws the figure. | `outputs/05b_hycom_comparison/` (see below) |

`03` also needs the Copernicus daily box means by model level
(`outputs/05_environmental_data/copernicus_phy_by_depth.csv`), which the
`05_environmental_data` scripts produce.

### Outputs of `03`

- `comparison_metrics.csv`: one row per comparison × variable × site ×
  species depth × value type. Columns: `comparison`, `variable`, `site`,
  `nominal_depth_m`, `depth_m` (the Copernicus level used), `values` (`raw` or
  `anomaly`), `n_days`, `bias` (mean of first product − second product),
  `rmse`, `r` (Pearson).
  - (a) Surface temperature: HYCOM vs MODIS and Copernicus vs MODIS. Both are
    scored on the same days, i.e. days with all three products.
  - (b) Temperature and salinity at every species depth: HYCOM vs Copernicus.
- `hycom_copernicus_profile.csv`: raw-value bias, RMSE and r at every
  Copernicus level down to 1684 m. It also gives the fraction of box cells
  that hold water in each model.
- `hycom_daily_box_means.csv`, `modis_daily_box_means.csv`: the daily box
  means that go into the comparison.
- `FigS3_hycom_copernicus.png` (600 dpi) and `.pdf`:
  - (a) surface temperature time series for each site;
  - (b) HYCOM − Copernicus bias and RMSE profiles, with species depths marked.

Metrics are computed on raw daily values and on deseasoned anomalies. Each
product's anomaly is its value minus its own site × calendar-month mean over
the matched days. The bias of anomalies is zero by construction, so it is
reported as NA. No hypothesis tests are applied. The earlier version used
Shapiro-Wilk, paired t, Wilcoxon and Diebold-Mariano tests; these are dropped
because daily values are strongly autocorrelated and the goal is to describe
agreement.

## Why the data were downloaded again (box-coverage problem)

The original comparison (`Environmental Data/Cop_vs_HYCOM.R`) used HYCOM `.mat`
files (`Environmental Data/HYCOM/Data/`) and MODIS files
(`Environmental Data/SST_ERDDAP/`). Those requests used ±0.5° longitude
windows that do not match the ±40 km boxes used for Copernicus:

| Site | Site longitude | ±40 km box longitude | Original HYCOM/MODIS window |
|---|---|---|---|
| EI  | −55.95 | −56.69 to −55.22 | −56.56 to −55.60 (partial) |
| KGI | −57.94 | −58.69 to −57.19 | −57.44 to −56.48 (**does not contain the site**) |
| CI  | −53.48 | −54.23 to −52.74 | −54.00 to −52.96 (partial) |

The original script also combined the series with `cbind()` by row position
rather than by date. This misaligns the series whenever one product is missing
a day. Both problems are fixed here: every product is re-extracted over the
config box and joined by date.

## Data sources

- **HYCOM GOFS 3.1**, GLBv0.08 (1/12° × 1/25°, 40 z-levels, 3-hourly),
  accessed through `https://tds.hycom.org/thredds/dodsC/GLBv0.08/`. No single
  experiment covers 2014–2016, so the experiments from the original download
  script (`Environmental Data/HYCOM/get_hycom.m`) are used:
  - EI: `expt_53.X` (reanalysis; `data/2014`)
  - KGI: `expt_56.3` (analysis, 2014-07-01 to 2016-09-30)
  - CI: `expt_56.3` until 2016-09-30, then `expt_57.2` (analysis, 2016-05-01
    to 2017-01-31)

  One snapshot per day (the time step closest to 12:00 UTC) is used. HYCOM's
  NetcdfSubset service (`ncss.hycom.org`) timed out or failed in testing, so
  OPeNDAP is used instead.
- **Aqua MODIS SST**: NOAA CoastWatch ERDDAP dataset
  `erdMH1sstd1dayR20190SQ_Lon0360`, field `sstMasked`. This is the daily,
  daytime, 4 km, quality-masked product (R2019.0). MODIS is a daytime skin
  temperature and is missing under cloud and sea ice. A day is used only if
  at least 20 % of the box's ocean pixels are valid.
- **Copernicus**: GLOBAL_MULTIYEAR_PHY_001_030 daily box means from
  `05_environmental_data`.

Caveats:

- HYCOM and Copernicus both assimilate satellite SST, so the surface
  comparison with MODIS is not fully independent.
- The depth comparison measures agreement between two models, not accuracy.

## In-situ data: limitation

No in-situ validation was possible. An availability check was run on
2026-11-16. It searched the World Ocean Database for 2014–2016, all probe
types (NCEI OPeNDAP), and Argo profiles (Ifremer ERDDAP). It found the
following profiles inside the ±40 km boxes during the recording periods:

| Site | Profiles in box during recording |
|---|---|
| EI  | 0 |
| KGI | 5 (2 Argo, 3 XBT) |
| CI  | 0 |

The wider region (60–63° S, 52–60° W) has about 900 profiles over the same
years. About 780 are from animal-borne (seal-tag) CTDs, about 73 from Argo
floats and about 55 from ship CTD or bottle casts. Five in-box profiles at a
single site are too few for a match-up validation, and the regional profiles
lie outside the boxes. The model comparison above is therefore the only check
of the hydrographic predictors.

## Items to confirm

Search the scripts for `[CONFIRM]`:

- which HYCOM experiment the original `.mat` files came from;
- the 12:00 UTC daily snapshot (instead of a daily mean of the 3-hourly steps);
- interpolating HYCOM to the Copernicus levels (instead of using the nearest
  HYCOM level);
- the 20 % MODIS coverage threshold;
- the anomaly definition. The main analysis deseasons against a 10-year
  climatology smoothed with a periodic spline. There is no 10-year HYCOM
  record here, so anomalies use within-period monthly means.
