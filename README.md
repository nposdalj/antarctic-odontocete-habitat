## Overview
Habitat modeling code for odontocete detections at three HARP sites off the Antarctic Peninsula (near the South Shetland Islands) from 2014 to 2016. Detected odontocetes include sperm whales, killer whales, 
long-finned pilot whales, southern bottlenose whales (BW29), and Gray's and strap-toothed whales (BW37 and BW58). Detected baleen whales will be habitat modeled separately. Models used are generalized additive models (GAMs).

## Procedure
### Step 1: Create timeseries of species detections.
*Input:* Species detections data (start/end time format for each encounter)
*Output:* Timeseries of species detections
Loaded detections (with start and end times) are in three files in the "Data" folder titled Antarc_*SITECODE*_01_Odontocetes.csv. Use TimeseriesPlots.R to create a site-by-site, species-by-species timeseries of detections.
Options include total hours of clicks detected per day and number of days with binary daily presence per week. Save all site detections to Antarc_Odontocetes.csv

### Step 2: Load environmental data and create timeseries.
*Input:* Environmental data (various sources, netCDF format)
*Output:* Environmental data timeseries & .csv sheets for each site
#### First, find data sources.
Look for satellite remotely-sensed data and model data to source environmental covariates such as temperature, salinity, chlorophyll-a, mixed layer thickness, mesoscale features, etc. Download netCDF files for the specified
bounding box (40 km square) and time range for each site to the specified folder within "Environmental Data".

*Elephant Island (EI):* 03/05/2014 to 07/17/2014; latitude range: -60.52602,-61.24778; longitude range: -56.69255, -55.21545; site coordinates: (-60.8869, -55.95400)

*King George Island (KGI):* 02/10/2015 to 01/29/2016; latitude range:  -61.09694, -61.8187; longitude range: -58.69396, -57.18988; site coordinates: ( -61.457817, -57.941917)

*Clarence Island (CI):* 02/04/2016 to 12/02/2016; latitude range: -60.89099, -61.61275; longitude range: -54.23054, -52.73632, site coordinates: (-61.251867, -53.483433)

*Downloaded/accessed data:* [Aqua MODIS ](https://cwcgom.aoml.noaa.gov/erddap/griddap/miamiModisAquaChlor.html)(chlrophyll-a), SMOS (salinity), [ERDDAP](https://coastwatch.noaa.gov/cwn/products/sea-surface-salinity-near-real-time-miras-smos.html) (SST, salinity), AVISO (FSLEs; data downloaded on working disk), [HYCOM](https://tds.hycom.org/thredds/catalogs/GLBv0.08/expt_53.X.html) (salinity, temperature at depths), OCNET (chlrophyll-a; did not use because out of bounding box),
[University of Bremen](https://seaice.uni-bremen.de/sea-ice-concentration/amsre-amsr2/) (sea ice, code does not work), Copernicus (from [ocean physics](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/download) and [biogeochemistr](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description)y: salinity, temperature, chlorophyll, net primary production, EKE, ice variables, mixed layer thickness, SSH; at depth when applicable),
[daily AAO index](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/aao.shtml)

#### Next, use the script for each data source.
Scripts are titled GetSOURCE.R. Run this to read the specified variable/s from the netCDF files, make a stacked timeseries of each variable at each site, and save the daily bounding-box averaged values of each variable
to a .csv format.

### Step 3: Determine ACF values.
*Input:* Sheet with all species/site encounters (Antarc_Odontocetes.csv)
*Output:* ACF table (acf_table.csv), PACF table (pacf_table.csv), detections table (binned_detections.csv), daily presence/absence dataframe (dailyDetections.csv, in "Data" folder)
In the "Autocorrelation" folder, use ACF_SiteSpecies.R to generate a table of ACF values (daily) and save it to acf_table.csv. Also, create a similar table of days with detected presence (binned_detections.csv), to 
quickly quantify which sites and species have enough data to develop a model. Also creates a dataframe with all days HARPs were active and whether or not specific species were detected at each site, saved to 
DailyDetections.csv in the "Data" folder.

### Step 4: Create final dataframe to use in models.
*Input:* dailyDetections.csv ("Data" folder), copernicus_40km.csv (in "Copernicus" folder of "Environmental Data"), *SITECODE*_fsle_40km (in "AVISO" folder of "Environmental Data"), Daily_AAO.csv (in "Environmental Data" folder)
*Output:* allData_40km.csv ("Data" folder)
In the main folder, use FinalData.R to create a dataframe that combines daily binary presence for each species with all the environmental variables (including at biologically relevant depth bins). This dataframe
will be used to create the GAMs. User can also specify species, sites, and variables of interest to create a stacked timeseries that visualizes temporal trends in presence.

### Step 5: Make GAMs to model site-by-site, species-by-species presence on environmental predictors.
*Input:* allData_40km.csv ("Data" folder)
*Output:* GAMs for chosen species
In the "GAMs" folder, there is one modeling script for each species that creates site-specific models for presence. The code automatically bins the final dataframe by ACF value (instead of daily bins) and filters for the relevant species (and depth levels). Then, it plots a timeseries for that species presence across all sites using the binned data. After that, a VIF analysis is conducted on all potential environmental variables to remove highly correlated variables. 

In order to build the final GAMs, a single regression is first made with each of the remaining environmental predictors after the VIF analysis. Then, variables are added to the model one at a time until a final model is determined with only significant predictors. Components of the model such as number of knots, smoothing parameter, weights, etc. are also changed as needed to settle on the best model.


### Step 6: Visualize GAMs

## Notes
- **Cop_vs_HYCOM.R:** This script (located in the "Environmental Data" folder) was used to determine whether Copernicus or HYCOM provided a better model for remotely sensed environmental data that had too many gaps (such as SST). Copernicus was determined to be better and was used as a primary data source for the GAMs.
- **"Site Map" folder:** This script (in MATLAB) generates a map of the Antarctic sites from which HARP data was gathered.
