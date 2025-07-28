## Overview
Habitat modeling code for odontocete detections at three HARP sites off the Antarctic Peninsula (near the South Shetland Islands) from 2014 to 2016. Detected odontocetes include sperm whales, killer whales, 
long-finned pilot whales, southern bottlenose whales (BW29), and Gray's and strap-toothed whales (BW37 and BW58). Detected baleen whales will be habitat modeled separately.

## Procedure
### Step 1: Create timeseries of species detections.
Loaded detections (with start and end times) are in three files in the "Data" folder titled Antarc_SITECODE_01_Odontocetes.csv. Use TimeseriesPlots.R to create a site-by-site, species-by-species timeseries of detections.
Options include total hours of clicks detected per day and number of days with binary daily presence per week.

### Step 2: Load environmental data and create timeseries.
#### First, find data sources.
Look for satellite remotely-sensed data and model data to source environmental covariates such as temperature, salinity, chlorophyll-a, mixed layer thickness, mesoscale features, etc. Download netCDF files for the specified
bounding box (40 km square) and time range for each site to the specified folder within "Environmental Data".

*Elephant Island (EI):* 03/05/2014 to 07/17/2014; latitude range: -60.52602,-61.24778; longitude range: -56.69255, -55.21545; site coordinates: (-60.8869, -55.95400)

*King George Island (KGI):* 02/10/2015 to 01/29/2016; latitude range:  -61.09694, -61.8187; longitude range: -58.69396, -57.18988; site coordinates: ( -61.457817, -57.941917)

*Clarence Island (CI):* 02/04/2016 to 12/02/2016; latitude range: -60.89099, -61.61275; longitude range: -54.23054, -52.73632, site coordinates: (-61.251867, -53.483433)

*Downloaded/accessed data:* Aqua MODIS (chlrophyll-a), SMOS (salinity), ERDDAP (SST, salinity), AVISO (FSLEs; data downloaded on working disk), HYCOM (salinity, temperature at depths), OCNET (chlrophyll-a),
University of Bremen (sea ice, code does not work), Copernicus (salinity, temperature, chlorophyll, net primary production, EKE, ice variables, mixed layer thickness, SSH; at depth when applicable),
daily AAO index

#### Next, use the script for each data source.
Scripts are titled GetSOURCE.R. Run this to read the specified variable/s from the netCDF files, make a stacked timeseries of each variable at each site, and save the daily bounding-box averaged values of each variable
to a .csv format.

### Step 3: Determine ACF values.

### Step 4: Create final dataframe to use in models.

### Step 5: Make GAMs to model site-by-site, species-by-species presence on environmental predictors.

## Notes
Cop v hycom script to choose which model
