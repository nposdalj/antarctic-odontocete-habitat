## Overview
Habitat modeling code for odontocete detections at three HARP sites off the Antarctic Peninsula (near the South Shetland Islands) from 2014 to 2016. Detected odontocetes include sperm whales, killer whales, 
long-finned pilot whales, southern bottlenose whales (BW29), and Gray's and strap-toothed whales (BW37 and BW58). Detected baleen whales will be habitat modeled separately. Models used are generalized additive models (GAMs). Use [this sheet](https://docs.google.com/spreadsheets/d/1OHM6oqbrj5t5F1CpTyd1P_VlXrbCY_-TRTlotMVn46c/edit?gid=887723042#gid=887723042) to keep track of the progress made on each model, methodologies tried, modeling results, and the environmental variables used.

#### DISCLAIMER: Repository was recently reorganized, so make sure to double-check the filepath of all data loaded into a script.

## Procedure
### Step 1: Create timeseries of species detections.
*Input:* Species detections data (start/end time format for each encounter)

*Output:* Timeseries of species detections

*Scripts Used:* TimeseriesPlots.R

Loaded detections (with start and end times) are in three files in the "Data" folder titled Antarc_*SITECODE*_01_Odontocetes.csv. Use TimeseriesPlots.R to create a site-by-site, species-by-species timeseries of detections.
Options include total hours of clicks detected per day and number of days with binary daily presence per week. Save all site detections to Antarc_Odontocetes.csv

### Step 2: Load environmental data and create timeseries.
*Input:* Environmental data (various sources, netCDF format, in various subfolders of "Environmental Data" or on the working disk for AVISO)

*Output:* Environmental data timeseries & .csv sheets for each site (timeseries stored in script, .csv sheets in various subfolders of "Environmental Data")

*Scripts Used:* GetAVISO.R, GetCopernicus.R, GetHYCOM.R, GetERDDAP_SST.R, GetERDDAP_SSS.R, GetERDDAP_Chla.R, TimeLaggedCopernicus.R

#### First, find data sources.
Look for satellite remotely-sensed data and model data to source environmental covariates such as temperature, salinity, chlorophyll-a, mixed layer thickness, mesoscale features, etc. Download netCDF files for the specified
bounding box (40 km square) and time range for each site to the specified folder within "Environmental Data".

*Elephant Island (EI):* 03/05/2014 to 07/17/2014; latitude range: -60.52602,-61.24778; longitude range: -56.69255, -55.21545; site coordinates: (-60.8869, -55.95400)

*King George Island (KGI):* 02/10/2015 to 01/29/2016; latitude range:  -61.09694, -61.8187; longitude range: -58.69396, -57.18988; site coordinates: ( -61.457817, -57.941917)

*Clarence Island (CI):* 02/04/2016 to 12/02/2016; latitude range: -60.89099, -61.61275; longitude range: -54.23054, -52.73632, site coordinates: (-61.251867, -53.483433)

*Downloaded/accessed data:* [Aqua MODIS ](https://cwcgom.aoml.noaa.gov/erddap/griddap/miamiModisAquaChlor.html)(chlrophyll-a), SMOS (salinity), [ERDDAP](https://coastwatch.noaa.gov/cwn/products/sea-surface-salinity-near-real-time-miras-smos.html) (SST, salinity), AVISO (FSLEs; data downloaded on working disk), [HYCOM](https://tds.hycom.org/thredds/catalogs/GLBv0.08/expt_53.X.html) (salinity, temperature at depths), OCNET (chlrophyll-a; did not use because out of bounding box),
[University of Bremen](https://seaice.uni-bremen.de/sea-ice-concentration/amsre-amsr2/) (sea ice, code does not work), Copernicus (from [ocean physics](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/download) and [biogeochemistry](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description): salinity, temperature, chlorophyll, net primary production, EKE, ice variables, mixed layer thickness, SSH; at depth when applicable),
[daily AAO index](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/aao.shtml)

#### Next, use the script for each data source.
Scripts are titled GetSOURCE.R. Run this to read the specified variable/s from the netCDF files, make a stacked timeseries of each variable at each site, and save the daily bounding-box averaged values of each variable to a .csv format.

#### Finally, determine the data sources to use for modeling.
Only Copernicus, AVISO, and AAO data was used for creating the GAMs. Remotely sensed data (Aqua MODIS, SMOS, ERDDAP sources) were found to have too many gaps in either the 40 km bounding box or daily timeseries to be used in the modeling code, but could be useful for additional ground-truthing Copernicus models. HYCOM was not used because Copernicus performed better and had more variables (analysis found in Cop_vs_HYCOM.R and [here](https://docs.google.com/document/d/11GEedOBmrZNMdO0CqTZY6Ltp4wF0kYFU6sATxuj63Ag/edit?tab=t.0)). OCNET did not have data for our sites. University of Bremen sea ice data is probably preferable to Copernicus, however the code to access this is not currently working, so Copernicus was used for convenience (but Bremen code and QGIS bounding box files can be found in the "Sea Ice Bremen" folder witihin "Environmental Data".

### Step 3: Determine ACF values.
*Input:* Sheet with all species/site encounters (Antarc_Odontocetes.csv)

*Output:* ACF table (acf_table.csv), PACF table (pacf_table.csv), detections table (binned_detections.csv), daily presence/absence dataframe (dailyDetections.csv, in "Data" folder)

*Scripts Used:* ACF_SiteSpecies.R

In the "Autocorrelation" folder, use ACF_SiteSpecies.R to generate a table of ACF values (daily) and save it to acf_table.csv. Also, create a similar table of days with detected presence (binned_detections.csv), to 
quickly quantify which sites and species have enough data to develop a model. Also creates a dataframe with all days HARPs were active and whether or not specific species were detected at each site, saved to 
DailyDetections.csv in the "Data" folder.

### Step 4: Create final dataframe to use in models.
*Input:* dailyDetections.csv ("Data" folder), copernicus_40km.csv (in "Copernicus" folder of "Environmental Data"), *SITECODE*_fsle_40km (in "AVISO" folder of "Environmental Data"), Daily_AAO.csv (in "Environmental Data" folder), copernicus_lagged.csv (in "Copernicus" folder of "Environmental Data"),

*Output:* allData_40km.csv ("Data" folder)

*Scripts Used:* FinalData.R

In the main folder, use FinalData.R to create a dataframe that combines daily binary presence for each species with all the environmental variables (including at biologically relevant depth bins and 1 to 6 month lags). This dataframe
will be used to create the GAMs. User can also specify species, sites, and variables of interest to create a stacked timeseries that visualizes temporal trends in presence.

### Step 5: Make GAMs to model odontocete presence based on environmental predictors.
*Input:* allData_40km.csv ("Data" folder)

*Output:* GAMs for chosen species

*Scripts Used:* Gm_surfaceGAM.R, Pm_surfaceGAM.R, BW29_surfaceGAM.R, BW37_surfaceGAM.R, Oo_surfaceGAM.R, allSite_surfaceGAM.R, Gm_depthsGAM.R, Pm_depthsGAM.R, BW29_depthsGAM.R, BW37_depthsGAM.R, Oo_depthsGAM.R (all in "GAMs" folder)

In the "GAMs" folder, there are modeling scripts for each species that create site-specific models for presence (as well as all site models). The code automatically bins the final dataframe by ACF value (instead of daily bins) and filters for the relevant species (and depth levels). Then, it plots a timeseries for that species presence across all sites using the binned data. After that, a VIF analysis is conducted on all potential environmental variables to remove highly correlated variables. 

In order to build the final GAMs, a single-variable GAM is first made with each of the remaining environmental predictors after the VIF analysis. Then, variables are added to the model one at a time until a final model is determined with only significant predictors. Components of the model such as number of knots and the smoothing parameters are also changed as needed to settle on the best model. Initial GAMs have been made for surface variables at for each site and a combined all site model. Each script also has code for visualizing the created models, though this may need to be modified to account for the three-dimensionality of the FSLE magnitude and orientation interaction term. Depth predictors have yet to be added in most cases, and each model varies slightly in terms of methodology and variables used. The status of all models and the variables they use can be found at: [Model and Environmental Variable Status Sheet](https://docs.google.com/spreadsheets/d/1OHM6oqbrj5t5F1CpTyd1P_VlXrbCY_-TRTlotMVn46c/edit?gid=887723042#gid=887723042)

*Note:* Framework for the depthsGAM scripts is written out (data loaded, timeseries created, VIF analysis, etc.), but the models themselves have not been created yet.


## Additional Visualizations
Located in "Data Visualizations" folder.
- **FSLE_GIF.R:** Makes a GIF of monthly FSLE across CI and KGI (Feb 2015 to Dec 2016) and GIF of FSLEs on a daily resolution for a chosen month. Images for the GIFs are created and stored in "FSLE" within "GIF Images". However, new dates can be chosen for GIFs (in the script) and new images generated for visualizations over any range (though this may take time to run).
- **SeaIceGIF.R:** Makes a GIF of monthly sea ice across all sites (with gaps for months without HARP data). Images were manually downloaded from [Universität Bremen](https://data.seaice.uni-bremen.de/databrowser/#p=sic).
- **siteMap.R:** Uses GMRT bathymetry data (stored in "GMRT" folder of "Environmental Data") to create a map of the South Shetland Islands region with all sites indicated.
- **DepthProfiles_Copernicus.R:** Uses Copernicus data to generate timeseries of salinity, temperature, oxygen concentration, and density (both value and deviation from mean) across depths. The depth profiles for density were derived using temperature, salinity, and a very rough estimate for pressure across depths (P = ρgh where ρ = 1023.6 kg/m3, g = 9.80665 m/s2, h = depth) that did not take into account the actual densities of the seawater above a depth bin.

## Notes
- **Cop_vs_HYCOM.R:** This script (located in the "Environmental Data" folder) was used to determine whether Copernicus or HYCOM provided a better model for remotely sensed environmental data that had too many gaps (such as SST). Copernicus was determined to be better and was used as a primary data source for the GAMs (analysis found in script and [here](https://docs.google.com/document/d/11GEedOBmrZNMdO0CqTZY6Ltp4wF0kYFU6sATxuj63Ag/edit?tab=t.0)).
- **OLDGm_GAM.R:** This script (located in the "GAMs" folder) was the first model constructed, trying to use data across all relevant depth bins. But, issues with the model output's confidence intervals led to this initial iteration being scrapped (and the surface only models being used as first steps for GAMs).
