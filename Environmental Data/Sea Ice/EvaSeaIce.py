import numpy as np
import rasterio as rs
from rasterio.mask import mask
import requests
from io import BytesIO
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

# path to shape file (locate within same folder of notebook)
mask_path = 'C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Sea Ice/EI_circle.shp' #your circular or polygon mask shape file created in qgis
site_path = 'C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Sea Ice/EI_point.shp' #your coordinates as a point shape file created in qgis

# For file saving:
maskLocation = 'SeaIce' # Name/prefix that will be used for naming the files that will be created
path = 'C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Sea Ice/Data' # Where you want your files to be saved. 

# Note: The full saving location will be as follows {path}/{maskLocation}/{year}_{maskLocation}_MeanIceConcentration.png'. Within your path, you need to have created a folder named exactly as your "maskLocation". Otherwise it will throw an error.


# Load shape files for the sites and masks
site_gpd = gpd.read_file(site_path)
mask_gpd = gpd.read_file(mask_path)
mask_gpd = mask_gpd.to_crs('EPSG:3413')  # Match the projections to the target one

# Grab the important stuff for projection
canarc_bs = site_gpd.iloc[0]
bs_mask = mask_gpd.iloc[0].geometry

# years for processing:
years = ['2014'] #Edit as needed. Process one year at a time (tried doing more, but it gets messy and the trade off wasn't worth it).
process_months = [3,4,5,6,7]  # Edit as needed. Specific months can be removed if not needed.

# These lines below were added to add marks as "dots" into the resulting plots, reflecting the first/last days of narwhal detections for each year. 

# specific_dates = {
#     '2013': ('2013-09-14','2013-10-07'),
#     '2014': ('2014-07-25','2014-09-09','2014-10-04'),
#     '2015': ('2015-06-05','2015-08-30','2015-09-12','2015-10-12'),
#     '2016': ('2016-07-03','2016-07-13'),
#     '2017': ('2017-09-03','2017-10-01'),
#     '2018': ('2018-07-04','2018-09-05')
# }


# Assuming non-leap year; adjust February to 29 for leap years.

months_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
months_str = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']

# Initialize a dictionary to store mean values
mean_values_by_date = {}

def reproject_raster(src, target_crs):
    '''
    Change CRS of imported data.
    Input: Source data (src), which is a raster. Every file that we are downloading, is a .tif.
    Output: Reprojected data (previously raw, now is still a raster with values but in the projection we want), a transformed matrix.
    Additionally it also gives you the target CRS to keep track.
    '''

    src_crs = src.crs
    src_transform = src.transform
    src_width = src.width
    src_height = src.height

    # Reproject the raster data to the target CRS
    reprojected_data, dst_transform = rs.warp.reproject(
        source=rs.band(src, 1),
        src_transform=src_transform,
        src_crs=src_crs,
        dst_crs=target_crs,
        resampling=rs.enums.Resampling.nearest)
    
    return reprojected_data, dst_transform, target_crs

# Loop through each year, month and day:
for year in years:
    for month_index, days_in_month in enumerate(months_days):
        month_int = month_index + 1 
        if month_int in process_months:  # Will filter for the specified months.
            month_str = months_str[month_index]
            for day in range(1, days_in_month + 1):
                date = f'{year}{month_int:02d}{day:02d}'
                url = f'https://data.seaice.uni-bremen.de/amsr2/asi_daygrid_swath/s3125/{year}/{month_str}/Antarctic3125/asi-AMSR2-s3125-{date}-v5.4.tif'
                print(f'Accessing data for {date}: {url}')

                # Request and process data for each day
                response = requests.get(url)
                if response.status_code == 200:
                    today_seaice_raw = rs.open(BytesIO(response.content))
                    today_seaice_3995, transform_3995, crs_3995 = reproject_raster(today_seaice_raw, 'EPSG:3995')

                    # Apply the masks
                    masked_bs, transform_bs = mask(today_seaice_raw, [bs_mask], pad=False, crop=True, nodata=255)

                    # Calculate mean sea ice concentration, excluding land and mask values
                    bs_mean = np.mean(masked_bs[0][(masked_bs[0] != 255) & (masked_bs[0] !=120)])
                    mean_values_by_date[date] = {'bs_mean': bs_mean}
                else:
                    print(f'Failed to download data for {date} | Status code: {response.status_code}')
                    mean_values_by_date[date] = {'bs_mean': np.nan}

    # Convert the dictionary to a DataFrame and plot
    df = pd.DataFrame.from_dict(mean_values_by_date, orient='index', columns=['bs_mean'])
    df.index = pd.to_datetime(df.index)

    # Save csv:
csv_path = f'{path}/{maskLocation}/{year}_{maskLocation}_Arctic3125_3125Km.csv'
df.to_csv(csv_path, index=True, index_label='Date')
print(f'Data saved to {csv_path}')

# Loop through each year to plot and save ice concentrations:

for year in years:
    # Filter dataframe by year
    df_year = df[df.index.year == int(year)]

    # Plot and save ice concentrations:
    plt.figure(figsize=(10, 6))
    plt.plot(df_year.index, df_year['bs_mean'], label='Daily mean', marker='x')
    plt.title(f'Daily Sea Ice Concentration in {year} for {maskLocation}')
    plt.xlabel('Date')
    plt.ylabel('Mean Sea Ice %')
    plt.legend()
    plt.xticks(rotation=45)

    # # Add dots for specific dates reflecting first/last days of detections:
    # for specific_year, specific_dates_list in specific_dates.items():  
    #     for specific_date in specific_dates_list:
    #         if specific_date in df_year.index:
    #             plt.plot(df_year.index[df_year.index == specific_date], df_year['bs_mean'][df_year.index == specific_date], 'ro')
    
    plt.tight_layout()
    #plt.savefig(f'{path}/{maskLocation}/CANARC_BS_{year}_{maskLocation}_MeanIceConcentration.png')
    #plt.savefig(f'{path}/{maskLocation}/{year}_{maskLocation}_MeanIceConcentration.png')
    plt.close()

plt.show()
