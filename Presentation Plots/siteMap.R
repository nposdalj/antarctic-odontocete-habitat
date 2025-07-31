library(tidyverse)
library(raster)
library(sp)
library(sf)
library(terra)
library(ggspatial)

bathymetry <- raster('/Users/trisha/R/antarctic-odontocete-habitat/Environmental Data/GMRT_bathymetry.grd')
# get projection
proj4string(bathymetry)
# change projection
chatham_crs <- CRS("+proj=tmerc +lat_0=-43.5 +lon_0=-176.5 +k=0.9996 
                    +x_0=1600000 +y_0=10000000 +datum=NZGD2000 +units=m +no_defs")
bathymetry_sp <- projectRaster(bathymetry, crs = chatham_crs)

ggplot() +
  ggspatial::geom_spatial(data = bathymetry) +
  coord_map()