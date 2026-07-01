# install.packages(c("robis", "dplyr", "sf", "readr"))

library(robis)
library(dplyr)
library(sf)
library(readr)

# --------------------------------------------------
# Western Antarctic Peninsula bounding box
# --------------------------------------------------

lon_min <- -71
lon_max <- -48
lat_min <- -69
lat_max <- -59

# --------------------------------------------------
# Create WKT polygon from bounding box
# --------------------------------------------------
# WKT = Well-Known Text, a standard way to describe spatial shapes.
# The polygon must close by repeating the first coordinate at the end.

wap_bbox_wkt <- paste0(
  "POLYGON ((",
  lon_min, " ", lat_min, ", ",
  lon_max, " ", lat_min, ", ",
  lon_max, " ", lat_max, ", ",
  lon_min, " ", lat_max, ", ",
  lon_min, " ", lat_min,
  "))"
)

# --------------------------------------------------
# Download fin whale records from OBIS
# --------------------------------------------------
# Scientific name for fin whale = Balaenoptera physalus

fin_whale_obis <- occurrence(
  scientificname = "Balaenoptera physalus",
  geometry = wap_bbox_wkt,
  absence = FALSE
)

# --------------------------------------------------
# Clean and keep useful columns
# --------------------------------------------------

fin_whale_clean <- fin_whale_obis %>%
  filter(
    !is.na(decimalLongitude),
    !is.na(decimalLatitude)
  ) %>%
  select(
    scientificName,
    decimalLongitude,
    decimalLatitude,
    eventDate,
    year,
    basisOfRecord,
    datasetName,
    institutionCode,
    catalogNumber,
    occurrenceID,
    individualCount,
    everything()
  )

# --------------------------------------------------
# Convert to spatial object
# --------------------------------------------------

fin_whale_sf <- st_as_sf(
  fin_whale_clean,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326,
  remove = FALSE
)

# --------------------------------------------------
# Save downloaded data
# --------------------------------------------------

output_dir <- "C:/Users/YourName/Documents/OBIS"

write_csv(
  fin_whale_clean,
  file.path(output_dir, "OBIS_fin_whale_WAP_bbox.csv")
)

st_write(
  fin_whale_sf,
  file.path(output_dir, "OBIS_fin_whale_WAP_bbox.gpkg"),
  delete_dsn = TRUE
)

# --------------------------------------------------
# Quick check
# --------------------------------------------------

nrow(fin_whale_clean)

head(fin_whale_clean)