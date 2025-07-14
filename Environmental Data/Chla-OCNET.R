library(ncdf4)
library(tidyverse)
library(raster)
library(gridExtra)

# ----------- Step 1: Extract Data ------------
EI_chla <- nc_open("FILEPATH")
KGI_chla <- nc_open("FILEPATH")
CI_chla <- nc_open("FILEPATH")

dfFromNC <- function(data) {
  ncvar_get
}
EI_df <- dfFromNC(EI_chla)
KGI_df <- dfFromNC(KGI_chla)
CI_df <- dfFromNC(CI_chla)

# ------------ Step 2: Make Timeseries -----------
# Compare to remotely sensed data

# ------------ Step 3: Save Dataframe ------------
EI_df$Site <- 'EI'
KGI_df$Site <- 'KGI'
CI_df$Site <- 'CI'

final <- rbind(EI_df,KGI_df,CI_df)
write.csv(final, "FILEPATH")