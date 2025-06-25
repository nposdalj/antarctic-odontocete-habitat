library(ncdf4)
library(lubridate)
library(CopernicusMarine)
library(stars)

CopernicusMarine_uid = "nposdaljian"
CopernicusMarine_pwd = "Phd0c3an0graphy!"

YEARS = c(2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020)

for (i in 1:length(YEARS)){
  year = YEARS[i]
  out_dir = paste("L:/My Drive/FourLoko/HabitatVariables/CCE_",as.character(year),".nc",sep = "")
  timerange1 = paste(as.character(year),"-01-01",sep="")
  timerange2 = paste(as.character(year),"-12-31",sep="")

copernicus_download_motu(
  username = CopernicusMarine_uid,
  password = CopernicusMarine_pwd,
  destination = out_dir,
  product = "GLOBAL_REANALYSIS_PHY_001_031",
  layer = "global-reanalysis-phy-001-031-grepv2-mnstd-daily",
  variable = "sea_surface_height",
  output = "netcdf",
  region = c(-130, 28, -113, 48),
  timerange = c(timerange1,timerange2),
  sub_variables = c("mlotst_mean","so_mean","thetao_mean","uo_mean",
  "vo_mean","zos_mean"),
  # sub_variables = c("mlotst_mean","mlotst_std","so_mean","so_std","thetao_mean","thetao_std","uo_mean",
  #                   "uo_std","vo_mean","vo_std","zos_mean","zos_std"),
  verticalrange = c(0.5057600140571594,1.56)
)
}