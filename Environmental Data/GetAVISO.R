library(rerddap)
library(rerddapXtracto)
library(ncdf4)
library(parsedate)
library(sp)
library(gganimate)
library(ggplot2)
library(plotdap)
library(abind)
library(dplyr)

# Only relevant AVISO variables are FSLE
# Ignore functions for other variables (data already obtained from Copernicus)

# --------------------Step 1: Subset and Save Relevant AVISO data--------------------------
subsetAVISO <- function(site) {
  directory <- "D:/FSLE/AVISO_Global_FSLE/dt_global_allsat_madt_fsle_"
  
  if (site == 'EI'){
    # site lat-long: -60.8869, -55.95400
    latitude = c(-59.3869,-61.3869)
    longitude = c(303.406, 304.406)
    years = c(2014)
    start <- as.Date('20140305',format='%Y%m%d')
    end <- as.Date('20140717',format='%Y%m%d') }
  if (site == 'KGI'){
    # site lat-long: -61.457817, -57.941917
    latitude = c(-60.957817,-61.957817)
    longitude = c(302.558083, 303.558083)
    years = c(2015,2016)
    start <- as.Date('20150210',format='%Y%m%d')
    end <- as.Date('20160129',format='%Y%m%d') }
  if (site == 'CI'){
    # site lat-long: -61.251867, -53.483433
    latitude = c(-60.751867,-61.751867)
    longitude = c(306.016567, 307.016567)
    years = c(2016)
    start <- as.Date('20160204',format='%Y%m%d')
    end <- as.Date('20161202',format='%Y%m%d') }
  
  dates <- seq(start, end, by = "day")
  df <- data.frame()
  for(current in dates) {
    date_str <- format(as.Date(current), "%Y%m%d")
    file_path <- paste(directory,date_str,'_20180704.nc',sep="")
    
    if (!file.exists(file_path)) {
      warning(paste("File does not exist:", file_path))
      next
    }
    nc <- nc_open(file_path)
    
    # Extract lat long and constrain to bounding box
    lon <- ncvar_get(nc, 'lon')
    lat <- ncvar_get(nc,'lat')
    time <- ncvar_get(nc, 'time')
    # Define bounding box
    lat_min <- latitude[2]
    lat_max <- latitude[1]
    lon_min <- longitude[1]
    lon_max <- longitude[2]
    # Find index ranges for the bounding box
    lat_idx <- which(lat >= lat_min & lat <= lat_max)
    lon_idx <- which(lon >= lon_min & lon <= lon_max)
    
    fsle <- ncvar_get(nc, "fsle_max", start = c(min(lon_idx), min(lat_idx),1),
                     count = c(length(lon_idx), length(lat_idx),1))
    
    grid <- expand.grid(lon = lon[lon_idx], lat = lat[lat_idx], time = time)
    grid$temp <- as.vector(fsle)
    grid$date <- as.Date(current)
    
    # Average over space for each depth
    daily_avg <- grid %>%
      group_by(date) %>%
      summarise(fsle = mean(fsle, na.rm = TRUE),
                .groups = "drop")
    
    df <- bind_rows(df, daily_avg)
    print(paste("File for ", current, "done."))
  }
  return(df)
}
EI_fsle <- subsetAVISO('EI')
write.csv(EI_fsle, "D:/FSLE/EI_fsle")
KGI_fsle <- subsetAVISO('KGI')
write.csv(EI_fsle, "D:/FSLE/KGI_fsle")
CI_fsle <- subsetAVISO('CI')
write.csv(EI_fsle, "D:/FSLE/CI_fsle")

LoadAVISO <- function(envDir){
  
  filenameStatAll = paste(envDir,"dt_global_allsat_madt_fsle_20140101_20180704",sep="")#load files as data frame
  names(AVISO$var)
}
envDir <- "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/AVISO/FSLE/"

readFSLE <- function(site, start, end) {
  start <- as.Date(start)
  end <- as.Date(end)
  dates <- seq(start, end, by = "day")
  
  df <- data.frame()
  
  for (current in dates) {
    # Format date into yyyymmdd
    date_str <- format(as.Date(current), "%Y%m%d")
    
    # Construct filename (you can make site-specific filename logic here if needed)
    file_path <- paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/AVISO/FSLE/",
                        year, "/dt_global_allsat_madt_fsle_", date_str,"_20210921.nc.filepart")
    
    if (!file.exists(file_path)) {
      warning(paste("File does not exist:", file_path))
      next
    }
    
    # Read the .mat file
    mat <- readMat(file_path)
    
    # Extract relevant variables (adapt these to match the actual structure of your .mat files)
    lon <- as.vector(mat$lon)
    lat <- as.vector(mat$lat)
    depth <- as.vector(mat$depth)
    temp <- mat$temp
    salt <- mat$salt
    
    # Assume temp and salt are 3D arrays: lon x lat x depth
    # Flatten arrays to data frame
    nlon <- length(lon)
    nlat <- length(lat)
    ndepth <- length(depth)
    
    grid <- expand.grid(lon = lon, lat = lat, depth = depth)
    grid$temp <- as.vector(temp)
    grid$salt <- as.vector(salt)
    grid$date <- as.Date(current)
    
    # Average over space for each depth
    daily_avg <- grid %>%
      group_by(depth, date) %>%
      summarise(temp = mean(temp, na.rm = TRUE),
                salt = mean(salt, na.rm = TRUE),
                .groups = "drop")
    
    df <- bind_rows(df, daily_avg)
  }
  #df$date <- as.POSIXct(df$date, origin = "1970-01-01", tz="UTC") # Changing from seconds to time
  
  return(df)
  
}
EI_df <- readHYCOM('EI','2014-03-05','2014-07-17')
KGI_df <- readHYCOM('KGI','2015-02-10','2016-01-29') 
CI_df <- readHYCOM('CI','2016-02-04','2016-12-02')

# # ---------------- Step NA: Functions for other variables--------------
# GetSSH <- function(envDir,Region){
#   inFilePaths = list.files(path=envDir, pattern=Region, full.names=TRUE) #look for all files
#   for (ss in 1:(length(inFilePaths)-1)){
#   AVISO = nc_open(inFilePaths[ss])
#   #zos - SSH
#   v6=AVISO$var[[3]]
#   SSHvar=ncvar_get(AVISO,v6)
#   SSH_lon=v6$dim[[1]]$vals
#   SSH_lat=v6$dim[[2]]$vals
#   SSH_dates=as.POSIXlt(v6$dim[[3]]$vals*60*60*24,origin='1950-01-01') #extract the date/time
#   SSH_dates <<- as.Date(SSH_dates, format = "%m/%d/%y") #get rid of the time
#   
#   #plotting timeseries
#   I=which(SSH_lon>=min(df1$long) & SSH_lon<= max(df1$long)) #only extract the region we care about
#   if (length(I) == 1){ #if the longitude only has 1 value, add a second
#     II = I:(I+1)
#   }else{
#     II = I
#   }
#   J=which(SSH_lat>=min(df1$lat) & SSH_lat<=max(df1$lat)) #only extract the region we care about
#   if (length(J) == 1){ #if the latitude only has 1 value, add a second
#     JJ = J:(J+1)
#   }else{
#     JJ = J
#   }
#     Ka = which(SSH_dates<startTime)
#     Kb = which(SSH_dates>endTime)
#   if (length(Ka) == 0 & length(Kb) == 0){
#     SSH2=SSHvar[II,JJ,] #index the original data frame to extract the lat, long, dates we care about
#   }else{
#     SSH2=SSHvar[II,JJ,-c(Ka,Kb)] #index the original data frame to extract the lat, long, dates we care about
#     SSH_dates = SSH_dates[-c(Ka,Kb)]
#   }
#   
#   n=dim(SSH2)[3] #find the length of time
#   
#   #take the mean
#   resSSH=rep(NA,n) 
#   for (i in 1:n){ 
#     resSSH[i]=mean(SSH2[,,i],na.rm=TRUE)}
#   
#     SSH_ddf <<- as.data.frame(SSH_dates)
#     SSH_ddf <- plyr::rename(SSH_ddf, replace = c("SSH_dates" = "time"))
# 
#   SSH_ddf$time=as.Date(SSH_ddf$time)
#   SSHdf<<- bind_cols(SSH_ddf,as.data.frame(resSSH))
#   assign(paste0("SSHdf_", ss), SSHdf)
#   }
#   
#   if(Region == 'CCE' | Region == 'CentralPac'){
#     SSH_list = list(SSHdf_1, SSHdf_2,SSHdf_3,SSHdf_4,SSHdf_5,SSHdf_6,SSHdf_7,
#                     SSHdf_8,SSHdf_9,SSHdf_10,SSHdf_11,SSHdf_12,SSHdf_13,SSHdf_14,SSHdf_15)
#   }else{
#     SSH_list = list(SSHdf_1, SSHdf_2,SSHdf_3,SSHdf_4,SSHdf_5,SSHdf_6,SSHdf_7,
#                       SSHdf_8,SSHdf_9,SSHdf_10)
#   }
#   SSHdf = bind_rows(SSH_list)
#   SSHdf[SSH < 0] <- NA
#   
#   #plot the time series
#   plot(SSHdf$time,SSHdf$resSSH,type='o',pch=20,xlab='Year',ylab='SSH',las = 3) 
# }
# 
# GetDEN <- function(AVISO){
# #mlotst - density ocean mixed layer thickness
# v1=AVISO$var[[2]]
# DENvar=ncvar_get(AVISO,v1)
# DEN_lon=v1$dim[[1]]$vals
# DEN_lat=v1$dim[[2]]$vals
# DEN_dates=as.POSIXlt(v1$dim[[3]]$vals*60*60*24,origin='1950-01-01') #extract the date/time
# DEN_dates = as.Date(DEN_dates, format = "%m/%d/%y") #get rid of the time
# 
# #plotting time series SAPTIN 
# I=which(DEN_lon>=min(df1$long) & DEN_lon<= max(df1$long)) #only extract the region we care about
# if (length(I) == 1){ #if the longitude only has 1 value, add a second
#   II = I:(I+1)
# }else{
#   II = I
# }
# J=which(DEN_lat>=min(df1$lat) & DEN_lat<=max(df1$lat)) #only extract the region we care about
# if (length(J) == 1){ #if the latitude only has 1 value, add a second
#   JJ = J:(J+1)
# }else{
#   JJ = J
# }
# K=which(DEN_dates>= startTime & DEN_dates<= endTime) #extract only the dates we care about
# DEN2=DENvar[II,JJ,K] #index the original data frame to extract the lat, long, dates we care about
# 
# n=dim(DEN2)[3] #find the length of time
# 
# #take the mean
# resDEN=rep(NA,n) 
# for (i in 1:n) 
#   resDEN[i]=mean(DEN2[,,i],na.rm=TRUE) 
# 
# DEN_ddf <- as.data.frame(DEN_dates[K])
# DEN_ddf = DEN_ddf %>% 
#   rename(
#     time = 'DEN_dates[K]',
#   )
# 
# DENdf<<- bind_cols(DEN_ddf,as.data.frame(resDEN))
# 
# #plot the time series
# plot(1:n,resDEN,axes=FALSE,type='o',pch=20,xlab='',ylab='Density',las = 3) 
# axis(2) 
# axis(1,1:n,format(DEN_dates[K]),las = 3) 
# box()
# }
# 
# GetSAL <- function(AVISO){
# #so - salinity
# v2=AVISO$var[[5]]
# SALvar=ncvar_get(AVISO,v2)
# SAL_lon=v2$dim[[1]]$vals
# SAL_lat=v2$dim[[2]]$vals
# SAL_dates=as.POSIXlt(v2$dim[[4]]$vals*60*60*24,origin='1950-01-01') #extract the date/time
# SAL_dates = as.Date(SAL_dates, format = "%m/%d/%y") #get rid of the time
# 
# #plotting time series SAPTIN 
# I=which(SAL_lon>=min(df1$long) & SAL_lon<= max(df1$long)) #only extract the region we care about
# if (length(I) == 1){ #if the longitude only has 1 value, add a second
#   II = I:(I+1)
# }else{
#   II = I
# }
# J=which(SAL_lat>=min(df1$lat) & SAL_lat<=max(df1$lat)) #only extract the region we care about
# if (length(J) == 1){ #if the latitude only has 1 value, add a second
#   JJ = J:(J+1)
# }else{
#   JJ = J
# }
# K=which(SAL_dates>= startTime & SAL_dates<= endTime) #extract only the dates we care about
# SAL2=SALvar[II,JJ,,K] #index the original data frame to extract the lat, long, dates we care about
# 
# n=dim(SAL2)[4] #find the length of time
# 
# #take the mean
# resSAL=rep(NA,n) 
# for (i in 1:n) 
#   resSAL[i]=mean(SAL2[,,,i],na.rm=TRUE)
# 
# SAL_ddf <- as.data.frame(SAL_dates[K])
# SAL_ddf = SAL_ddf %>% 
#   rename(
#     time = 'SAL_dates[K]',
#   )
# 
# SALdf<<- bind_cols(SAL_ddf,as.data.frame(resSAL))
# 
# #plot the time series
# plot(1:n,resSAL,axes=FALSE,type='o',pch=20,xlab='',ylab='Salinity',las = 3) 
# axis(2) 
# axis(1,1:n,format(SAL_dates[K]),las = 3) 
# box()
# }
# 
# # ----------------------Step 2: SST-------------------
# GetTEMP <- function(AVISO){
# #thetao - temperature
# v3=AVISO$var[[1]]
# TEMPvar=ncvar_get(AVISO,v3)
# TEMP_lon=v3$dim[[1]]$vals
# TEMP_lat=v3$dim[[2]]$vals
# TEMP_dates=as.POSIXlt(v3$dim[[4]]$vals*60*60*24,origin='1950-01-01') #extract the date/time
# TEMP_dates = as.Date(TEMP_dates, format = "%m/%d/%y") #get rid of the time
# 
# #plotting time series SAPTIN 
# I=which(TEMP_lon>=min(df1$long) & TEMP_lon<= max(df1$long)) #only extract the region we care about
# if (length(I) == 1){ #if the longitude only has 1 value, add a second
#   II = I:(I+1)
# }else{
#   II = I
# }
# J=which(TEMP_lat>=min(df1$lat) & TEMP_lat<=max(df1$lat)) #only extract the region we care about
# if (length(J) == 1){ #if the latitude only has 1 value, add a second
#   JJ = J:(J+1)
# }else{
#   JJ = J
# }
# K=which(TEMP_dates>= startTime & TEMP_dates<= endTime) #extract only the dates we care about
# TEMP2=TEMPvar[II,JJ,,K] #index the original data frame to extract the lat, long, dates we care about
# 
# n=dim(TEMP2)[4] #find the length of time
# 
# #take the mean
# resTEMP=rep(NA,n) 
# for (i in 1:n) 
#   resTEMP[i]=mean(TEMP2[,,,i],na.rm=TRUE)
# 
# TEMP_ddf <- as.data.frame(TEMP_dates[K])
# TEMP_ddf = TEMP_ddf %>% 
#   rename(
#     time = 'TEMP_dates[K]',
#   )
# 
# TEMPdf<<- bind_cols(TEMP_ddf,as.data.frame(resTEMP))
# 
# #plot the time series
# plot(1:n,resTEMP,axes=FALSE,type='o',pch=20,xlab='',ylab='Temperature',las = 3) 
# axis(2) 
# axis(1,1:n,format(TEMP_dates[K]),las = 3) 
# box()
# }

# ------------------- Step 3: EKE ---------------------
# Function for EKE from velocity:
# GetEKE <- function(AVISO){
# #uo - eastward velocity
# v4=AVISO$var[[4]]
# EASTVvar=ncvar_get(AVISO,v4)
# EASTV_lon=v4$dim[[1]]$vals
# EASTV_lat=v4$dim[[2]]$vals
# EAST_dates=as.POSIXlt(v4$dim[[4]]$vals*60*60*24,origin='1950-01-01') #extract the date/time
# EAST_dates = as.Date(EAST_dates, format = "%y/%m/%d") #get rid of the time
# EASTdf <- as.data.frame(EASTVvar)
# 
# #plotting time series SAPTIN 
# I=which(EASTV_lon>=min(df1$long) & EASTV_lon<= max(df1$long)) #only extract the region we care about
# if (length(I) == 1){ #if the longitude only has 1 value, add a second
#   II = I:(I+1)
# }else{
#   II = I
# }
# J=which(EASTV_lat>=min(df1$lat) & EASTV_lat<=max(df1$lat)) #only extract the region we care about
# if (length(J) == 1){ #if the latitude only has 1 value, add a second
#   JJ = J:(J+1)
# }else{
#   JJ = J
# }
# K=which(EAST_dates>= startTime & EAST_dates<= endTime) #extract only the dates we care about
# EASTV2=EASTVvar[II,JJ,,K] #index the original data frame to extract the lat, long, dates we care about
# 
# n=dim(EASTV2)[4] #find the length of time
# 
# #take the mean
# resEV=rep(NA,n) 
# for (i in 1:n) 
#   resEV[i]=mean(EASTV2[,,,i],na.rm=TRUE) 
# 
# EV_ddf <- as.data.frame(EAST_dates[K])
# EV_ddf = EV_ddf %>% 
#   rename(
#     time = 'EAST_dates[K]',
#   )
# 
# EVdf<- bind_cols(EV_ddf,as.data.frame(resEV))
# 
# #plot the time series
# plot(1:n,resEV,axes=FALSE,type='o',pch=20,xlab='',ylab='Eastward Velocity',las = 3) 
# axis(2) 
# axis(1,1:n,format(EAST_dates[K]),las = 3) 
# box()
# 
# #vo - northward velocity
# v5=AVISO$var[[6]]
# NORVvar=ncvar_get(AVISO,v5)
# NORV_lon=v5$dim[[1]]$vals
# NORV_lat=v5$dim[[2]]$vals
# NOR_dates=as.POSIXlt(v5$dim[[4]]$vals*60*60*24,origin='1950-01-01') #extract the date/time
# NOR_dates = as.Date(NOR_dates, format = "%m/%d/%y") #get rid of the time
# NORdf <- as.data.frame(NORVvar)
# 
# #plotting time series SAPTIN 
# I=which(NORV_lon>=min(df1$long) & NORV_lon<= max(df1$long)) #only extract the region we care about
# if (length(I) == 1){ #if the longitude only has 1 value, add a second
#   II = I:(I+1)
# }else{
#   II = I
# }
# J=which(NORV_lat>=min(df1$lat) & NORV_lat<=max(df1$lat)) #only extract the region we care about
# if (length(J) == 1){ #if the latitude only has 1 value, add a second
#   JJ = J:(J+1)
# }else{
#   JJ = J
# }
# K=which(NOR_dates>= startTime & NOR_dates<= endTime) #extract only the dates we care about
# NORV2=NORVvar[II,JJ,,K] #index the original data frame to extract the lat, long, dates we care about
# 
# n=dim(NORV2)[4] #find the length of time
# 
# #take the mean
# resNV=rep(NA,n) 
# for (i in 1:n) 
#   resNV[i]=mean(NORV2[,,,i],na.rm=TRUE)
# 
# NV_ddf <- as.data.frame(NOR_dates[K])
# NV_ddf = NV_ddf %>% 
#   rename(
#     time = 'NOR_dates[K]',
#   )
# 
# NVdf<- bind_cols(NV_ddf,as.data.frame(resNV))
# 
# #plot the time series
# plot(1:n,resNV,axes=FALSE,type='o',pch=20,xlab='',ylab='Northward Velocity',las = 3) 
# axis(2) 
# axis(1,1:n,format(NOR_dates[K]),las = 3) 
# box()
# 
# ####
# 
# #Calculate EKE
# u <- (resEV)^2
# v <- (resNV)^2
# velocity = u + v
# EKE_meters = 0.5 * velocity
# EKE_cm = EKE_meters * 10000 
# 
# EKE <- bind_cols(SSH_ddf,as.data.frame(EKE_cm))
# colnames(EKE)[1] <- "time"
# EKE$time=as.Date(EKE$time)
# EKE <<- EKE
# }
# 
# # Function to combine EKE files into one timeseries & plot
# accessEKE <- function() {
#   
# }
