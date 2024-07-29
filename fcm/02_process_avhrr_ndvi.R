#' ############################################################################################
#' Author: Benson Kenduiywo
#' ############################################################################################

rm(list=ls(all=TRUE))
g <- gc(reset = T); 
list.of.packages <- c("luna","terra","geodata")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)
path <-  '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/raw/avhrr/'
afiles <- list.files(path=path, pattern = "AVHRR-Land_*.*.nc", full.names = T)
vfiles <- list.files(path=path, pattern = "VIIRS-Land_*.*.nc", full.names = T)

processAVHRR <- function(yr, files){
  names <- gsub('_c.*', "", files)
  dates <- substr(names, nchar(names)-7, nchar(names))
  dd <- as.Date(dates, format = "%Y%m%d")
  years <- as.numeric(format(as.Date(dd, format="%Y-%m-%d"),"%Y"))
  df <- data.frame(files=files, dates=dates, years=years)
  df <- df[df$years==yr,]
  ff <- df$files
  ndvi <- terra::rast(ff, subds  = "NDVI")
  qa <- terra::rast(ff, subds = "QA")
 }

ndvi <- terra::rast(afiles[1:10], subds = "NDVI")
#How to unpack bits for masking out areas
