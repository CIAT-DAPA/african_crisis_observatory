##Download MODIS NDVI and EVI https://modis.gsfc.nasa.gov/data/dataprod/mod13.php
#' ############################################################################################
#' Author: Benson Kenduiywo
#' ############################################################################################

#https://viirsland.gsfc.nasa.gov/Products/NASA/LAI_FparESDR.html
#https://lpdaac.usgs.gov/products/myd15a2hv006/
rm(list=ls(all=TRUE))
g <- gc(reset = T); 
root <-  '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/raw/modis/'
dir.create(root, recursive = T)
remotes::install_github("rspatial/luna")
list.of.packages <- c("luna","terra", "meteor")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)

prod <- getProducts("^MOD13Q1|^MYD13Q1|")
product <- 'MYD13Q1'
#productInfo(product)
#login information

#saveRDS(df, paste0(root,'earthdata.rds'))
df <- readRDS(paste0(root,'earthdata.rds'))

iso  <- 'KEN'
start <- "2008-11-08" #2008-11-08
end <- "2008-12-31"
aoi <- vect('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/raw/admin/KEN/Kenya_county_dd.shp')
path <- paste0(root, iso, '/ndvi/')
dir.create(path, recursive = T)
mf <- luna::getNASA(product, start, end, aoi=aoi, download=TRUE,
                     path=path, username=df$user, password=df$pwd)


tt <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/raw/modis/KEN/ndvi/MYD13Q1.A2008313.h22v08.061.2021111213053.hdf'



library(luna)
library(meteor)
selectModisFiles <- function(files, startdate, enddate) {
  base_names <- basename(files)
  dates <- luna::modisDate(base_names)
  i <- (dates >= as.Date(startdate)) & (dates <= as.Date(enddate) )
  files[i]
}

rlist <- purrr::map(mf, rast)
rlist <- terra::sprc(rlist)
im <- merge(rlist)
plot(im['Fpar_500m'])









#XXXXXXX
#Check https://flograttarola.com/post/modis-downloads/
library(remotes)
install_github("ropensci/MODIStsp")

library(MODIStsp)
library(terra)
library(sf)
library(tidyverse)

