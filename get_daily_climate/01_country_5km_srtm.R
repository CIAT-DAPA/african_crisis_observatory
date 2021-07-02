# Up-sampling of SRTM raster from 1 km to 5 km
# SRTM (Shuttle Radar Topography Mission): elevation model
# H. Achicanoy
# Alliance Bioversity-CIAT, 2021

# R options and load packages
options(warn = -1, scipen = 999)
suppressPackageStartupMessages(suppressMessages(pacman::p_load(tidyverse, raster, sf, RSAGA)))
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

get_5km_srtm <- function(iso = 'KEN', country = 'Kenya'){
  
  # CHIRPS data (ref raster 5 km)
  chr_tmp <- raster::raster('//catalogue/BaseLineDataCluster01/observed/gridded_products/chirps/daily/chirps-v2.0.1981.01.01.tif')
  chr_tmp[which(chr_tmp[] == -9999)] <- NA
  chr_tmp[which(chr_tmp[] < 0)] <- 0
  
  # Get country shapefile
  shp <- raster::getData(name = 'GADM', country = iso, level = 0, download = TRUE)
  shp <- lowest_gadm(iso = iso, out = NULL)
  
  # SRTM
  srtm <- raster::getData(name = 'alt', country = iso, mask = TRUE)
  srtm <- raster::crop(x = srtm, y = raster::extent(shp))
  srtm_5km <- raster::resample(srtm, chr_tmp, method = 'bilinear')
  srtm_5km <- srtm_5km %>%
    raster::crop(., extent(shp)) %>%
    raster::mask(., mask = shp)
  out <- paste0(root,'/cpc_data/srtm/',tolower(country),'/srtm_5km.tif')
  if(!dir.exists(dirname(out))){dir.create(dirname(out), FALSE, TRUE)}
  raster::writeRaster(srtm_5km, out, overwrite = T)
  
  return(cat('Done\n'))
  
}
