# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain agricultural area raster at 1 km resolution
# Original data source: http://www.earthstat.org/cropland-pasture-area-2000/
# Steps:
# 1. Download manually the rasters of cropland and pasture area from the data source
# 2. Execute this script to obtain a single rasters file with the total cropland and
#    pasture area cropped by the country of interest
# Author: Harold Achicanoy
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, terra, vegan, FactoMineR))

# Define root directory
root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'

crp <- terra::rast(paste0(root,'/data/_global/cropland/Cropland2000_5m.tif')) # Read global cropland raster
pas <- terra::rast(paste0(root,'/data/_global/cropland/Pasture2000_5m.tif'))  # Read global pasture area raster

smm <- crp + pas # Compute cropland and pasture area raster

isos <- 'PHL' #c('KEN','MLI','NGA','SDN','SEN','UGA','ZWE', 'PHL', 'GTM')

isos %>%
  purrr::map(.f = function(iso){
    shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
    shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
    
    smm <- smm %>%
      terra::crop(shpr) %>%
      terra::resample(shpr) %>%
      terra::mask(shpr)
    
    out <- paste0(root,'/data/',iso,'/agricultural_area/ag_area.tif')
    
    dir.create(dirname(out), showWarnings = F, recursive = T)
    if(!file.exists(out)){ terra::writeRaster(x = smm, out, overwrite=TRUE) }
  })

