# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain pasture at 1 km resolution
# Original data source: http://www.earthstat.org/cropland-pasture-area-2000/
# Steps:
# 1. Download manually the raster of pasture area from the data source
# 2. Execute this script to obtain a single rasters file with the pasture area cropped
#    by the country of interest
# Author: Harold Achicanoy
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, terra, vegan, FactoMineR))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

pas <- terra::rast(paste0(root,'/data/_global/cropland/Pasture2000_5m.tif'))

isos <- c('KEN','MLI','NGA','SDN','SEN','UGA','ZWE')

isos %>%
  purrr::map(.f = function(iso){
    shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
    shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
    
    smm <- pas %>%
      terra::crop(shpr) %>%
      terra::resample(shpr) %>%
      terra::mask(shpr)
    
    out <- paste0(root,'/data/',iso,'/pasture_area/pasture.tif')
    
    dir.create(dirname(out), showWarnings = F, recursive = T)
    if(!file.exists(out)){ terra::writeRaster(x = smm, out) }
  })
