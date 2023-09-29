# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain irrigation raster at 1 km resolution: 
# Original data source: 
# Steps:
# 1. Download manually the irrigation raster:
#    irrigation.tif
# 2. Execute this script to obtain:
#    Irrigation raster cropped by country
# Author: Harold Achicanoy
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, terra, vegan, FactoMineR))

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'

irr <- terra::rast(paste0(root,'/data/_global/irrigation/irrigation.tif'))

isos <- c('KEN','MLI','NGA','SDN','SEN','UGA','ZWE')

isos %>%
  purrr::map(.f = function(iso){
    shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
    shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
    
    irr <- irr %>%
      terra::crop(shpr) %>%
      terra::resample(shpr) %>%
      terra::mask(shpr)
    
    out <- paste0(root,'/data/',iso,'/irrigation/irrigation.tif')
    
    dir.create(dirname(out), showWarnings = F, recursive = T)
    if(!file.exists(out)){ terra::writeRaster(x = irr, out, overwrite = T) }
  })
