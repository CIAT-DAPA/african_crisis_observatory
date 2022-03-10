# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain recent migration variable at 1 km resolution
# Original data source: https://sedac.ciesin.columbia.edu/data/set/popdynamics-global-est-net-migration-grids-1970-2000
# Steps:
# 1. Download manually the education rasters:
#    30arcsec-net-migration-1970-1980.tif
#    30arcsec-net-migration-1980-1990.tif
#    30arcsec-net-migration-1990-2000.tif
# 2. Execute this script to obtain:
#    Recent migration variable
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

mgr <- terra::rast(paste0(root,'/data/_global/migration/1990-2000/30arcsec-net-migration-1990-2000.tif'))

isos <- c('KEN','MLI','NGA','SDN','SEN','UGA','ZWE')

isos %>%
  purrr::map(.f = function(iso){
    shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
    shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
    
    smm <- mgr %>%
      terra::crop(shpr) %>%
      terra::resample(shpr) %>%
      terra::mask(shpr)
    
    out <- paste0(root,'/data/',iso,'/migration/rcnt_migration.tif')
    
    dir.create(dirname(out), showWarnings = F, recursive = T)
    if(!file.exists(out)){ terra::writeRaster(x = smm, out) }
  })
