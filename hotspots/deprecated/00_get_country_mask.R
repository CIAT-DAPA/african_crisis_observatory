# --------------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain country mask raster at 1 km resolution: croplands + livestock + human settlements
# Original data source: https://developers.google.com/earth-engine/datasets/catalog/ESA_WorldCover_v100?hl=en
# Filters: 30 (Grassland), 40 (Cropland), 50 (Built-up)
# Steps:
# 1. Filter and download the country cropped raster at 10 m resolution using https://code.earthengine.google.com/b3710312e7fe6bbd615bfb90bac9160d
# 2. Merge the mosaic downloaded pieces using any GIS software (e.g. QGIS) to produce a
# single raster file at 10 m resolution
# 3. Use this script to count the number of pixels within a 1 km pixel where some of the
# filters defined are accomplished
# 4. The output corresponds to a single raster file at 1 km resolution where values higher
#    than 0 corresponds of areas of interest
#    [ISO3]_mask.tif
# Author: Harold Achicanoy
# Alliance Bioversity International - CIAT, 2022
# --------------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Emptying the garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, exactextractr))
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'
iso <- 'GTM'

# Binary land cover map: human settlements + grasslands + croplands
r <- terra::rast('D:/temp.tif') # This is a temporal result from QGIS merge processing

# Aggregate results to 1 km
diversity <- terra::aggregate(x = r, fact = 100, fun = 'sum', na.rm = T)
gc(reset =  T)
diversity[diversity == 0] <- NA
plot(diversity)

# shp template
shp <- paste0(root,'/data/',iso,'/_shps/',iso,'.shp')
if(!dir.exists(dirname(shp))){dir.create(dirname(shp),F,T)}
if(!file.exists(shp)){
  lowest_gadm(iso = iso, out = shp)
  shp <- terra::vect(shp)
} else {
  shp <- terra::vect(shp)
}

# Raster template
tmp <- terra::rast(paste0(root,'/data/_global/masks/mask_world_1km.tif'))
tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(mask = shp)

diversity <- terra::resample(x = diversity, y = tmp)

out <- paste0(root,'/data/',iso,'/mask/',iso,'_mask.tif')
dir.create(path = dirname(out), F, T)
terra::writeRaster(diversity, filename = out, overwrite = T)
