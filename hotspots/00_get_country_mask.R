# R options
g <- gc(reset = T); rm(list = ls()) # Emptying the garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, exactextractr))
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
iso <- 'GTM'

# Binary land cover map: human settlements + grasslands + croplands
r <- terra::rast('D:/temp.tif')

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
