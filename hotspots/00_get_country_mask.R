# R options
g <- gc(reset = T); rm(list = ls()) # Emptying the garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, exactextractr))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
iso <- 'KEN'

# Binary land cover map: human settlements + grasslands + croplands
r <- terra::rast('D:/kenya10m.tif')

# Aggregate results to 1 km
diversity <- terra::aggregate(x = r, fact = 100, fun = 'sum', na.rm = T)
gc(reset =  T)
plot(diversity)
diversity[diversity == 0] <- NA

# shp template
shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))

# Raster template
tmp <- terra::rast(paste0(root,'/data/_global/masks/mask_world_1km.tif'))
tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(mask = shp)
plot(tmp)

diversity <- terra::resample(x = diversity, y = tmp)

out <- paste0(root,'/data/',iso,'/mask/',iso,'_mask.tif')
dir.create(path = dirname(out), F, T)
terra::writeRaster(diversity, filename = out, overwrite = T)
