options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, terra, vegan, FactoMineR))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

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
    if(!file.exists(out)){ terra::writeRaster(x = irr, out) }
  })
