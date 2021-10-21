options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, vegan))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

pth <- paste0(root,'/data/_global/livestock')
iso <- 'UGA'
shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
lvst <- list.files(path = pth, pattern = '*_2010_Da.tif$', full.names = T, recursive = T) %>%
  terra::rast() %>%
  terra::crop(shpr) %>%
  terra::resample(shpr) %>%
  terra::mask(shpr)
lvst <- round(lvst)

lvst_crp <- terra::app(x = lvst, fun = function(x){
  x <- as.numeric(na.omit(x))
  if(length(x) > 1){
    y <- vegan::diversity(x = x, index = "shannon", MARGIN = 1)
  } else {
    y <- NA
  }
  return(y)
})
out <- paste0(root,'/data/',iso,'/livestock/lvst_diver.tif')
dir.create(dirname(out), showWarnings = F, recursive = T)
if(!file.exists(out)){ terra::writeRaster(x = lvst_crp, out) }
