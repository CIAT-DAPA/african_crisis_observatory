# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain flooding frequency raster at 1 km resolution: 
# Original data source: https://preview.grid.unep.ch/index.php?preview=data&events=floods&lang=eng
# Steps:
# 1. Download manually the flooding raster:
#    fl_frequency.tif
# 2. Execute this script to obtain:
#    Average number of flooding events per 100 years raster cropped by country
# Author: Harold Achicanoy
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, terra))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

fld <- terra::rast(paste0(root,'/data/_global/flooding/fl_frequency.tif'))

isos <- c('KEN','MLI','NGA','SDN','SEN','UGA','ZWE')

isos %>%
  purrr::map(.f = function(iso){
    shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    ref <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
    shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
    
    crd <- shpr %>% terra::as.data.frame(xy = T, cells = T, na.rm = T)
    crd <- cbind(crd, terra::extract(x = fld, y = crd[,c('x','y')]))
    names(crd)[ncol(crd)] <- 'Flood'
    crd <- crd[,c('x','y','Flood')]
    
    rst <- terra::rast(x = crd, type = 'xyz', crs = terra::crs(shp))
    
    out <- paste0(root,'/data/',iso,'/flooding/flood.tif')
    
    dir.create(dirname(out), showWarnings = F, recursive = T)
    if(!file.exists(out)){ terra::writeRaster(x = rst, out) }
  })
