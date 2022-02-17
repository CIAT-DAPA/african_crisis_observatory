# R options
g <- gc(reset = T); rm(list = ls()) # Emptying the garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, exactextractr))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory' # Change for your own working path
iso <- 'KEN'

source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

# Load country shapefile
# if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
#   shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
# } else {
#   shp <- raster::shapefile(paste0(root,'/data/',iso,"/_shps/",iso,".shp" ))
# }
shp <- lowest_gadm(iso = iso)
shp <- shp %>% sf::st_as_sf() %>% dplyr::mutate(id = 1:nrow(.))

# Generate a country grid of 20 km^2
grd <- sf::st_make_grid(sf::st_bbox(raster::extent(shp)), cellsize = 0.2, square = T) %>% sf::st_as_sf(.) %>% dplyr::mutate(id = 1:nrow(.))
grd <- as(grd, 'Spatial')
raster::shapefile(grd, paste0('D:/',iso,'_grid20km.shp'))
