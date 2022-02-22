# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain Relative Wealth Index at 5 km resolution
# Original data source: https://data.humdata.org/dataset/relative-wealth-index
# Steps:
# 1. Download manually the relative wealth index files for each country of interest:
#    [Country]_relative_wealth_index.csv
# 2. Execute this script to obtain:
#    [ISO3]_RWI.tif
# Author: Andres Mendez
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(raster,tidyverse,readxl,sf,sp))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
ISO3 <- 'UGA'

wealth_dir  <- paste0(root,'/data/_global/wealth_index')
rwi_out_dir <- paste0(root,'/data/',ISO3,'/wealth_index' )

if(!dir.exists(rwi_out_dir)){dir.create(rwi_out_dir,F,T)}

wealth_df   <- readr::read_csv(list.files(wealth_dir, pattern = paste0("^",ISO3), full.names = T))
mask        <- raster::raster(paste0(root,'/data/_global/masks/mask_world_1km.tif'))
shp_country <- raster::shapefile(paste0(root,'/data/',ISO3,'/_shps/',ISO3,'.shp'))

coordinates(wealth_df) <- ~longitude+latitude
crs(wealth_df) <- "+proj=longlat +datum=WGS84 +no_defs"

r <- raster(resolution = 0.04166667, ext = extent(shp_country))
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
r <- raster::mask(r, shp_country)
r_f <- rasterize(wealth_df, r, fun = mean, field = "rwi")
r_f <- raster::resample(r_f, mask %>% raster::crop(., extent(shp_country)) )

raster::writeRaster(r_f, paste0(rwi_out_dir,"/",ISO3,"_rwi.tif"), overwrite = T)
