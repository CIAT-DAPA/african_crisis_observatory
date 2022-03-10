# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain livestock diversity at 1 km resolution
# Original data source: http://www.fao.org/livestock-systems/global-distributions/en/
# Steps:
# 1. Download manually the livestock distribution rasters:
#    5_Bf_2010_Da.tif (buffaloes)
#    5_Ct_2010_Da.tif (cattle)
#    5_Ch_2010_Da.tif (chickens)
#    5_Dk_2010_Da.tif (ducks)
#    5_Gt_2010_Da.tif (goats)
#    5_Ho_2010_Da.tif (horses)
#    5_Pg_2010_Da.tif (pigs)
#    5_Sh_2010_Da.tif (sheeps)
# 2. Execute this script to obtain:
#    lvst_diver.tif (Shannon diversity index of livestock units)
# Author: Harold Achicanoy
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, vegan))

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'
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
