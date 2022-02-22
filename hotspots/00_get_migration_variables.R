# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain migration variables at 1 km resolution
# Original data source: https://sedac.ciesin.columbia.edu/data/set/popdynamics-global-est-net-migration-grids-1970-2000
# Steps:
# 1. Download manually the education rasters:
#    30arcsec-net-migration-1970-1980.tif
#    30arcsec-net-migration-1980-1990.tif
#    30arcsec-net-migration-1990-2000.tif
# 2. Execute this script to obtain:
#    Estimated Net Migration (multi-annual median)
#    Estimated Net Migration (multi-annual coefficient of variation)
#    Estimated Net Migration (multi-annual trend: Sen's slope)
# Author: Harold Achicanoy
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

# -------------------------------------- #
# Obtain migration variables
# -------------------------------------- #

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA')

get_migr_vars <- function(iso = 'SDN'){
  
  # Load the country lowest administrative level shapefile
  if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
    dir.create(path = dirname(paste0(root,'/data/',iso,'/_shps/',iso,'.shp')), recursive = TRUE)
    shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  } else {
    shp <- raster::shapefile(x = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  }
  
  # Reference raster at 1 km
  ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
  shpr <- terra::rasterize(x = shp, y = ref, field = 'key')
  
  # ----------------------------------------------------------------- #
  # Migration: mean of percentage (time series 1970-2000)
  # ----------------------------------------------------------------- #
  # Migration
  mgrt <- list.files(path = paste0(root,'/data/_global/migration'), pattern = '*.tif$', full.names = T, recursive = T) %>% terra::rast(.)
  names(mgrt) <- paste0('yr',c(1975,1985,1995))
  mgrt <- terra::crop(x = mgrt, terra::ext(shp))
  
  # ------------------------------------ #
  # Median
  # ------------------------------------ #
  
  mmgrt <- median(mgrt)
  mmgrt_crp <- terra::resample(x = mmgrt, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/migration/medn_migration.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mmgrt_crp, out) }
  
  # ------------------------------------ #
  # Coefficient of variation
  # ------------------------------------ #
  
  vmgrt <- terra::stdev(mgrt)/mean(mgrt)
  vmgrt_crp <- terra::resample(x = vmgrt, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/migration/cvar_migration.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vmgrt_crp, out) }
  
  # Crop Trend
  tmgrt <- terra::app(x = mgrt, fun = function(x){
    x <- as.numeric(na.omit(x))
    if(length(x) > 1){
      y <- trend::sens.slope(x)$estimates
    } else {
      y <- NA
    }
    return(y)
  })
  tmgrt_crp <- terra::resample(x = tmgrt, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/migration/trnd_migration.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = tmgrt_crp, out) }
  
}
isos %>%
  purrr::map(.f = function(iso){
    get_migr_vars(iso = iso)
  })
