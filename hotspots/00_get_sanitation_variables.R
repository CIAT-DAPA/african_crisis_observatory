# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain sanitation variables at 1 km resolution
# Original data source: https://cloud.ihme.washington.edu/s/bkH2X2tFQMejMxy?path=%2F
# Steps:
# 1. Download manually the sanitation rasters:
#    IHME_LMIC_WASH_2000_2017_S_PIPED_PERCENT_MEAN
#    IHME_LMIC_WASH_2000_2017_W_PIPED_PERCENT_MEAN
# 2. Execute this script to obtain:
#    Piped water (percentage) (multi-annual median)
#    Piped water (percentage) (multi-annual coefficient of variation)
#    Piped water (percentage) (multi-annual trend: Sen's slope)
#    Sanitation facilities (multi-annual median)
#    Sanitation facilities (multi-annual coefficient of variation)
#    Sanitation facilities (multi-annual trend: Sen's slope)
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
# Obtain sanitation facilities variables
# -------------------------------------- #

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA')

get_sanit_vars <- function(iso = 'SDN'){
  
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
  # Sanitation facilities: mean of percentage (time series 2000-2017)
  # ----------------------------------------------------------------- #
  # Piped water
  pp_wtr <- list.files(path = paste0(root,'/data/_global/piped_water'), pattern = '*.TIF$', full.names = T) %>% terra::rast(.)
  names(pp_wtr) <- paste0('yr',2000:2017)
  # Sanitation facilities
  st_flt <- list.files(path = paste0(root,'/data/_global/sanitation_facilities'), pattern = '*.TIF$', full.names = T) %>% terra::rast(.)
  names(st_flt) <- paste0('yr',2000:2017)
  
  # ------------------------------------ #
  # Median
  # ------------------------------------ #
  
  mpp_wtr <- median(pp_wtr)
  mst_flt <- median(st_flt)
  
  # ------------------------------------ #
  # Coefficient of variation
  # ------------------------------------ #
  
  vpp_wtr <- terra::stdev(pp_wtr)/mean(pp_wtr)
  vst_flt <- terra::stdev(st_flt)/mean(st_flt)
  
  # Crop Median
  mpp_wtr_crp <- terra::crop(x = mpp_wtr, terra::ext(shp))
  mpp_wtr_crp <- terra::resample(x = mpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/medn_piped_water.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mpp_wtr_crp, out) }
  
  mst_flt_crp <- terra::crop(x = mst_flt, terra::ext(shp))
  mst_flt_crp <- terra::resample(x = mst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/medn_sanitation_facilities.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mst_flt_crp, out) }
  
  # Crop Coefficient of variation
  vpp_wtr_crp <- terra::crop(x = vpp_wtr, terra::ext(shp))
  vpp_wtr_crp <- terra::resample(x = vpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/cvar_piped_water.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vpp_wtr_crp, out) }
  
  vst_flt_crp <- terra::crop(x = vst_flt, terra::ext(shp))
  vst_flt_crp <- terra::resample(x = vst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/cvar_sanitation_facilities.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vst_flt_crp, out) }
  
  # Crop Trend
  pp_wtr_crp <- terra::crop(x = pp_wtr, terra::ext(shp))
  tpp_wtr_crp <- terra::app(x = pp_wtr_crp, fun = function(x){
    x <- as.numeric(na.omit(x))
    if(length(x) > 1){
      y <- trend::sens.slope(x)$estimates
    } else {
      y <- NA
    }
    return(y)
  })
  tpp_wtr_crp <- terra::resample(x = tpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/trnd_piped_water.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = tpp_wtr_crp, out) }
  
  st_flt_crp <- terra::crop(x = st_flt, terra::ext(shp))
  tst_flt_crp <- terra::app(x = st_flt_crp, fun = function(x){
    x <- as.numeric(na.omit(x))
    if(length(x) > 1){
      y <- trend::sens.slope(x)$estimates
    } else {
      y <- NA
    }
    return(y)
  })
  tst_flt_crp <- terra::resample(x = tst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/trnd_sanitation_facilities.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = tst_flt_crp, out) }
  
}
isos %>%
  purrr::map(.f = function(iso){
    get_sanit_vars(iso = iso)
  })
