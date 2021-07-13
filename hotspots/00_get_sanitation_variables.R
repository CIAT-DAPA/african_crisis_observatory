options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

# -------------------------------------------------------- #
# Piped water: mean of percent
# -------------------------------------------------------- #
pp_wtr <- list.files(path = paste0(root,'/data/_global/piped_water'), pattern = '*.TIF$', full.names = T) %>%
  terra::rast(.)

# Assign years names to each layer
names(pp_wtr) <- paste0('yr',2000:2017)

# ------------------------------------ #
# Median
# ------------------------------------ #

mpp_wtr <- median(pp_wtr)

# ------------------------------------ #
# Coefficient of variation
# ------------------------------------ #

vpp_wtr <- terra::stdev(pp_wtr)/mean(pp_wtr)

# ------------------------------------ #
# Process by country
# ------------------------------------ #

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA')
iso <- 'SDN'

if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
  dir.create(path = dirname(paste0(root,'/data/',iso,'/_shps/',iso,'.shp')), recursive = TRUE)
  shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
  shp <- as(shp, 'SpatVector')
} else {
  shp <- terra::vect(x = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
}

# Crop Median
mpp_wtr_crp <- terra::crop(x = mpp_wtr, terra::ext(shp))
# What is missing:
#  1. Load/create a reference raster of 1 km
#  2. Rasterize shp object using the reference raster
#  3. Resampling the raster of interest to 1 km
#  4. Masking the raster of interest using the rasterized shp

# Crop Coefficient of variation
vpp_wtr_crp <- terra::crop(x = vpp_wtr, terra::ext(shp))

# Crop Trend
pp_wtr_crp <- terra::crop(x = pp_wtr, terra::ext(shp))
tpp_wtr_crp <- terra::app(x = pp_wtr_crp, fun = function(x){
  x <- as.numeric(na.omit(x))
  if(length(x) > 0){
    y <- trend::sens.slope(x)$estimates
  } else {
    y <- NA
  }
  return(y)
})

# -------------------------------------------------------- #
# Sanitation facilities: mean of percent
# -------------------------------------------------------- #
st_flt <- list.files(path = paste0(root,'/data/_global/sanitation_facilities'), pattern = '*.TIF$', full.names = T) %>%
  terra::rast(.)

# Assign years names to each layer
names(st_flt) <- paste0('yr',2000:2017)

# ------------------------------------ #
# Median
# ------------------------------------ #

mst_flt <- median(st_flt)

# ------------------------------------ #
# Coefficient of variation
# ------------------------------------ #

vst_flt <- terra::stdev(st_flt)/mean(st_flt)

# ------------------------------------ #
# Process by country
# ------------------------------------ #

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA')
iso <- 'SDN'

if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
  dir.create(path = dirname(paste0(root,'/data/',iso,'/_shps/',iso,'.shp')), recursive = TRUE)
  shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
  shp <- as(shp, 'SpatVector')
} else {
  shp <- terra::vect(x = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
}

# Crop Median
mst_flt_crp <- terra::crop(x = mst_flt, terra::ext(shp))
# What is missing:
#  1. Load/create a reference raster of 1 km
#  2. Rasterize shp object using the reference raster
#  3. Resampling the raster of interest to 1 km
#  4. Masking the raster of interest using the rasterized shp

# Crop Coefficient of variation
vst_flt_crp <- terra::crop(x = vst_flt, terra::ext(shp))

# Crop Trend
st_flt_crp <- terra::crop(x = st_flt, terra::ext(shp))
tst_flt_crp <- terra::app(x = st_flt_crp, fun = function(x){
  x <- as.numeric(na.omit(x))
  if(length(x) > 0){
    y <- trend::sens.slope(x)$estimates
  } else {
    y <- NA
  }
  return(y)
})
