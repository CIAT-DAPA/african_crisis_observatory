options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

# -------------------------------------------------------- #
# Child growth failure. Stunting prevalence: mean of percent
# -------------------------------------------------------- #
stntg <- list.files(path = paste0(root,'/data/_global/prevalence_stunting'), pattern = '*.TIF$', full.names = T) %>%
  terra::rast(.)

# Assign years names to each layer
names(stntg) <- paste0('yr',2000:2019)

# ------------------------------------ #
# Median
# ------------------------------------ #

mstntg <- median(stntg)

# ------------------------------------ #
# Coefficient of variation
# ------------------------------------ #

vstntg <- terra::stdev(stntg)/mean(stntg)

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

ref <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008983153,0.008983153))
shpr <- terra::rasterize(x = shp, y = ref, field = grep(pattern = '^NAME_[0-9]', x = names(shp), value = T) %>% .[length(.)])

# Crop Median
mstntg_crp <- terra::crop(x = mstntg, terra::ext(shp))
# What is missing:
#  1. Load/create a reference raster of 1 km
#  2. Rasterize shp object using the reference raster
#  3. Resampling the raster of interest to 1 km
#  4. Masking the raster of interest using the rasterized shp

# Crop Coefficient of variation
vstntg_crp <- terra::crop(x = vstntg, terra::ext(shp))

# Crop Trend
stntg_crp <- terra::crop(x = stntg, terra::ext(shp))
tstntg_crp <- terra::app(x = stntg_crp, fun = function(x){
  x <- as.numeric(na.omit(x))
  if(length(x) > 0){
    y <- trend::sens.slope(x)$estimates
  } else {
    y <- NA
  }
  return(y)
})

# -------------------------------------------------------- #
# Child growth failure. Underweight prevalence: mean of percent
# -------------------------------------------------------- #
uwght <- list.files(path = paste0(root,'/data/_global/prevalence_underweight'), pattern = '*.TIF$', full.names = T) %>%
  terra::rast(.)

# Assign years names to each layer
names(uwght) <- paste0('yr',2000:2019)

# ------------------------------------ #
# Median
# ------------------------------------ #

muwght <- median(uwght)

# ------------------------------------ #
# Coefficient of variation
# ------------------------------------ #

vuwght <- terra::stdev(uwght)/mean(uwght)

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
muwght_crp <- terra::crop(x = muwght, terra::ext(shp))
# What is missing:
#  1. Load/create a reference raster of 1 km
#  2. Rasterize shp object using the reference raster
#  3. Resampling the raster of interest to 1 km
#  4. Masking the raster of interest using the rasterized shp

# Crop Coefficient of variation
vuwght_crp <- terra::crop(x = vuwght, terra::ext(shp))

# Crop Trend
uwght_crp <- terra::crop(x = uwght, terra::ext(shp))
tuwght_crp <- terra::app(x = uwght_crp, fun = function(x){
  x <- as.numeric(na.omit(x))
  if(length(x) > 0){
    y <- trend::sens.slope(x)$estimates
  } else {
    y <- NA
  }
  return(y)
})

# -------------------------------------------------------- #
# Child growth failure. Wasting prevalence: mean of percent
# -------------------------------------------------------- #
wstng <- list.files(path = paste0(root,'/data/_global/prevalence_wasting'), pattern = '*.TIF$', full.names = T) %>%
  terra::rast(.)

# Assign years names to each layer
names(wstng) <- paste0('yr',2000:2019)

# ------------------------------------ #
# Median
# ------------------------------------ #

mwstng <- median(wstng)

# ------------------------------------ #
# Coefficient of variation
# ------------------------------------ #

vwstng <- terra::stdev(wstng)/mean(wstng)

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
mwstng_crp <- terra::crop(x = mwstng, terra::ext(shp))
# What is missing:
#  1. Load/create a reference raster of 1 km
#  2. Rasterize shp object using the reference raster
#  3. Resampling the raster of interest to 1 km
#  4. Masking the raster of interest using the rasterized shp

# Crop Coefficient of variation
vwstng_crp <- terra::crop(x = vwstng, terra::ext(shp))

# Crop Trend
wstng_crp <- terra::crop(x = wstng, terra::ext(shp))
twstng_crp <- terra::app(x = wstng_crp, fun = function(x){
  x <- as.numeric(na.omit(x))
  if(length(x) > 0){
    y <- trend::sens.slope(x)$estimates
  } else {
    y <- NA
  }
  return(y)
})
