options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

# ------------------------------------ #
# Process by country
# ------------------------------------ #

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA')
iso <- 'SDN'

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

# shp of reference preparation
ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
shpr <- terra::rasterize(x = shp, y = ref, field = 'key')

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

# Crop Median
mstntg_crp <- terra::crop(x = mstntg, terra::ext(shp))
mstntg_crp <- terra::resample(x = mstntg_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

# Crop Coefficient of variation
vstntg_crp <- terra::crop(x = vstntg, terra::ext(shp))
vstntg_crp <- terra::resample(x = vstntg_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

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
tstntg_crp <- terra::resample(x = tstntg_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

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

# Crop Median
muwght_crp <- terra::crop(x = muwght, terra::ext(shp))
muwght_crp <- terra::resample(x = muwght_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

# Crop Coefficient of variation
vuwght_crp <- terra::crop(x = vuwght, terra::ext(shp))
vuwght_crp <- terra::resample(x = vuwght_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

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
tuwght_crp <- terra::resample(x = tuwght_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

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

# Crop Median
mwstng_crp <- terra::crop(x = mwstng, terra::ext(shp))
mwstng_crp <- terra::resample(x = mwstng_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

# Crop Coefficient of variation
vwstng_crp <- terra::crop(x = vwstng, terra::ext(shp))
vwstng_crp <- terra::resample(x = vwstng_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

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
twstng_crp <- terra::resample(x = twstng_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
