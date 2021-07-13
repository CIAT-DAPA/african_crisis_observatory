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

# Crop Median
mpp_wtr_crp <- terra::crop(x = mpp_wtr, terra::ext(shp))
mpp_wtr_crp <- terra::resample(x = mpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

# Crop Coefficient of variation
vpp_wtr_crp <- terra::crop(x = vpp_wtr, terra::ext(shp))
vpp_wtr_crp <- terra::resample(x = vpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

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
tpp_wtr_crp <- terra::resample(x = tpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

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

# Crop Median
mst_flt_crp <- terra::crop(x = mst_flt, terra::ext(shp))
mst_flt_crp <- terra::resample(x = mst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

# Crop Coefficient of variation
vst_flt_crp <- terra::crop(x = vst_flt, terra::ext(shp))
vst_flt_crp <- terra::resample(x = vst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

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
tst_flt_crp <- terra::resample(x = tst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
