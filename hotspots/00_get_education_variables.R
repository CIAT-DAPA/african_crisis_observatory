options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

# -------------------------------------------------------- #
# Education: mean of years of education
# -------------------------------------------------------- #
# By male
mle <- terra::rast(paste0(root,'/data/_global/education/IHME_LMIC_EDU_2000_2017_MEAN_15_49_MALE_MEAN_Y2019M12D24.TIF'))
# By female
fle <- terra::rast(paste0(root,'/data/_global/education/IHME_LMIC_EDU_2000_2017_MEAN_15_49_FEMALE_MEAN_Y2019M12D24.TIF'))

# Difference of years of education
if(!file.exists(paste0(root,'/data/_global/education/difference.TIF'))){
  dff <- mle - fle
  terra::writeRaster(dff, paste0(root,'/data/_global/education/difference.TIF'))
} else {
  dff <- terra::rast(paste0(root,'/data/_global/education/difference.TIF'))
}

# Assign years names to each layer
names(mle) <- names(fle) <- names(dff) <- paste0('yr',2000:2017)

# ------------------------------------ #
# Median
# ------------------------------------ #

mmle <- median(mle)
mfle <- median(fle)
mdff <- median(dff)

# ------------------------------------ #
# Coefficient of variation
# ------------------------------------ #

vmle <- terra::stdev(mle)/mean(mle)
vfle <- terra::stdev(fle)/mean(fle)
vdff <- terra::stdev(dff)/mean(dff)

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

# Crop Median
mmle_crp <- terra::crop(x = mmle, terra::ext(shp))
mmle_crp <- terra::resample(x = mmle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
mfle_crp <- terra::crop(x = mfle, terra::ext(shp))
mfle_crp <- terra::resample(x = mfle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
mdff_crp <- terra::crop(x = mdff, terra::ext(shp))
mdff_crp <- terra::resample(x = mdff_crp, y = ref) %>% terra::mask(x = ., mask = shpr)

# Crop Coefficient of variation
vmle_crp <- terra::crop(x = vmle, terra::ext(shp))
vmle_crp <- terra::resample(x = vmle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
vfle_crp <- terra::crop(x = vfle, terra::ext(shp))
vfle_crp <- terra::resample(x = vfle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
# vdff_crp <- terra::crop(x = vdff, terra::ext(shp))

# Crop Trend
mle_crp <- terra::crop(x = mle, terra::ext(shp))
tmle_crp <- terra::app(x = mle_crp, fun = function(x){
  x <- as.numeric(na.omit(x))
  if(length(x) > 0){
    y <- trend::sens.slope(x)$estimates
  } else {
    y <- NA
  }
  return(y)
})
tmle_crp <- terra::resample(x = tmle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
fle_crp <- terra::crop(x = fle, terra::ext(shp))
tfle_crp <- terra::app(x = fle_crp, fun = function(x){
  x <- as.numeric(na.omit(x))
  if(length(x) > 0){
    y <- trend::sens.slope(x)$estimates
  } else {
    y <- NA
  }
  return(y)
})
tfle_crp <- terra::resample(x = tfle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
dff_crp <- terra::crop(x = dff, terra::ext(shp))
tdff_crp <- terra::app(x = dff_crp, fun = function(x){
  x <- as.numeric(na.omit(x))
  if(length(x) > 0){
    y <- trend::sens.slope(x)$estimates
  } else {
    y <- NA
  }
  return(y)
})
tdff_crp <- terra::resample(x = tdff_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
