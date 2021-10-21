options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, vegan))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

pth <- paste0(root,'/data/_global/livestock')
iso <- 'KEN'
shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern', SIMPLIFY = T)
lvst <- list.files(path = pth, pattern = '*_2010_Da.tif$', full.names = T, recursive = T) %>%
  grep2(pattern = c('cattle','sheep','goats','chicken','pig','horse'), x = ., value = T) %>%
  terra::rast() %>%
  terra::crop(shpr) %>%
  terra::resample(shpr) %>%
  terra::mask(shpr)

lvst_crp <- terra::app(x = lvst, fun = function(x){
  x <- as.numeric(x)
  x[which(is.na(x))] <- 0
  y <- 0.7*x[1] + 0.1*x[2] + 0.1*x[3] + 0.01*x[4] + 0.2*x[5] + 0.8*x[6]
  y <- ifelse(y == 0, NA, y)
  return(y)
})
out <- paste0(root,'/data/',iso,'/livestock/lvst_tlu.tif')
dir.create(dirname(out), showWarnings = F, recursive = T)
if(!file.exists(out)){ terra::writeRaster(x = lvst_crp, out) }

# TLU: tropical livestock units
# TLU = Cattle * 0.7 + Sheep * 0.1 + Goat * 0.1 + Chicken * 0.01 + Pig * 0.2 + Horse * 0.8
# Percentile: 0.1

# + vacas - recursos
# - vacas + pobreza
