# Up-sampling of SRTM raster from 1 km to 5 km
# SRTM (Shuttle Radar Topography Mission): elevation model
# H. Achicanoy
# Alliance Bioversity-CIAT, 2021

# R options and load packages
options(warn = -1, scipen = 999)
suppressPackageStartupMessages(suppressMessages(pacman::p_load(tidyverse, raster, sf, RSAGA)))

# Root directory
root <- '//dapadfs.cgiarad.org/.../african_crisis_observatory'

# CHIRPS data (ref raster 5 km)
chr_tmp <- raster('//catalogue/BaseLineDataCluster01/observed/gridded_products/chirps/daily/chirps-v2.0.1981.01.01.tif')
chr_tmp[which(chr_tmp[] == -9999)] <- NA
chr_tmp[which(chr_tmp[] < 0)] <- 0

# Country info
iso     <- 'SOM'
country <- 'Somalia'

shp_rdc <- raster::getData(name = 'GADM', country = iso, level = 0, download = T)

# SRTM
srtm <- raster('//dapadfs.cgiarad.org/data_cluster_4/observed/gridded_products/srtm/SRTM_v41_30s/srtm_v41_30s')
srtm <- raster::crop(x = srtm, y = raster::extent(shp_rdc))
srtm_5km <- raster::resample(srtm, chr_tmp, method = 'bilinear')
srtm_5km <- srtm_5km %>%
  raster::crop(., extent(shp_rdc)) %>%
  raster::mask(., mask = shp_rdc)
raster::writeRaster(srtm_5km, paste0(root,'/cpc_data/srtm/',tolower(country),'/srtm_5km.tif'), overwrite = T)
