suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))

source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

zscore <- function(x){
  return((x-mean(x, na.rm=T))/(sd(x, na.rm=T)))
}

get_5_dry_spell <- function(prec){
  
  n_days <- rle(prec)$lengths[rle(prec)$values]
  ret <-  sum(round(n_days[n_days >= 5]/5), na.rm = T)
  return(ret)
}


coef_var_index <- function(prec){
  
  ret <- sd(prec, na.rm=T)/mean(prec, na.rm = T)
  return(ret)
}

iso <- 'KEN'

## ROI: regions of interest
shp <- geoda  

## Daily files
# Precipitation
chr_pth <- '//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/Chirps'
chr_fls <- gtools::mixedsort(list.files(chr_pth, pattern = '*.tif$', full.names = T))
chr_dts <- strsplit(x = chr_fls, split = 'chirps-v2.0.', fixed = T) %>% purrr::map(2) %>% unlist()
chr_dts <- strsplit(x = chr_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
chr_dts <- as.Date(gsub('.', '-', chr_dts, fixed = T))

# Tmax
era5Dir <- '//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/ERA5'
tmx_pth <- paste0(era5Dir,'/2m_temperature-24_hour_maximum')
tmx_fls <- gtools::mixedsort(list.files(tmx_pth, pattern = '*.nc$', full.names = T))
tmx_dts <- strsplit(x = tmx_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
tmx_dts <- strsplit(x = tmx_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
tmx_dts <- as.Date(tmx_dts, "%Y%m%d")

# Tmin
tmn_pth <- paste0(era5Dir,'/2m_temperature-24_hour_minimum')
tmn_fls <- gtools::mixedsort(list.files(tmn_pth, pattern = '*.nc$', full.names = T))
tmn_dts <- strsplit(x = tmn_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
tmn_dts <- strsplit(x = tmn_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
tmn_dts <- as.Date(tmn_dts, "%Y%m%d")

# Tmean
tav_pth <- paste0(era5Dir,'/2m_temperature-24_hour_mean')
tav_fls <- gtools::mixedsort(list.files(tav_pth, pattern = '*.nc$', full.names = T))
tav_dts <- strsplit(x = tav_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
tav_dts <- strsplit(x = tav_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
tav_dts <- as.Date(tav_dts, "%Y%m%d")

# Solar radiation
srd_pth <- paste0(era5Dir,'/solar_radiation_flux')
srd_fls <- gtools::mixedsort(list.files(srd_pth, pattern = '*.nc$', full.names = T))
srd_dts <- strsplit(x = srd_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
srd_dts <- strsplit(x = srd_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
srd_dts <- as.Date(srd_dts, "%Y%m%d")

# Relative humidity
rhy_pth <- paste0(era5Dir,'/2m_relative_humidity')
rhy_fls <- gtools::mixedsort(list.files(rhy_pth, pattern = '*.nc$', full.names = T))
rhy_dts <- strsplit(x = rhy_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
rhy_dts <- strsplit(x = rhy_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
rhy_dts <- as.Date(rhy_dts, "%Y%m%d")

# Filtering days within the season
yrs <- lubridate::year(tmx_dts)
yrs <- names(table(yrs)[table(yrs) %in% 365:366])
yrs <- yrs[!yrs %in% c("2020", "2021")]

tmx_fls <- tmx_fls[lubridate::year(tmx_dts) %in% yrs]
tmn_fls <- tmn_fls[lubridate::year(tmn_dts) %in% yrs]
tav_fls <- tav_fls[lubridate::year(tav_dts) %in% yrs]
srd_fls <- srd_fls[lubridate::year(srd_dts) %in% yrs]
rhy_fls <- rhy_fls[lubridate::year(rhy_dts) %in% yrs]

tmx_dts <- tmx_dts[lubridate::year(tmx_dts) %in% yrs]
tmn_dts <- tmn_dts[lubridate::year(tmn_dts) %in% yrs]
tav_dts <- tav_dts[lubridate::year(tav_dts) %in% yrs]
srd_dts <- srd_dts[lubridate::year(srd_dts) %in% yrs]
rhy_dts <- rhy_dts[lubridate::year(rhy_dts) %in% yrs]





