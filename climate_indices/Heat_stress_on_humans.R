#********
#*This script compute heat stress on human
#*Author:Brenda Chepngetich, 2024
#*******

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))

source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')
shp <- terra::vect("")

#Compute heat stress on humans, tmax greater than 41 degress and relative humidity greater than 50
calc_HSI <- function(TMAX, RH){
  HSI <- ifelse(TMAX >= 41 & RH > 50, 1, 0)
  return(HSI)}
calc_HSI <- compiler::cmpfun(calc_HSI)

Human_HSI <- function(shp){
  # Tmax
  #era5Dir <- '//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/ERA5'
  era5Dir <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5'
  tmx_pth <- paste0(era5Dir,'/2m_temperature-24_hour_maximum')
  tmx_fls <- gtools::mixedsort(list.files(tmx_pth, pattern = '*.nc$', full.names = T))
  tmx_dts <- strsplit(x = tmx_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  tmx_dts <- strsplit(x = tmx_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  tmx_dts <- as.Date(tmx_dts, "%Y%m%d")
  
  # Relative humidity
  rhy_pth <- paste0(era5Dir,'/2m_relative_humidity')
  rhy_fls <- gtools::mixedsort(list.files(rhy_pth, pattern = '*.nc$', full.names = T))
  rhy_dts <- strsplit(x = rhy_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  rhy_dts <- strsplit(x = rhy_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  rhy_dts <- as.Date(rhy_dts, "%Y%m%d")
  
  yrs <- lubridate::year(tmx_dts)
  yrs <- names(table(yrs)[table(yrs) %in% 365:366])
  yrs <- yrs[as.integer(yrs) >= 1991 & as.integer(yrs) <= 2020]
  
  tmx_fls <- tmx_fls[lubridate::year(tmx_dts) %in% yrs]
  rhy_fls <- rhy_fls[lubridate::year(rhy_dts) %in% yrs]
  tmx_dts <- tmx_dts[lubridate::year(tmx_dts) %in% yrs]
  rhy_dts <- rhy_dts[lubridate::year(rhy_dts) %in% yrs]
  
  yrs <- lubridate::year(tmx_dts)
  grp <- with(rle(yrs), rep(seq_along(values), lengths))
  yrs_dts <<- split(tmx_dts, grp)
  
  SHI <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      rhy <- terra::rast(rhy_fls[rhy_dts %in% yrs_dts[[i]]])
      rhy <- rhy %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      SHI <- terra::lapp(x = terra::sds(tmx, rhy), fun = calc_HSI)
      SHI <- sum(SHI)
      names(SHI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(SHI)
    }) %>% terra::rast()
  SHI <- SHI %>% terra::mask(shp)
  return (SHI)
}

HSI <- Human_HSI(shp)
tmp_path <- "//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif"
tmp <- terra::rast(tmp_path)
tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
tmp[!is.na(tmp)] <- 1
#resampling
HSI_resampled <- HSI %>% purrr::map(.f = function(r){r <- r %>% terra::resample(x = ., y = tmp) %>% terra::mask(shp); return(r)})
terra::writeRaster(HSI_resampled, filename="",overwrite = T)
