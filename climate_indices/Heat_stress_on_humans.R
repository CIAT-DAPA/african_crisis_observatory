#********
#*This script compute heat stress on human
#*Author:Brenda Chepngetich, 2024
#*******
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))

shp <- terra::vect("C:/Users/bchepngetich/Documents/Brenda/karamoja_shp/IGAD_cluster_1_buffer_4km.shp")

era5Dir <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/AgERA5'
# Constants
c1 = -8.78469475556
c2 =  1.61139411
c3 =  2.33854883889
c4 = -0.14611605
c5 = -0.012308094
c6 = -0.0164248277778
c7 =  2.211732 * 10^(-3)
c8 =  7.2546 * 10^(-4)
c9 = -3.582 * 10^(-6)

calc_HSH <- function(tmean, rhum){
  hi <- ifelse(tmean >= 25,
               c1 + (c2*tmean) + (c3*rhum) + (c4*tmean*rhum) + (c5*tmean^2) + (c6*rhum^2) + (c7*tmean^2*rhum) + (c8*tmean*rhum^2) + (c9*tmean^2*rhum^2),
               tmean)
  return(hi)
}

Human_HSI <- function(shp){
  # Tmax
  era5Dir <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5'
  tmx_pth <- paste0(era5Dir,'/2m_temperature-24_hour_maximum')
  tmx_fls <- gtools::mixedsort(list.files(tmx_pth, pattern = '*.nc$', full.names = T))
  tmx_dts <- strsplit(x = tmx_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  tmx_dts <- strsplit(x = tmx_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  tmx_dts <- as.Date(tmx_dts, "%Y%m%d")


  # Tmean
  tav_pth <- paste0(era5Dir,'/2m_temperature-24_hour_mean')
  tav_fls <- gtools::mixedsort(list.files(tav_pth, pattern = '*.nc$', full.names = T))
  tav_dts <- strsplit(x = tav_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  tav_dts <- strsplit(x = tav_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  tav_dts <- as.Date(tav_dts, "%Y%m%d")
  
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
  
  yrs <- lubridate::year(tav_dts)
  yrs <- names(table(yrs)[table(yrs) %in% 365:366])
  yrs <- yrs[as.integer(yrs) >= 1991 & as.integer(yrs) <= 2020]

  tav_fls <- tav_fls[lubridate::year(tav_dts) %in% yrs]
  rhy_fls <- rhy_fls[lubridate::year(rhy_dts) %in% yrs]
  tav_dts <- tav_dts[lubridate::year(tav_dts) %in% yrs]
  rhy_dts <- rhy_dts[lubridate::year(rhy_dts) %in% yrs]

  yrs <- lubridate::year(tmx_dts)
  grp <- with(rle(yrs), rep(seq_along(values), lengths))
  yrs_dts <<- split(tmx_dts, grp)
  
  yrs <- lubridate::year(tav_dts)
  grp <- with(rle(yrs), rep(seq_along(values), lengths))
  yrs_dts <<- split(tav_dts, grp)
  print(length(yrs_dts))
  HSI <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      tav <- terra::rast(tav_fls[tav_dts %in% yrs_dts[[i]]])
      tav <- tav %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tav <- tav - 273.15
      rhy <- terra::rast(rhy_fls[rhy_dts %in% yrs_dts[[i]]])
      rhy <- rhy %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      HI <- terra::lapp(x = terra::sds(tmx, rhy), fun = calc_HSH)
      HI <- sum(HI)/terra::nlyr(tmx)
      HI <- terra::lapp(x = terra::sds(tav, rhy), fun = calc_HSH)
      HI <- sum(HI)/terra::nlyr(tav)
      names(HI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(HI)
    }) %>% terra::rast()
  SHI <- HSI %>% terra::mask(shp)
  return (SHI)
}

Human_heat <- Human_HSI(shp)
HH <- mean(Human_heat)
tmp_path <- "//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif"
tmp <- terra::rast(tmp_path)
tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
tmp[!is.na(tmp)] <- 1
HH <- terra::resample(HH, tmp)
terra::writeRaster(HH, filename="C:/Users/bchepngetich/Documents/Brenda/Heat stress/Human_HSI.tif",overwrite = T)

