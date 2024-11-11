#********
#*This script compute heat stress on human
#*Author:Brenda Chepngetich, 2024
#*******

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))

shp <- terra::vect("C:/Users/bchepngetich/Documents/Brenda/karamoja_shp/IGAD_cluster_1_buffer_4km.shp")

#Compute heat stress on humans, tmax greater than 41 degress and relative humidity greater than 50
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
  cat('..... tmax\n')
  
  # Relative humidity
  rhy_pth <- paste0(era5Dir,'/2m_relative_humidity')
  rhy_fls <- gtools::mixedsort(list.files(rhy_pth, pattern = '*.nc$', full.names = T))
  rhy_dts <- strsplit(x = rhy_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  rhy_dts <- strsplit(x = rhy_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  rhy_dts <- as.Date(rhy_dts, "%Y%m%d")
  cat('..... rh\n')
  yrs <- lubridate::year(tmx_dts)
  yrs <- names(table(yrs)[table(yrs) %in% 365:366])
  yrs <- yrs[as.integer(yrs) >= 1991 & as.integer(yrs) <= 2020]
  cat('..... yrs\n')
  tmx_fls <- tmx_fls[lubridate::year(tmx_dts) %in% yrs]
  rhy_fls <- rhy_fls[lubridate::year(rhy_dts) %in% yrs]
  tmx_dts <- tmx_dts[lubridate::year(tmx_dts) %in% yrs]
  rhy_dts <- rhy_dts[lubridate::year(rhy_dts) %in% yrs]
  cat('..... rhy\n')
  yrs <- lubridate::year(tmx_dts)
  grp <- with(rle(yrs), rep(seq_along(values), lengths))
  yrs_dts <<- split(tmx_dts, grp)
  cat("...\n")
  print(length(yrs_dts))
  cat('..... \n')
  HSI <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      cat("new\n")
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      cat('rast\n')
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      cat('crop\n')
      tmx <- tmx - 273.15
      cat('kel;vin\n')
      rhy <- terra::rast(rhy_fls[rhy_dts %in% yrs_dts[[i]]])
      cat('rhy rast\n')
      rhy <- rhy %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      cat('rhy crop\n')
      HI <- terra::lapp(x = terra::sds(tmx, rhy), fun = calc_HSH)
      cat('hi\n')
      HI <- sum(HI)
      cat('sum\n')
      names(HI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      cat('name\n')
      return(HI)
    }) %>% terra::rast()
  cat('ahem\n')
  SHI <- HSI %>% terra::mask(shp)
  cat('..... done\n')
  return (SHI)
}

Human_heat <- Human_HSI(shp)
tmp_path <- "//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif"
tmp <- terra::rast(tmp_path)
tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
tmp[!is.na(tmp)] <- 1
#resampling
HSI_resampled <- Human_heat %>% purrr::map(.f = function(r){r <- r %>% terra::resample(x = ., y = tmp) %>% terra::mask(shp); return(r)})
terra::writeRaster(HSI_resampled, filename="C:/Users/bchepngetich/Documents/Brenda/Heat stress/Human_HSI.tif",overwrite = T)
