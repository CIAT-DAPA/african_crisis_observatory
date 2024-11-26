#********
#*This script compute heat stress on cattle
#*Author:Brenda Chepngetich, 2024
#*******
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))

shp <- terra::vect("C:/Users/bchepngetich/Documents/Brenda/karamoja_shp/IGAD_cluster_1_buffer_4km.shp")

Calc_HSC <- function(tmax, rhum){
  thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
  return(thi)
}
# 

Cattle_HSC <- function(shp){
  # Tmax
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
  print(length(yrs_dts))
  HSC <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      rhy <- terra::rast(rhy_fls[rhy_dts %in% yrs_dts[[i]]])
      rhy <- rhy %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      HI <- terra::lapp(x = terra::sds(tmx, rhy), fun = Calc_HSC)
      HI <- sum(HI)/terra::nlyr(tmx)
      names(HI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(HI)
    }) %>% terra::rast()
  SHI <- HSC %>% terra::mask(shp)
  return (SHI)
}

Cattle_heat <- Cattle_HSC(shp)
Cattle_heat <- mean(Cattle_heat)
tmp_path <- "//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif"
tmp <- terra::rast(tmp_path)
tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
tmp[!is.na(tmp)] <- 1
Cattle_heat <- terra::resample(Cattle_heat, tmp)
terra::writeRaster(Cattle_heat, filename="C:/Users/bchepngetich/Documents/Brenda/Heat stress/Cattle_HSI.tif",overwrite = T)


#future computations
rhy_pth <- paste0(era5Dir,'/2m_relative_humidity')
rhy_fls <- gtools::mixedsort(list.files(rhy_pth, pattern = '*.nc$', full.names = T))
rhy_dts <- strsplit(x = rhy_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
rhy_dts <- strsplit(x = rhy_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
rhy_dts <- as.Date(rhy_dts, "%Y%m%d")

yrs <- lubridate::year(rhy_dts)
yrs <- yrs[yrs >= "1991" & yrs <= "2020"]
rhy_fls <- rhy_fls[lubridate::year(rhy_fls) %in% yrs]
rhy_dts <- rhy_dts[lubridate::year(rhy_dts) %in% yrs]
monthly_averages <- list()

for (year in 1991:2020){
  for (month in 1:12){
    month_files <- rhy_fls[format(rhy_dts, "%Y") == as.character(year) & format(rhy_dts, "%m") == sprintf("%02d", month)]
    if (length(month_files) > 0) {
      rasters <- rast(month_files)
      monthly_mean <- mean(rasters, na.rm = TRUE)
      monthly_averages[[paste(year, month, sep = "_")]] <- monthly_mean
    }
  }}

  
rhy_monthly <- terra::rast(monthly_averages) 
terra::writeRaster(rhy_monthly,filename="C:/Users/bchepngetich/Documents/Brenda/RH_Monthly.tif",overwrite = T)

rh <- terra::rast("C:/Users/bchepngetich/Documents/Brenda/RH_Monthly.tif")
tmax <- tmax
future_HSC <- lapp(list(tmax, rh), fun = Calc_HSC)
future_HSC <- mean(future_HSC)

