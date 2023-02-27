g <- gc(reset = T); rm(list = ls())
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))

source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

get_5_dry_spell <- function(prec){
  
  n_days <- rle(prec)$lengths[rle(prec)$values]
  ret <-  sum(round(n_days[n_days >= 5]/5), na.rm = T)
  return(ret)
}


coef_var_index <- function(prec){
  
  ret <- sd(prec, na.rm=T)/mean(prec, na.rm = T)
  return(ret)
}


calc_SHIMP <- function(TMAX, RH){SHI <- ifelse(TMAX >= 29 & RH > 50, 1, 0); return(SHI)}; calc_SHIMP <- compiler::cmpfun(calc_SHIMP)
calc_THIMP <- function(TMAX, RH){THI <- (1.8 * TMAX + 32) - ((0.55 - 0.0055 * RH) * (1.8 * TMAX - 26.8)); THI_n <- ifelse(test = THI >= 79 & THI < 89 | THI >= 89, yes = 1, no = 0); return(THI_n)}; calc_THIMP <- compiler::cmpfun(calc_THIMP)
calc_HSIMP <- function(TMAX, RH){
  # 2 ~ danger, 3 ~ emergency
  # HSI_n <- ifelse((TMAX <= 26 & RH > 70 |
  #                 TMAX <= 27 & RH >= 40 & RH < 85 |
  #                 TMAX <= 28 & RH < 85 |
  #                 TMAX <= 29 & RH < 60 |
  #                 TMAX <= 30 & RH < 40), 1, 0)
  HSI_n <- ifelse((TMAX <= 27 & RH >= 85 |
                     TMAX <= 28 & RH >= 85 |
                     TMAX <= 29 & RH >= 60 |
                     TMAX <= 30 & RH >= 40 |
                     TMAX > 30), 1, 0)
  return(HSI_n)
}; calc_HSIMP <- compiler::cmpfun(calc_HSIMP)


# Function to compute basic Agro-climatic indices
calc_AgrClm <- function(season = season, shp_fl = shp_fl){
  
  ## ROI: regions of interest
  shp <- terra::vect(shp_fl)
  
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
  
  if(length(season) < 12){
    cnd <- lubridate::month(tmx_dts) %in% season # Days within the season
    yrs_dts <<- split(tmx_dts[cnd],cumsum(c(1,diff(tmx_dts[cnd])!=1)))
  } else {
    yrs <- lubridate::year(tmx_dts)
    grp <- with(rle(yrs), rep(seq_along(values), lengths))
    yrs_dts <<- split(tmx_dts, grp)
  }
  
  cat('..... Computing: AT. Average temperature.\n')
  AT <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tav <- terra::rast(tav_fls[tav_dts %in% yrs_dts[[i]]])
      tav <- tav %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tav <- tav - 273.15
      AT <- terra::app(x = tav, fun = function(x){ y = mean(x, na.rm = T); return(y) })
      names(AT) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(AT)
    }) %>% terra::rast()
  AT <- AT %>% terra::mask(shp)
  
  cat('..... Computing: TR. Total rainfall.\n')
  TR <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      cat(i, "\n")
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      TR <- terra::app(x = prc, fun = function(x){ y = sum(x, na.rm = T); return(y) })
      names(TR) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(TR)
    }) %>% terra::rast()
  TR <- TR %>% terra::mask(shp)
  
  cat('..... Computing: spell_5D. Frequency of 5 dry spell.\n')
  spell_5D <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      TR <- terra::app(x = prc, fun = function(x){ y = get_5_dry_spell(x); return(y) })
      names(TR) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(TR)
    }) %>% terra::rast()
  spell_5D <- spell_5D %>% terra::mask(shp)
  
  
  cat('..... Computing: Coef var within rainy seasons.\n')
  CV_rainfall <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      TR <- terra::app(x = prc, fun = function(x){ y = coef_var_index(x); return(y) })
      names(TR) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(TR)
    }) %>% terra::rast()
  CV_rainfall <- CV_rainfall %>% terra::mask(shp)
  
  cat('..... Computing: NTx35. Number of days with Tmax above 35C.\n')
  NTx35 <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      NTx35 <- terra::app(x = tmx, fun = function(x){ y = calc_htsCMP(tmax = x, t_thresh = 35); return(y) })
      names(NTx35) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(NTx35)
    }) %>% terra::rast()
  NTx35 <- NTx35 %>% terra::mask(shp)
  
  cat('..... Computing: CDD. Maximum number of consecutive dry days.\n')
  CDD <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      CDD <- terra::app(x = prc, fun = function(x){ y = calc_cddCMP(PREC = x, p_thresh = 1); return(y) })
      names(CDD) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(CDD)
    }) %>% terra::rast()
  CDD <- CDD %>% terra::mask(shp)
  
  cat('..... Computing: P5D. Rolling precipitation average of 5 days.\n')
  P5D <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      P5D <- terra::app(x = prc, fun = function(x){ y = calc_p5dCMP(PREC = x); return(y) })
      names(P5D) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(P5D)
    }) %>% terra::rast()
  P5D <- P5D %>% terra::mask(shp)
  
  cat('..... Computing: P95. Percentile 95% of daily precipitation.\n')
  P95 <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      P95 <- terra::app(x = prc, fun = function(x){ y = calc_p95CMP(PREC = x); return(y) })
      names(P95) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(P95)
    }) %>% terra::rast()
  P95 <- P95 %>% terra::mask(shp)
  
  cat('..... Computing: CSDI. Cold spell duration Index.\n')
  CSDI <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmn <- terra::rast(tmn_fls[tmn_dts %in% yrs_dts[[i]]])
      tmn <- tmn %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmn <- tmn - 273.15
      CSDI <- terra::app(x = tmn, fun = function(x){ y = calc_csdiMP(TMIN = x); return(y) })
      names(CSDI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(CSDI)
    }) %>% terra::rast()
  CSDI <- CSDI %>% terra::mask(shp)
  
  cat('..... Computing: SHI. Number of days with maximum temperature > 29C and relative humidity > 50% (sheep)\n')
  SHI <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      rhy <- terra::rast(rhy_fls[rhy_dts %in% yrs_dts[[i]]])
      rhy <- rhy %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      SHI <- terra::lapp(x = terra::sds(tmx, rhy), fun = calc_SHIMP)
      SHI <- sum(SHI)
      names(SHI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(SHI)
    }) %>% terra::rast()
  SHI <- SHI %>% terra::mask(shp)
  
  cat('..... Computing: THI. Thermal humidity index (cattle).\n')
  THI <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      rhy <- terra::rast(rhy_fls[rhy_dts %in% yrs_dts[[i]]])
      rhy <- rhy %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      THI <- terra::lapp(x = terra::sds(tmx, rhy), fun = calc_THIMP)
      THI <- sum(THI)/terra::nlyr(tmx)
      names(THI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(THI)
    }) %>% terra::rast()
  THI <- THI %>% terra::mask(shp)
  
  cat('..... Computing: HSI. Heat stress index (pig).\n')
  HSI <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      rhy <- terra::rast(rhy_fls[rhy_dts %in% yrs_dts[[i]]])
      rhy <- rhy %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      HSI <- terra::lapp(x = terra::sds(tmx, rhy), fun = calc_HSIMP)
      HSI <- sum(HSI)/terra::nlyr(tmx)
      names(HSI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(HSI)
    }) %>% terra::rast()
  HSI <- HSI %>% terra::mask(shp)
  
  cat('..... End.\n')
  return(list(AT    = AT,
              TR    = TR,
              NTx35 = NTx35,
              CDD   = CDD,
              P5D   = P5D,
              P95   = P95,
              CSDI  = CSDI,
              SHI   = SHI,
              THI   = THI,
              HSI   = HSI,
              spell_5D = spell_5D,
              CV = CV_rainfall))
  
}

iso <- "GTM"

seasons <- switch(iso, "KEN" = list(season_type_1 = 1:6, season_type_2 = 7:12),
                  "SEN" = list(season_type_1 = 6:12),
                  "ETH" = list(season_type_1 = 4:10),
                  "GTM" = list(season_type_1 = 5:11),
                  "MLI" = list(season_type_1 = 6:12),
                  "NGA" = list(season_type_1 = 6:12),
                  "PHL" = list(season_type_1 = 6:12),
                  "SDN" = list(season_type_1 = 6:10),
                  "UGA" = list(season_type_1 = 3:10),
                  "ZMB" = list(season_type_1 = 1:6, season_type_2 = 7:12),
                  "ZWE" = list(season_type_1 = 1:6, season_type_2 = 7:12)
                  )

shp_fl <- paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/_shps/", iso, ".shp")
tmp_path <- "//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif"
out_root_dir <- paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/temp/")

1:length(seasons) %>%
  purrr::map(.f = function(s){
    cat(paste0('Processing season ',names(seasons)[s],':\n'))
    # Indices calculation
    indices <- calc_AgrClm(seasons[[s]], shp_fl)
    # Load 5 km raster template
    tmp <- terra::rast(tmp_path)
    shp <- terra::vect(shp_fl)
    tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
    tmp[!is.na(tmp)] <- 1
    # Indices resampling
    indices <- indices %>% purrr::map(.f = function(r){r <- r %>% terra::resample(x = ., y = tmp) %>% terra::mask(shp); return(r)})
    # Saving results
    out <- paste0(out_root_dir, names(seasons)[s]); if(!dir.exists(out)){dir.create(out,F,T)}
    1:length(names(indices)) %>%
      purrr::map(.f = function(j){
        terra::writeRaster(x = indices[[j]], filename = paste0(out,'/',names(indices)[j],'.tif'), overwrite = T)
      })
    return(cat('Process finished successfully!\n'))
  })






