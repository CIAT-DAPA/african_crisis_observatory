suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))

source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

# Function to compute complex Agro-climatic indices
calc_AgrClm_cmplx <- function(season = season, shp_fl = shp_fl,
                              soil_cp_pth = "",
                              soil_sat_pth = ""){
  
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
  
  # Raster template
  tmp <- terra::rast('//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif')
  shp <- terra::vect(shp_fl)
  tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
  tmp[!is.na(tmp)] <- 1
  
  if(length(season) < 12){
    cnd <- lubridate::month(tmx_dts) %in% season # Days within the season
    yrs_dts <<- split(tmx_dts[cnd],cumsum(c(1,diff(tmx_dts[cnd])!=1)))
  } else {
    yrs <- lubridate::year(tmx_dts)
    grp <- with(rle(yrs), rep(seq_along(values), lengths))
    yrs_dts <<- split(tmx_dts, grp)
  }
  
  cat('..... Computing water balance model.\n')
  WTBL <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmn <- terra::rast(tmn_fls[tmn_dts %in% yrs_dts[[i]]])
      tmn <- tmn %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmn <- tmn - 273.15
      tmn <- tmn %>% terra::resample(x = ., y = tmp) %>% terra::mask(tmp)
      tav <- terra::rast(tav_fls[tav_dts %in% yrs_dts[[i]]])
      tav <- tav %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tav <- tav - 273.15
      tav <- tav %>% terra::resample(x = ., y = tmp) %>% terra::mask(tmp)
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      tmx <- tmx %>% terra::resample(x = ., y = tmp) %>% terra::mask(tmp)
      srd <- terra::rast(srd_fls[srd_dts %in% yrs_dts[[i]]])
      srd <- srd %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      srd <- srd/1000000
      srd <- srd %>% terra::resample(x = ., y = tmp) %>% terra::mask(tmp)
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- 0
      
      # Maximum evapotranspiration
      ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)
      
      # Soil data
      scp <- terra::rast(soil_cp_pth)
      sst <- terra::rast(soil_sat_pth)
      scp <- scp %>% terra::resample(tmp) %>% terra::mask(tmp) # Soil water capacity
      sst <- sst %>% terra::resample(tmp) %>% terra::mask(tmp) # Soil water saturation point
      
      # Compute water balance model
      AVAIL <<- tmp
      AVAIL[!is.na(AVAIL)] <- 0
      watbal <- 1:terra::nlyr(ETMAX) %>%
        purrr::map(.f = function(i){
          water_balance <- eabyep_calc(soilcp  = scp,
                                       soilsat = sst,
                                       avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                       rain    = prc[[i]],
                                       evap    = ETMAX[[i]])
          AVAIL <<- water_balance$Availability
          return(water_balance)
        })
      ERATIO  <- watbal %>% purrr::map('Eratio') %>% terra::rast()
      LOGGING <- watbal %>% purrr::map('Logging') %>% terra::rast()
      IRR     <- ETMAX - prc
      GDAY    <- terra::lapp(x = terra::sds(tav, ERATIO), fun = function(TAV, ERATIO){ifelse(TAV >= 6 & ERATIO >= 0.35, 1, 0)})
      
      NDWS    <- terra::app(x = ERATIO, fun = function(ERATIO){ifelse(ERATIO < 0.5, 1, 0)}) %>% sum()
      NWLD    <- terra::app(x = LOGGING, fun = function(LOGGING){ifelse(LOGGING > 0, 1, 0)}) %>% sum()
      NWLD50  <- sum(LOGGING > (sst*0.5))
      NWLD90  <- sum(LOGGING >= sst)
      IRR     <- sum(IRR)
      
      names(NDWS) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      names(NWLD) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      names(NWLD50) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      names(NWLD90) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      names(IRR) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      
      out <- list(NDWS   = NDWS,
                  NWLD   = NWLD,
                  NWLD50 = NWLD50,
                  NWLD90 = NWLD90,
                  IRR    = IRR)
      
      return(out)
    })
  NDWS   <- WTBL %>% purrr::map('NDWS') %>% terra::rast()
  NWLD   <- WTBL %>% purrr::map('NWLD') %>% terra::rast()
  NWLD50 <- WTBL %>% purrr::map('NWLD50') %>% terra::rast()
  NWLD90 <- WTBL %>% purrr::map('NWLD90') %>% terra::rast()
  IRR    <- WTBL %>% purrr::map('IRR') %>% terra::rast()
  
  cat('..... End.\n')
  return(list(NDWS   = NDWS,
              NWLD   = NWLD,
              NWLD50 = NWLD50,
              NWLD90 = NWLD90,
              IRR    = IRR))
  
}



iso <- "ZMB"

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


shp_fl <- paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/" ,iso, "/_shps/" , iso, ".shp")
tmp_path <- "//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif"
out_root_dir <- paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/temp/")
#soil_pth <- "//catalogue/workspace_cluster_14/WFP_ClimateRiskPr/1.Data/soil/KEN"


# Loop through seasons
1:length(seasons) %>%
  purrr::map(.f = function(s){
    cat(paste0('Processing season ',names(seasons)[s],':\n'))
    # Indices calculation
    indices <- calc_AgrClm_cmplx(seasons[[s]], shp_fl,
                                 soil_cp_pth =  paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/temp/soilcp.tif"),
                                 soil_sat_pth = paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/temp/soilsat.tif"))
    # Saving results
    out <- paste0(out_root_dir,names(seasons)[s]); if(!dir.exists(out)){dir.create(out,F,T)}
    1:length(names(indices)) %>%
      purrr::map(.f = function(j){
        terra::writeRaster(x = indices[[j]], filename = paste0(out,'/',names(indices)[j],'.tif'), overwrite = T)
      })
    return(cat('Process finished successfully!\n'))
  })
