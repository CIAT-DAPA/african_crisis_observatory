# Extract daily climate data from precipitation (CHIRPS) and temperature (CPC)
# H. Achicanoy
# Alliance Bioversity-CIAT, 2021

options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(terra,plyr,tidyverse,countrycode,raster,ncdf4,sf,future,furrr,lubridate,glue,cowsay,fst))
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Main functions

get_table_wrap <- function(iso = 'SDN', country = 'Sudan'){
  # Root path
  out_dir <- paste0('//dapadfs.cgiarad.org/workspace_cluster_9/Sustainable_Food_System/',tolower(country))
  dir.create(path = out_dir, F, T)
  
  # Download shapefile 
  shp <- lowest_gadm(iso = iso, out = NULL) # Define to put out
  mxL <- grep(pattern = '^NAME_', x = names(shp@data), value = T) %>% .[length(.)] %>% readr::parse_number(.)
  shp@data <- shp@data[,grep(pattern = '^NAME_', x = names(shp@data), value = T)]
  
  path_prec <- '//catalogue/BaseLineDataCluster01/observed/gridded_products/chirps/daily'
  path_temp <- paste0('//dapadfs.cgiarad.org/workspace_cluster_9/Sustainable_Food_System/Grazia/cpc_data/5km/',tolower(country))
  
  # date <- "1989.01.01"
  # This function is reading spatial layer by layer.
  do_raster_to_table <- function(date,path_prec,path_temp,path_fst)
  {
    
    # Read data for each date: precipitation and temperatures
    a  <- raster(glue::glue('{path_prec}/chirps-v2.0.{date}.tif')) %>% raster::crop(shp)
    if(file.exists(glue::glue('{path_temp}/tmax.{date}.tif')) & file.exists(glue::glue('{path_temp}/tmin.{date}.tif')))
    {
      bc <- raster::stack(glue::glue('{path_temp}/tmax.{date}.tif'),
                          glue::glue('{path_temp}/tmin.{date}.tif')) %>% raster::crop(shp)
      abc <- raster::stack(a, bc) %>% raster::mask(shp); rm(a, bc)
      names(abc) <- c('prec', 'tmax', 'tmin')
      abc2 <- as(abc, "SpatRaster")
      ref  <- terra::rast('//dapadfs.cgiarad.org/data_cluster_2/gcm/cmip5/downscaled/ensemble/rcp26/global_30s/2020_2049/bio_1')
      shp2 <- as(shp, 'SpatVector')
      ref  <- ref %>% terra::crop(x = ., y = terra::ext(shp2))
      ref  <- as(ref, 'SpatRaster')
      shp_r <- shp_r <- terra::rasterize(x = shp2, y = ref, field = paste0('NAME_',mxL))
      abc2 <- terra::resample(x = abc2, y = ref, method = "bilinear")
      abc2 <- abc2 %>% terra::mask(mask = shp_r)
    } else {
      abc <- a
      names(abc) <- c('prec')
      abc2 <- as(abc, "SpatRaster")
      ref  <- terra::rast('//dapadfs.cgiarad.org/data_cluster_2/gcm/cmip5/downscaled/ensemble/rcp26/global_30s/2020_2049/bio_1')
      shp2 <- as(shp, 'SpatVector')
      ref  <- ref %>% terra::crop(x = ., y = terra::ext(shp2))
      ref  <- as(ref, 'SpatRaster')
      shp_r <- shp_r <- terra::rasterize(x = shp2, y = ref, field = paste0('NAME_',mxL))
      abc2 <- terra::resample(x = abc2, y = ref, method = "bilinear")
      abc2 <- abc2 %>% terra::mask(mask = shp_r)
    }
    
    # This transform stack to tibble
    points <- abc2 %>%
      terra::as.data.frame(xy = T, cell = F) %>%
      tibble::as.tibble() %>%
      dplyr::mutate(Date = date %>% lubridate::as_date()) %>%
      dplyr::select(Date, everything()) %>%
      dplyr::mutate_each(funs(replace(., . == -9999, NA_real_)), -Date, -x, -y)
    
    if(sum(colSums(is.na(points)) %in% nrow(points)) > 1){
      points <- points[,!(colSums(is.na(points)) %in% nrow(points))]
      points <- points %>%
        tidyr::drop_na() %>%
        dplyr::mutate(id = 1:nrow(.)) %>%
        dplyr::select(id, everything())
    } else {
      points <- points %>%
        tidyr::drop_na() %>%
        dplyr::mutate(id = 1:nrow(.)) %>%
        dplyr::select(id, everything())
    }
    
    pnt <- points %>% dplyr::select('x','y') %>% sp::SpatialPoints(coords = .)
    raster::crs(pnt) <- raster::crs(shp)
    
    for(var in paste0('NAME_', 0:mxL)){
      eval(parse(text = paste0('points$',var,' <- sp::over(pnt, shp, returnList = F)$',var)))
    }; rm(var)
    adm <- paste0('NAME_', 0:mxL)
    points$key <- do.call(paste, c(points[adm], sep="-"))
    
    if(c('tmax','tmin') %in% names(points)){
      points <- points %>%
        dplyr::select('Date','prec','tmax','tmin','key',adm) %>%
        dplyr::group_by('key',adm) %>%
        dplyr::summarise(tprec = sum(prec, na.rm = T),
                         aprec = mean(prec, na.rm = T),
                         mprec = median(prec, na.rm = T),
                         tmax  = mean(tmax, na.rm = T),
                         tmin  = mean(tmin, na.rm = T))
    } else {
      points <- points %>%
        dplyr::select('Date','prec','key',adm) %>%
        dplyr::group_by('key',adm) %>%
        dplyr::summarise(tprec = sum(prec, na.rm = T),
                         aprec = mean(prec, na.rm = T),
                         mprec = median(prec, na.rm = T))
      points$tmax <- NA
      points$tmin <- NA
    }
    
    return(points)
    
  }
  
  cores <- 10
  plan(cluster, workers = cores)
  options(future.globals.maxSize = 891289600)
  tbl <- tibble::tibble(Date = seq(ymd('1981-01-01'),ymd('2020-12-31'),by='day') %>%
                          stringr::str_replace_all('-', '.')) %>%
    dplyr::mutate(Climate = furrr::future_map(.x = Date, .f = do_raster_to_table, path_prec = path_prec, path_temp = path_temp, path_fst = NULL))
  gc()
  gc(reset = T)
  
  tbl2 <- tbl %>%
    tidyr::unnest()
  
  vroom::vroom_write(tbl2, paste0(out_dir,'/',tolower(country),'_daily_data_1981-2020.csv'), delim = ',')
  
}
