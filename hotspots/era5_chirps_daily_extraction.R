# 00

require(pacman)
pacman::p_load(tidyverse, raster, terra,sf, stars, fst, stringi, stringr, lubridate, furrr, purrr, future, ncdf4)

base_dir <- switch(Sys.info()[1],
                   "Linux" = "/cluster01/Workspace/ONECGIAR/Data/")

### select country to process
shp <- raster::shapefile(paste0(base_dir, "Africa_shp/African_continet.shp"))

country_name <- "Zimbabwe"

shp_c <- shp[shp@data$ADM0_NAME == country_name, ]

#### define raster paths
chirps_paht <- paste0(base_dir, "chirps/global_daily/tifs/p05/")
era5_path <- paste0(base_dir, "climate/era5/sis-agromet/nc/")

layer_ref <- raster('/cluster01/Workspace/ONECGIAR/Data/chirps/global_daily/tifs/p05/1982/chirps-v2.0.1982.01.01.tif') %>% 
  raster::crop(., extent(shp_c))%>% 
  raster::mask(., shp_c)
layer_ref[!is.na(layer_ref[])] <- 1

layer_ref_ids <- raster::as.data.frame(layer_ref, xy = T) %>%
  dplyr::mutate(id = raster::cellFromXY(layer_ref, .[, 1:2])) %>%
  drop_na()


era5_file_names <- data.frame(dir_name =c("2m_relative_humidity", "solar_radiation_flux", "2m_temperature", "2m_temperature", "2m_temperature", "10m_wind_speed"  ) ,
                              harold_name = c('rh', 'srad','tmax', 'tmean', 'tmin', 'wind'),
                              pattern = c("2m", "Flux", "Max", "Mean", "Min", "Mean"))

era5_file_info <- apply(era5_file_names, 1, function(i){
  
  f_pths <- paste0(era5_path, i[1])
  r_pths <- list.files(f_pths, full.names = T, pattern = ".nc$")
  
  pos <- grepl(i[3],   r_pths )
  r_pths <- r_pths[pos]
  r_names <- list.files(f_pths, full.names = F, pattern = ".nc$")
  r_names <- r_names[pos]
  
  cat(i[1], "\n")
  stopifnot("number of objects in era5 are less than 364 " = length(r_names)>= 364)
  
  r_date <- stringr::str_extract(r_names, '_[0-9]+') %>% 
    stringr::str_replace_all(., "_", "") %>% 
    as.Date(., format = "%Y%m%d")
  cat(i, "\n")
  
  
  ret <- tibble(file_path = r_pths, 
                file_name = r_names, 
                file_date = r_date,
                file_year = lubridate::year(r_date),
                file_month = lubridate::month(r_date),
                file_day = lubridate::day(r_date),
                file_yday = lubridate::yday(r_date),
                file_type = i[2])
  return(ret)
}) %>% 
  dplyr::bind_rows()

chirps_dirs <- list.dirs(chirps_paht, recursive = F, full.names = F)
chirps_dirs <- chirps_dirs[-length(chirps_dirs)]

chirps_file_info <- lapply(chirps_dirs, function(i){
  
  f_pths <- paste0(chirps_paht, i)
  r_pths <- list.files(f_pths, full.names = T, pattern = ".tif$")
  r_names <- list.files(f_pths, full.names = F, pattern = ".tif$")
  
  cat(i, "\n")
  stopifnot("number of objects in chirps is less than 364 " = length(r_names)>= 364)
  
  r_date <- stringr::str_extract(r_names, '[0-9]{4}\\.[0-9]{2}\\.[0-9]{2}') %>% 
    stringr::str_replace_all(., "\\.", "") %>% 
    as.Date(., format = "%Y%m%d")
  
  ret <- tibble(file_path = r_pths, 
                file_name = r_names, 
                file_date = r_date,
                file_year = lubridate::year(r_date),
                file_month = lubridate::month(r_date),
                file_day = lubridate::day(r_date),
                file_yday = lubridate::yday(r_date))
  
  return(ret)
  
})%>% 
  dplyr::bind_rows() %>%
  dplyr::mutate(file_type = "prec")



era5_chirps_file_info <- bind_rows(era5_file_info,
                                   chirps_file_info)




chirp_dates <- chirps_file_info %>% 
  dplyr::filter(file_year == 2016) %>%
  pull(file_date) %>% unique()

data_extracted_era5 <- era5_chirps_file_info %>% 
  dplyr::filter(file_year == 2016) %>%
  dplyr::filter(file_date %in% chirp_dates) %>%
  dplyr::mutate(data_extracted = furrr::future_pmap(.options = furrr_options(seed = TRUE), .l = list(.x = file_path, 
                                                                                                     .y = file_type, 
                                                                                                     .z = file_date), .f = function(.x,.y,.z){
                                                                                                       
                                                                                                       #cat("Extracting values for: ", .x, "-", paste0(.y,"_",.z ), "\n" )
                                                                                                       layer_r <- raster(.x) %>%
                                                                                                         raster::crop(x = ., y = extent(shp_c)) 
                                                                                                       
                                                                                                       if(.y == "srad"){
                                                                                                         layer_r <- layer_r/1000000 
                                                                                                       }else if(.y %in% c('tmax', 'tmean', 'tmin')){
                                                                                                         layer_r <- layer_r-273.15
                                                                                                       }
                                                                                                       
                                                                                                       if(.y != "prec"){
                                                                                                         layer_r <- raster::resample(layer_r, layer_ref)   
                                                                                                       }
                                                                                                       
                                                                                                       
                                                                                                       r_df <- raster::as.data.frame(layer_r, xy = T) 
                                                                                                       r_df$id <- raster::cellFromXY(layer_r, r_df[, 1:2])
                                                                                                       
                                                                                                       r_df <- r_df %>% 
                                                                                                         tibble::as_tibble() %>% 
                                                                                                         dplyr::filter(id %in% layer_ref_ids$id ) %>% 
                                                                                                         dplyr::select( x, y,id, everything(.)) %>% 
                                                                                                         dplyr::na_if(-9999)
                                                                                                       
                                                                                                       names(r_df) <- c("x", "y", "id",  as.character(.z))
                                                                                                       
                                                                                                       
                                                                                                       return(r_df)
                                                                                                       
                                                                                                     }) ) %>%
  dplyr::group_by(file_type) %>%
  dplyr::group_split()

era5_chirps_file_info %>% 
  filter(file_year == 1982, file_type == "prec") %>%
  head()


#valid_years <- era5_chirps_file_info %>% 
#filter(file_year >= "1981") %>%
#pull(file_year)%>%
#unique()
valid_years <- 2016:2020

future::plan(future::multicore, workers = 20)
system.time({
  
  for(k in valid_years){
    cat("getting rastv values for year: ", k, "\n")
    
    ### anotación: para algunos años chirps tiene 364 días lo que afecta el procesamiento puesto que era5 tiene los dias completos
    chirp_dates <- chirps_file_info %>% 
      dplyr::filter(file_year == k) %>%
      pull(file_date) %>% unique()
    
    
    data_extracted_era5 <- era5_chirps_file_info %>% 
      dplyr::filter(file_year == k) %>%
      dplyr::filter(file_date %in% chirp_dates) %>%
      dplyr::mutate(data_extracted = furrr::future_pmap(.options = furrr_options(seed = TRUE), .l = list(.x = file_path, 
                                                                                                         .y = file_type, 
                                                                                                         .z = file_date), .f = function(.x,.y,.z){
                                                                                                           
                                                                                                           #cat("Extracting values for: ", .x, "-", paste0(.y,"_",.z ), "\n" )
                                                                                                           layer_r <- raster(.x) %>%
                                                                                                             raster::crop(x = ., y = extent(shp_c)) 
                                                                                                           
                                                                                                           if(.y == "srad"){
                                                                                                             layer_r <- layer_r/1000000 
                                                                                                           }else if(.y %in% c('tmax', 'tmean', 'tmin')){
                                                                                                             layer_r <- layer_r-273.15
                                                                                                           }
                                                                                                           
                                                                                                           if(.y != "prec"){
                                                                                                             layer_r <- raster::resample(layer_r, layer_ref)   
                                                                                                           }
                                                                                                           
                                                                                                           
                                                                                                           r_df <- raster::as.data.frame(layer_r, xy = T) 
                                                                                                           r_df$id <- raster::cellFromXY(layer_r, r_df[, 1:2])
                                                                                                           
                                                                                                           r_df <- r_df %>% 
                                                                                                             tibble::as_tibble() %>% 
                                                                                                             dplyr::filter(id %in% layer_ref_ids$id ) %>% 
                                                                                                             dplyr::select( x, y,id, everything(.)) %>% 
                                                                                                             dplyr::na_if(-9999)
                                                                                                           
                                                                                                           names(r_df) <- c("x", "y", "id",  as.character(.z))
                                                                                                           
                                                                                                           
                                                                                                           return(r_df)
                                                                                                           
                                                                                                         }) ) %>%
      dplyr::group_by(file_type) %>%
      dplyr::group_split()
    
    final_tbl <-   lapply(data_extracted_era5, function(l){
      var  <- unique(l$file_type)
      
      ret <- l %>%
        dplyr::pull(data_extracted) %>%
        purrr::reduce(., left_join, by = c("x", "y", "id")) %>%
        dplyr::select( x, y, id , everything(.))%>%
        tidyr::pivot_longer(cols = -c(x,y,id), names_to = "date", values_to = var)
      
      return(ret)  
      
    }) 
    
    nrow_check <- sapply(final_tbl, nrow)
    stopifnot("Differents row numbers. " = all(x))
    
    final_tbl<- final_tbl %>% 
      purrr::reduce(., left_join, by =  c("x", "y", "id", "date"))%>%
      dplyr::select( x, y, id, date , everything(.), prec)
    
    fst::write_fst(x = final_tbl, path = paste0( "/home/acmendez/era5_extracted/", "climate_", k, "_mod.fst"))
    
  }#end for
  
  
})#end system.time
future::plan(future::sequential)


