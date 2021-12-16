######
## @author : Andres Camilo Mendez
## script to process chirps raster to calculate metrics
## 

require(pacman)
pacman::p_load(tidyverse, raster, terra,sf, stars, fst, stringi, stringr, lubridate, furrr, purrr, future)


#######################################################
############## EXTRACT VALUES FROM CHIRPS ############
#####################################################

base_dir <- switch(Sys.info()[1],
                       "Linux" = "/cluster01/Workspace/ONECGIAR/Data/")

### select country to process
shp <- raster::shapefile(paste0(base_dir, "Africa_shp/African_continet.shp"))

country_name <- "Zimbabwe"

shp_c <- shp[shp@data$ADM0_NAME == country_name, ]

#### define raster paths
chirps_paht <- paste0(base_dir, "chirps/global_daily/tifs/p05/")
era5_path <- paste0(base_dir, "climate/era5/sis-agromet/nc/")


## get raster file names and dates

## get all info from chirps files
chirps_dirs <- list.dirs(chirps_paht, recursive = F, full.names = F)
chirps_dirs <- chirps_dirs[-length(chirps_dirs)]

chirps_file_info <- lapply(chirps_dirs, function(i){
  f_pths <- paste0(chirps_paht, i)
  r_pths <- list.files(f_pths, full.names = T, pattern = ".tif$")
  r_names <- list.files(f_pths, full.names = F, pattern = ".tif$")
  
  cat("Gettin info from year: ", i ,"\n")
  stopifnot("number of objects less than 364 " = length(r_names)>= 364)

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
  dplyr::bind_rows()


#define a raster mask from chirps data


future::plan(future::multicore, workers = 15)

chirps_file_info <- chirps_file_info %>% 
  dplyr::mutate(data_extracted = furrr::future_pmap(.l = list(.x = file_path, 
                                                              .y = file_year, 
                                                              .z = file_yday), .f = function(.x,.y,.z){
                                                                
                                                                cat("Extracting values for: ", .x, "-", paste0(.y,"_",.z ), "\n" )
                                                                layer_r <- raster(.x) %>%
                                                                  raster::crop(x = ., y = extent(shp_c)) %>% 
                                                                  raster::mask(x= ., mask = shp_c)
                                                                
                                                                r_df <- raster::as.data.frame(layer_r, xy = T)
                                                                r_df$cellid <- raster::cellFromXY(layer_r, r_df[, 1:2])
                                                                
                                                                r_df <- r_df %>% 
                                                                  tibble::as_tibble() %>% 
                                                                  dplyr::select(cellid, x, y, everything(.)) %>% 
                                                                  dplyr::na_if(-9999)
                                                                
                                                                names(r_df) <- c("cellid","x", "y", paste0(.y,"_", .z))
                                                                
                                                               return(r_df)
                                                                
                                                              }) )

saveRDS(data_extracted, "/cluster01/Workspace/ACO/1.Data/ZWE_chirps_extracted.rds")

future::plan(future::sequential)

## calcular la suma de los 120 dias mas lluviosos y caluclarles el coef var

future::plan(future::multicore, workers = 15)
temp_list <- furrr::future_map(.x = unique(data_extracted$file_year), function(.x){
  
to_process <- data_extracted %>% 
  dplyr::filter(file_year == .x) %>%
  dplyr::pull(data_extracted)%>% 
  purrr::reduce(., .f = function(x, y){
      x %>% 
      dplyr::select(-x,-y) %>% 
      right_join(., y, by = "cellid")}, .dir = "forward") %>% 
  dplyr::select(cellid, x, y, everything(.))


ncols <- to_process %>% 
  dplyr::select(-cellid,-x,-y) %>% 
  ncol()

to_process$coef_var <- to_process %>% 
  dplyr::select(-cellid,-x,-y) %>%
  apply(., 1, function(i){
    start <- 1:ncols
    end <- start + 120
    stop <- which(end == ncols)
    start <- start[1:stop]
    end <- end[1:stop]
    sums <- c()
    for(k in 1:stop){
      sums[k] <- sum(i[seq(start[k], end[k])], na.rm = T)
    }
   
   coef_var <- sd(i[seq(start[which.max(sums)], end[which.max(sums)])], na.rm = T)/mean(i[seq(start[which.max(sums)], end[which.max(sums)])], na.rm = T)
   return(coef_var)
   })

names(to_process)[names(to_process) == "coef_var"] <- paste0("coef_var_", .x)
return(to_process %>% 
         dplyr::select(cellid, x, y, starts_with("coef_var")))

}) %>% 
  purrr::reduce(., .f = function(x,y){
    x %>% 
      dplyr::select(-x,-y) %>% 
      right_join(., y, by = "cellid")
  }, .dir = "forward")%>%
  dplyr::select(cellid, x, y , everything(.)) %>% 
  write_csv(., "/home/acmendez/ZWE_chirps_rainfall_index1.csv")

future::plan(future::sequential)


## calcular la suma de los 120 dias mas lluviosos y caluclarles el la frecuencia de dias consecutivos secos


future::plan(future::multicore, workers = 15)
furrr::future_map(.x = unique(data_extracted$file_year), function(.x){
  
  to_process <- data_extracted %>% 
    dplyr::filter(file_year == .x) %>%
    dplyr::pull(data_extracted)%>% 
    purrr::reduce(., .f = function(x, y){
      x %>% 
        dplyr::select(-x,-y) %>% 
        right_join(., y, by = "cellid")}, .dir = "forward") %>% 
    dplyr::select(cellid, x, y, everything(.))
  
  
  ncols <- to_process %>% 
    dplyr::select(-cellid,-x,-y) %>% 
    ncol()
  
  to_process$sum_5dry <- to_process %>% 
    dplyr::select(-cellid,-x,-y) %>%
    apply(., 1, function(i){
      start <- 1:ncols
      end <- start + 120
      stop <- which(end == ncols)
      start <- start[1:stop]
      end <- end[1:stop]
      sums <- c()
      for(k in 1:stop){
        sums[k] <- sum(i[seq(start[k], end[k])], na.rm = T)
      }
      
      to_eval <- as.vector(i[seq(start[which.max(sums)], end[which.max(sums)])] <= 1)
      
      n_days <- rle(to_eval)$lengths[rle(to_eval)$values]
      ret <-  sum(round(n_days[n_days >= 5]/5))
      return(ret)
    })
  
  names(to_process)[names(to_process) == "sum_5dry"] <- paste0("sum_5dry_", .x)
  return(to_process %>% 
           dplyr::select(cellid, x, y, starts_with("sum_5dry")))
  
}) %>% 
  purrr::reduce(., .f = function(x,y){
    x %>% 
      dplyr::select(-x,-y) %>% 
      right_join(., y, by = "cellid")
  }, .dir = "forward")%>%
  dplyr::select(cellid, x, y , everything(.)) %>% 
  write_csv(., "/home/acmendez/ZWE_chirps_rainfall_index2.csv")

future::plan(future::sequential)





### extraer indicadores para CHIRPS
r_df <- raster::as.data.frame(vals = layer_r, xy = T)
r_df$cellid <- raster::cellFromXY(layer_r, r_df[, 1:2])

r_df <- r_df %>% 
  tibble::as_tibble() %>% 
  dplyr::select(cellid, x, y, everything(.)) %>% 
  dplyr::na_if(-9999)

names(r_df) <- c("cellid","x", "y", paste0(file_year,"_", file_yday))




#### calcular otros indicadores de para era5

r_e <- raster(paste0(era5_path,"/10m_wind_speed/Wind-Speed-10m-Mean_C3S-glob-agric_AgERA5_19991217_final-v1.0.nc")) %>% 
  raster::crop(x = ., y = extent(shp_c)) %>% 
  raster::mask(x= ., mask = shp_c)

if(file_type == "solar_radiation_flux"){
  r_e <- r_e/1000000
}else if(file_type %in% c("Temperature_Air_2m_Max_24h_", 'Temperature_Air_2m_Mean_24h_', 'Temperature_Air_2m_Min_24h_')){
  r_e <- r_e - 273.15
}

 r_e <-  raster::resample(x = r_e, y = layer_r)
 
 ret <- raster::as.data.frame(r_e, xy = T)
 names(ret) <- c("x", "y", paste0(file_year,"_", file_yday)) 

future::plan(future::multicore, workers = 15)

furrr::future_map(.x = unique(era5_file_info$file_year), .f= function(.x){ 
  
  era5_file_info %>% 
    dplyr::filter(file_year == .x) %>% 
    dplyr::pull(file_path) %>% 
    raster::Stack() %>% 
    raster::crop(., extent(shp_c))
    
  
  r_raw <- raster::raster(x)
  
  
  })


####################################################
########### EXTRAER DATOS DE ERA5   ###############
##################################################
 
base_dir <- switch(Sys.info()[1],
                   "Linux" = "/cluster01/Workspace/ONECGIAR/Data/")

### select country to process
shp <- raster::shapefile(paste0(base_dir, "Africa_shp/African_continet.shp"))

country_name <- "Zimbabwe"

shp_c <- shp[shp@data$ADM0_NAME == country_name, ]

#### define raster paths
chirps_paht <- paste0(base_dir, "chirps/global_daily/tifs/p05/")
era5_path <- paste0(base_dir, "climate/era5/sis-agromet/nc/")

layer_ref <- raster('/cluster01/Workspace/ONECGIAR/Data/chirps/global_daily/tifs/p05/1982//chirps-v2.0.1982.01.01.tif') %>% 
  raster::crop(., extent(shp_c))


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
  print(sum(pos))  
  r_date <- stringr::str_extract(r_names, '_[0-9]+') %>% 
    stringr::str_replace_all(., "_", "") %>% 
    as.Date(., format = "%Y%m%d")
  
  
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

valid_years <- era5_file_info %>% 
  filter(file_year >= "1981") %>%
  pull(file_year)%>%
  unique()

future::plan(future::multicore, workers = 20)
system.time({
  
  for(k in valid_years){
    cat("getting rastv values for year: ", k, "\n")
    
    
    
    data_extracted_era5 <- era5_file_info %>% 
      filter(file_year == k) %>%
      dplyr::mutate(data_extracted = furrr::future_pmap(.options = furrr_options(seed = TRUE), .l = list(.x = file_path, 
                                                                                                         .y = file_type, 
                                                                                                         .z = file_date), .f = function(.x,.y,.z){
                                                                                                           
                                                                                                           #cat("Extracting values for: ", .x, "-", paste0(.y,"_",.z ), "\n" )
                                                                                                           layer_r <- raster(.x) %>%
                                                                                                             raster::crop(x = ., y = extent(shp_c)) %>% 
                                                                                                             raster::mask(x= ., mask = shp_c) 
                                                                                                           
                                                                                                           if(.y == "srad"){
                                                                                                             layer_r <- layer_r/1000000 
                                                                                                           }else if(.y %in% c('tmax', 'tmean', 'tmin')){
                                                                                                             layer_r <- layer_r-273.15
                                                                                                           }
                                                                                                           layer_r <- raster::resample(layer_r, layer_ref)   
                                                                                                           
                                                                                                           r_df <- raster::as.data.frame(layer_r, xy = T) %>%
                                                                                                             drop_na()
                                                                                                           r_df$id <- raster::cellFromXY(layer_r, r_df[, 1:2])
                                                                                                           
                                                                                                           r_df <- r_df %>% 
                                                                                                             tibble::as_tibble() %>% 
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
        purrr::reduce(., .f = function(x,y){
          x %>% 
            dplyr::select(-x,-y) %>% 
            right_join(., y, by = "id")
        }, .dir = "forward")%>%
        dplyr::select( x, y, id , everything(.))%>%
        tidyr::pivot_longer(cols = -c(x,y,id), names_to = "date", values_to = var)
      
      return(ret)  
      
    }) %>% purrr::reduce(., .f = function(x,y){
      x %>% 
        dplyr::select(-x,-y, -id, -date) %>% 
        bind_cols(., y)
    }, .dir = "forward")%>%
      dplyr::select( x, y, id, date , everything(.))
    
    fst::write_fst(x = final_tbl, path = paste0( "/home/acmendez/era5_extracted/", "climate_", k, "_mod.fst"))
    
  }#end for
  
  
})#end system.time
future::plan(future::sequential)
 