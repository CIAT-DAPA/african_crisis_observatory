
pacman::p_load(raster, tidyverse, sf, lubridate, stringr)




shp <- raster::shapefile(paste0("/cluster01/Workspace/ONECGIAR/Data/Africa_shp/African_continet.shp"))

country_name <- "Zimbabwe"

shp_c <- shp[shp@data$ADM0_NAME == country_name, ]

layer_ref <- raster('/cluster01/Workspace/ONECGIAR/Data/chirps/global_daily/tifs/p05/1982/chirps-v2.0.1982.01.01.tif') %>% 
  raster::crop(., extent(shp_c))%>% 
  raster::mask(., shp_c)
layer_ref[!is.na(layer_ref[])] <- 1

clim_data2 <- readRDS("/home/acmendez/climate_idex_halfway.rds")
clim_data2 <- clim_data2 %>%
  dplyr::select(-Climate)


median_idx <- clim_data2 %>%
  dplyr::mutate(idx = purrr::map(.x =  climatic_index , function(.x){
    tbl <- .x
    if( !all(is.na(tbl$time_spam)) ){
      tbl <- tbl %>%
        dplyr::select(-time_spam)
      
      ret <- apply(tbl, 2, median, na.rm = T)
      names(ret) <- names(tbl)
      
      return(ret)
      
    }
  })) %>% 
  tidyr::unnest_wider(., idx) %>%
  dplyr::select(-climatic_index)


for(i in 4:ncol(median_idx)){
  
  df <- median_idx[, c(2:3, i)]
  cat("rasterizing: ", names(df)[3], "\n")
  rst <- rasterFromXYZ(df)
  plot(rst) 
  
  writeRaster(rst, paste0("/home/acmendez/climatic_index/", names(df)[3], "_median.tif" ), overwrite = T)
  
  #coordinates(df) <-  ~x+y
  #crs(df) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  #r <- raster(resolution = 0.04166667, ext = extent(shp_country))
  #crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
  #r <- raster::mask(r, shp_country)
  #r_f <- rasterize(country_rwi, r, fun = mean, field = "AWE")
  #r_f <- raster::resample(r_f, mask %>% raster::crop(., extent(shp_country)) )
  
  #writeRaster(r_f, paste0(rwi_out_dir, "/",ISO3, "_AWE.tif"), overwrite = T)
  
}



coef_var_idx <- clim_data2 %>%
  dplyr::mutate(idx = purrr::map(.x =  climatic_index , function(.x){
    tbl <- .x
    if( !all(is.na(tbl$time_spam)) ){
      tbl <- tbl %>%
        dplyr::select(-time_spam)
      
      ret <- apply(tbl, 2, function(i){
        cv <- sd(i, na.rm=T)/mean(i, na.rm = T)
        return(cv)
      })
      names(ret) <- names(tbl)
      
      return(ret)
      
    }
  })) %>% 
  tidyr::unnest_wider(., idx) %>%
  dplyr::select(-climatic_index) 
print(coef_var_idx)



for(i in 4:ncol(coef_var_idx)){
  
  df <- coef_var_idx[, c(2:3, i)]
  cat("rasterizing: ", names(df)[3], "\n")
  rst <- rasterFromXYZ(df)
  #plot(rst) 
  
  writeRaster(rst, paste0("/home/acmendez/climatic_index/coef_var/", names(df)[3], "_cv.tif" ), overwrite = T)
  
  #coordinates(df) <-  ~x+y
  #crs(df) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  #r <- raster(resolution = 0.04166667, ext = extent(shp_country))
  #crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
  #r <- raster::mask(r, shp_country)
  #r_f <- rasterize(country_rwi, r, fun = mean, field = "AWE")
  #r_f <- raster::resample(r_f, mask %>% raster::crop(., extent(shp_country)) )
  
  #writeRaster(r_f, paste0(rwi_out_dir, "/",ISO3, "_AWE.tif"), overwrite = T)
  
}



trend_idx <- clim_data2 %>%
  dplyr::mutate(idx = purrr::map(.x =  climatic_index , function(.x){
    tbl <- .x
    if( !all(is.na(tbl$time_spam)) ){
      tbl <- tbl %>%
        dplyr::select(-time_spam)
      
      ret <- apply(tbl, 2, function(i){
        
        i <- na.omit(i)
        if(length(i) >= 2){
          y <- trend::sens.slope(i)$estimates
        }else{
          y <- NA
        }
        
        
        return(y)
      })
      
      names(ret) <- names(tbl)
      
      return(ret)
      
    }
  })) %>% 
  tidyr::unnest_wider(., idx) %>%
  dplyr::select(-climatic_index)  

print(trend_idx)


for(i in 4:ncol(trend_idx)){
  
  df <- trend_idx[, c(2:3, i)]
  cat("rasterizing: ", names(df)[3], "\n")
  rst <- rasterFromXYZ(df)
  #plot(rst) 
  
  writeRaster(rst, paste0("/home/acmendez/climatic_index/trend/", names(df)[3], "_trnd.tif" ), overwrite = T)
  
  #coordinates(df) <-  ~x+y
  #crs(df) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  #r <- raster(resolution = 0.04166667, ext = extent(shp_country))
  #crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
  #r <- raster::mask(r, shp_country)
  #r_f <- rasterize(country_rwi, r, fun = mean, field = "AWE")
  #r_f <- raster::resample(r_f, mask %>% raster::crop(., extent(shp_country)) )
  
  #writeRaster(r_f, paste0(rwi_out_dir, "/",ISO3, "_AWE.tif"), overwrite = T)
  
}



spi <- tidyft::parse_fst("/home/acmendez/climatic_index/ZWE_spis.fst") %>%
  as_tibble()

res <- spi %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(median = median(SPI, na.rm= T),
                   coef_var = sd(SPI, na.rm = T)/mean(SPI, na.rm = T),
                   tren = trend::sens.slope(na.omit(SPI))$estimates) %>%
  dplyr::left_join(., median_idx %>% dplyr::select(id, x, y), by = c("id"))

rst <- rasterFromXYZ(res %>% dplyr::select(id, x, y, SPI))







