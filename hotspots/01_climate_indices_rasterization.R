rm(list=ls(all=TRUE))
pacman::p_load(tidyverse, terra, raster)

iso <- "SOM"

n_dirs <- list.files(paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/",iso, "/climatic_indexes/temp/"), pattern = "season_type") %>% 
  length()

c(paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/",iso, "/climatic_indexes/median") ,
paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/",iso, "/climatic_indexes/coef_var"),
paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/",iso, "/climatic_indexes/trend")) %>% 
  lapply(., function(pth){
    if(!dir.exists(pth)){dir.create(pth)}
  })


if(n_dirs == 1){
  fls_df <- tibble(season_type_1 = list.files(paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/",iso, "/climatic_indexes/temp/season_type_1"), recursive = T, full.names= T, pattern = ".tif$"),
                   file_name = list.files(paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso,"/climatic_indexes/temp/season_type_1"), full.names = F,  pattern = ".tif$"))
  
  
  
  fls_df %>% 
    dplyr::mutate(rst = purrr::map2(.x = season_type_1, .y = file_name, .f = function(.x, .y){
      
      cat("Processing file: ", .y, "\n")
      
      x_season_x <- terra::rast(.x )
      x_season_df <- terra::as.data.frame(x_season_x, xy = T, na.rm = T)
      
      
      medn <-  apply(x_season_df[, 3:ncol(x_season_df)], 1, function(i){median(i, na.rm = T)} )
      cv <-  apply(x_season_df[,  3:ncol(x_season_df)], 1, function(i){
        if(mean(i, na.rm = T) == 0){
          cv_fn <- 0
        }else{
          cv_fn <- sd(i, na.rm = T)/mean(i, na.rm = T)
        }
        
        return(cv_fn)
        }  )
      
      trnd <- apply(x_season_df[,  3:ncol(x_season_df)], 1, function(i){trend::sens.slope(i)$estimates})
      
      if("NWLD.tif" == .y){
        
        avg <- apply(x_season_df[, 3:ncol(x_season_df)], 1, function(i){ as.numeric(quantile(i, probs = 0.90, na.rm = T)) } )
        
        
        terra::rast(x = as.matrix(data.frame(x_season_df[, c("x", "y")], avg = avg)), type="xyz") %>% 
          terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/median/", stringr::str_replace(.y, ".tif", "_p90.tif")),overwrite=TRUE)
        
        
      }
      if("CV_prec.tif" == .y|"NDWS.tif" == .y | "THI.tif" == .y| "NTx35.tif" == .y| "CV.tif" == .y){
        avg <- apply(x_season_df[, 3:ncol(x_season_df)], 1, function(i){mean(i, na.rm = T)} )
        
        
        terra::rast(x = as.matrix(data.frame(x_season_df[, c("x", "y")], avg = avg)), type="xyz") %>% 
          terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/median/", stringr::str_replace(.y, ".tif", "_avg.tif")),overwrite=TRUE)
        
        
      }
      if("spell_5D.tif" == .y){
        avg <- apply(x_season_df[, 3:ncol(x_season_df)], 1, function(i){mean(i, na.rm = T)} )
        
        
        terra::rast(x = as.matrix(data.frame(x_season_df[, c("x", "y")], avg = avg) ), type="xyz") %>% 
          terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/median/", stringr::str_replace(.y, ".tif", "_avg.tif")),overwrite=TRUE)
        
        
      }
      
      
      
      x_season_df$medn <- medn
      x_season_df$cv <-cv
      x_season_df$trnd <- trnd
      
      terra::rast(x = as.matrix(x_season_df[, c("x", "y", "medn")]), type="xyz") %>% 
        terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/median/", stringr::str_replace(.y, ".tif", "_median.tif")),overwrite=TRUE)
      
      
      terra::rast(x = as.matrix(x_season_df[, c("x", "y", "cv")]), type="xyz") %>% 
        terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/coef_var/", stringr::str_replace(.y, ".tif", "_cv.tif")), overwrite=TRUE)
      
      terra::rast(x = as.matrix(x_season_df[, c("x", "y", "trnd")]), type="xyz") %>% 
        terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/trend/", stringr::str_replace(.y, ".tif", "_trnd.tif")), overwrite=TRUE)
      
      
      
      return("Done")
      
      
      
    }))
  
  
}else if(n_dirs == 2){
  
  fls_df <- tibble(season_type_1 = list.files(paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/",iso, "/climatic_indexes/temp/season_type_1"), recursive = T, full.names= T, pattern = ".tif$"),
                   season_type_2 = list.files(paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/temp/season_type_2"), recursive = T, full.names= T, pattern = ".tif$"),
                   file_name = list.files(paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso,"/climatic_indexes/temp/season_type_1"), full.names = F, pattern = ".tif$"))
  
  
  
  
  ##### para dos temporadas
  fls_df %>% 
    dplyr::mutate( rst =   purrr::pmap(.l = list(x = season_type_1, y = season_type_2, fl_name = file_name), function(x,y, fl_name){
      cat(">>> Processing : ", fl_name, "\n")
      .y = fl_name
      
      x_season_x <- terra::rast(purrr::map(list(x, y), terra::rast) )
      x_season_df <- terra::as.data.frame(x_season_x, xy = T, na.rm = T)
      
      
      new_order <- tibble(old_date = names(x_season_df)[-c(1,2)] ) %>% 
        dplyr::mutate(new_date =  stringr::str_replace(old_date , "X", "") %>% 
                        stringr::str_replace(., "\\.1", ".2")) %>% 
        dplyr::arrange(new_date) %>% 
        dplyr::pull(old_date)
      
      x_season_df <- x_season_df[, c("x","y", new_order)]
      
      medn <-  apply(x_season_df[, 3:ncol(x_season_df)], 1, function(i){median(i, na.rm = T)} )
      cv <-  apply(x_season_df[,  3:ncol(x_season_df)], 1, function(i){
        if(mean(i, na.rm = T) == 0){
          cv_fn <- 0
        }else{
          cv_fn <- sd(i, na.rm = T)/mean(i, na.rm = T)
        }
        
        return(cv_fn)
      }  )
      trnd <- apply(x_season_df[,  3:ncol(x_season_df)], 1, function(i){trend::sens.slope(i)$estimates})
      
      
      if("NWLD.tif" == .y){
        
        avg <- apply(x_season_df[, 3:ncol(x_season_df)], 1, function(i){ as.numeric(quantile(i, probs = 0.90, na.rm = T)) } )
        
        
        terra::rast(x = as.matrix(data.frame(x_season_df[, c("x", "y")], avg = avg)), type="xyz") %>% 
          terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/median/", stringr::str_replace(.y, ".tif", "_p90.tif")),overwrite=TRUE)
        
        
      }
      if("CV_prec.tif" == .y|"NDWS.tif" == .y | "NTx35.tif" == .y | "THI.tif" == .y | "CV.tif" == .y){
        avg <- apply(x_season_df[, 3:ncol(x_season_df)], 1, function(i){mean(i, na.rm = T)} )
        
        
        terra::rast(x = as.matrix(data.frame(x_season_df[, c("x", "y")], avg = avg)), type="xyz") %>% 
          terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/median/", stringr::str_replace(.y, ".tif", "_avg.tif")),overwrite=TRUE)
        
        
      }
      if("spell_5D.tif" == .y){
        avg <- apply(x_season_df[, 3:ncol(x_season_df)], 1, function(i){mean(i, na.rm = T)} )
        
        
        terra::rast(x = as.matrix(data.frame(x_season_df[, c("x", "y")], avg = avg) ), type="xyz") %>% 
          terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/median/", stringr::str_replace(.y, ".tif", "_avg.tif")),overwrite=TRUE)
        
        
      }
      
      
      x_season_df$medn <- medn
      x_season_df$cv <-cv
      x_season_df$trnd <- trnd
      
      terra::rast(x = as.matrix(x_season_df[, c("x", "y", "medn")]), type="xyz") %>% 
        terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/median/", stringr::str_replace(fl_name, ".tif", "_median.tif")),overwrite=TRUE)
      
      
      terra::rast(x = as.matrix(x_season_df[, c("x", "y", "cv")]), type="xyz") %>% 
        terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/coef_var/", stringr::str_replace(fl_name, ".tif", "_cv.tif")), overwrite=TRUE)
      
      terra::rast(x = as.matrix(x_season_df[, c("x", "y", "trnd")]), type="xyz") %>% 
        terra::writeRaster(., paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/trend/", stringr::str_replace(fl_name, ".tif", "_trnd.tif")), overwrite=TRUE)
      
      return("Done")
    })
    
    )
  
  
  
}else{
stop("more than two season dirs found")
}


### para una sola temporada





