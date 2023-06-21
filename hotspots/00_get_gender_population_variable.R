
pacman::p_load(raster)

root_dir <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data"  
iso <- "NER"

####################################################
######### female population #######################
##################################################

file_names <- paste0(tolower(iso),"_f_", c(0,1, seq(5,80, by = 5)), "_2020_constrained.tif" )

for(i in file_names){
  
  download.file(paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_Constrained/2020/",toupper(iso),"//",i), 
                
                paste0(root_dir, "/", iso, "/gender_population/", i),
                method = "wget"
  )
  Sys.sleep(30)
}

stk <- list.files(paste0(root_dir, "/", iso, "/gender_population"), pattern = "_f_", full.names = T) %>%
  purrr::map(., raster)

sum_r <- sum(raster::stack(stk), na.rm = T)
sum_r[sum_r[] <= 1] <- NA
writeRaster(sum_r,paste0(root_dir, "/", iso, "/gender_population/", iso, "female_population.tif"), overwirte = T)


#####################################
##### Male population ##############
###################################

file_names_male <- paste0(tolower(iso),"_m_", c(0,1, seq(5,80, by = 5)), "_2020_constrained.tif" )

for(i in file_names_male){
  
  download.file(paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_Constrained/2020/",toupper(iso),"//",i), 
                
                paste0(root_dir, "/", iso, "/gender_population/", i),
                method = "wget"
  )
  Sys.sleep(30)
}



stk <- list.files(paste0(root_dir, "/", iso, "/gender_population"), pattern = "_m_", full.names = T) %>%
  purrr::map(., raster)

sum_r <- sum(raster::stack(stk), na.rm = T)
sum_r[sum_r[] <= 1] <- NA
writeRaster(sum_r,paste0(root_dir, "/", iso, "/gender_population/",iso , "_male_population.tif"), overwirte = T)





