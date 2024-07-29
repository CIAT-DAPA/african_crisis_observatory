###########################################################
### script to move GEE downloaded unzipped files to their
### respective folder in alliancedfs AFO folder
### 1. download and unzip GEE google drive folder
### 2. define ISO3
### 3. define GEE downloaded unzipped folder path
### 4. Run and wait



ISO <- "SOM"
GEE_data_path <-  "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/GEE data"



source_folders <- data.frame(folders =  list.dirs("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/KEN", recursive = T))
source_folders <- base::subset(source_folders,!grepl("old|season_type|_results|_shps|climatic_indexes|KEN/$|conflict|mask|livelihood|wealth_index",source_folders$folders) )
source_files <- lapply(source_folders$folders, function(i){
  unlist(list.files(i, pattern = ".tif$", full.names = T))
})
source_files <- data.frame(files_path = unlist(source_files) )
source_files$files_code = basename(source_files$files_path)


source_files$files_path <- gsub("KEN", ISO, source_files$files_path)
source_files$files_code <- gsub("KEN", ISO, source_files$files_code )


av_files <- data.frame( from = list.files(GEE_data_path, pattern = ".tif$", full.names = T))
av_files$files_code <- basename(av_files$from)



to_move <- base::merge(source_files, av_files, by.x = "files_code", by.y = "files_code", all.x = T )
to_move <- base::subset(to_move, !is.na(to_move$from))


for(i in 1:nrow(to_move)){
  cat("Moving file", to_move$files_code[i], " to:",  to_move$files_path[i], "\n")
  r <- raster::raster(to_move$from[i])
  raster::writeRaster(r, to_move$files_path[i], overwrite = T)
  
}















