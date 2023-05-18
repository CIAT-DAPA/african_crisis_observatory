
pacman::p_load(terra, ncdf4, parallel)

years <- 1958:2022


msk <- terra::rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/_global/masks/mask_world_1km.tif")
out_dir <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/_global/climate_water_deficit_terraclimate/raw_files"

r_lst <- list()

for(i in years){
  cat("downloading year:", i, "\n")
  
  ulr <- paste0('http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data/TerraClimate_def_',i,'.nc')
  out_path <- paste0(out_dir,"/cwd_", i, ".nc")
  #temp_file_out <- paste0(tempdir() ,"/cwd_", i, ".nc")
  download.file(ulr, out_path, method = "wget" )
  Sys.sleep((1+runif(1, 1,1.5)))
  
  
  
}


pacman::p_load(terra, ncdf4, parallel)

unpack <- function(file_pth, msk = msk, out_dir){
  require(ncdf4)
  cat("processing: ", file_pth, "\n")
  r <- terra::rast(file_pth)
  r <- terra::resample(r, msk )
  gc()
  out_name <- basename(terra::sources(r))
  out_name <- gsub(".nc$",".tif",out_name)
  terra::writeRaster(r, paste0(out_dir, out_name))
  terra::tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
  return(r)
}

raw_dir <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/_global/climate_water_deficit_terraclimate/raw_files"

msk <- terra::rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/_global/masks/mask_world_1km.tif")
out_dir <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/_global/climate_water_deficit_terraclimate/"
fls <- list.files(raw_dir, pattern = ".nc$", full.names = T)

terraOptions(tempdir = "W:/1.Data/Palmira/CSO/data/_global/climate_water_deficit_terraclimate/raw_files/tempfiles")

rasts <- list()
for(i in 14:65){
  rasts[[i]] <- unpack(fls[i], msk, out_dir)
}


#rasts <-base::Map(unpack, fls, msk, out_dir)

