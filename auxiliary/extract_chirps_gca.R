# ***********
# Author: Benson kenduiywo 
# Date: 11/11/2024
# Extract Daily CHIRPS rainfall and aggregate per admin unit

rm(list=ls(all=TRUE))
library(terra)
library(sf)
path <- "//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/Chirps/"
kenya <- vect("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/raw/admin/KEN/Kenya_Ward_WGS.shp")
county <- terra::subset(kenya, kenya$COUNTY_NAM=="TRANS NZOIA")
files <- list.files(path, pattern=paste0("^chirps-v2.*",".+.",".+.","*.tif"), full.names=TRUE)
#Filter dates 
start <- "1995.01.01"
files <- files[gsub("chirps-v2.0.","", gsub(".tif", "" ,basename(files))) >= start]
head(files)
tail(files)

loadimg <- function(file){
  temp <- terra::rast(file)
  temp <- terra::crop(temp, county)
  return(temp)
}
start.time <- Sys.time()
imgs <- purrr::map(files, loadimg)
#imgs <- terra::lapp(as.vector(files), loadimg, cores=parallel::detectCores()-2)
imgs <- terra::rast(imgs)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#Create grids of approximately 10km2 and extract rainfall
library(sf)
county <- sf::st_as_sf(county)
grd <- st_make_grid(st_bbox(terra::ext(county)), cellsize = 0.09, square =  T) %>% 
  st_as_sf(.) %>%
  dplyr::mutate(id = 1:nrow(.))
st_crs(grd) <- st_crs(county)
temp_g <- st_intersection(grd, county[,c("NAME", "COUNTY_NAM")])
temp_g <- na.omit(temp_g)
index <- unique(temp_g$id)
temp_g$Ward <- ""
for(i in 1:length(index)){
  temp_g$Ward[temp_g$id==index[i]] <- paste(temp_g$NAME[temp_g$id==index[i]], collapse = '/ ')
}
temp_g <- temp_g[!duplicated(temp_g$id),]
#library(tidyverse)

system.time(df <- terra::extract(imgs, grd, fun=mean, bind=T))
dff <- as.data.frame(df)

#Transform to long
dff <- reshape(dff,
               direction = "long",
               varying = list(names(df)[-1]), #
               v.names = "Rainfall",
               idvar = "id",
               timevar = "Date",
               times= names(df)[-1]
)

dff$Date <- gsub("chirps-v2.0.", "", dff$Date)
rownames(dff) <- NULL
temp  <- as.data.frame(temp_g[,c("id", "Ward")])

index <- unique(temp$id)
dff$Ward <- ""
for(i in 1:length(index)) {
  mindex <- which(dff$id == index[i])
  if(length(mindex) > 0) {
    dff$Ward[mindex] <- rep(temp$Ward[temp$id==index[i]], length(mindex))
  }
}

head(dff)
# #Purrmapmin paralle
# ncores <- 20
# future::plan(cluster, workers = ncores, gc = T)
# 1:nrow(stp) |>
#   furrr::future_map(.f = function(i){calc_mprec(yr = stp$yrs[i], mn = stp$mns[i]); gc(T)})
# future:::ClusterRegistry('stop')


#write.csv(dff, "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/results/gca/Trans-Nzoia_Chirps_Rainfall.csv")           
saveRDS(dff, "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/results/gca/Trans-Nzoia_Chirps_Rainfall.rds")

#library(foreign)
#write.dta(dff, "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/results/gca/Trans-Nzoia_Chirps_Rainfall.dta")           


#Tes parallel processing

# 
# library(doParallel)
# start.time <- Sys.time()
# cl <- makeCluster(20)
# registerDoParallel(cl)
# imgs <- foreach(i=files,
#                     .packages = "terra") %dopar% {
#                       crp <- wrap(terra::rast(i))
#                       #crp <- wrap(terra::crop(crp, county))
#                       return(crp)
#                     }
# 
# imgs <- lapply(imgs, unwrap)
# parallel::stopCluster(cl)
# end.time <- Sys.time()
# end.time - start.time

