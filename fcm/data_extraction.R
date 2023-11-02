rm(list=ls(all=TRUE))
library(terra)
library(sp)
library(mapview)
root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/'
mig <- read.csv(paste0(root,'Data/IOM_gps_dates.csv'))
head(mig)
mig <- na.omit(mig)
coordinates(mig) <- ~longitude+latitude
mig <- mig[mig$city!='Kurri',]
proj4string(mig) <-  CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") #CRS("+proj=longlat +datum=WGS84 +no_defs")
#plot(mig)
#text(mig, mig$city)
mapview(mig, cex = 3, alpha = .5, popup = NULL)

library(sf)
library(raster)
grd <- st_make_grid(st_bbox(extent(mig)+2), cellsize = 0.46, square = TRUE)
grd <- st_as_sf(grd)
grd$id <- 1:nrow(grd)
plot(grd)
mig <- st_as_sf(mig)
st_crs(grd) <- st_crs(mig)

#grd$mcount <- lengths(st_intersects(grd, mig))
#st_intersects(conflict_sf, shp), function(i){ifelse(length(i)==0,NA, i)})
test < - sapply(mig$index, function(i){lengths(st_intersects(i))})

grid$year <- NA
grid$city <- NA
grid$country <- NA
grid$mcount <- NA
for(year in unique(mig$year)){
  temp <- mig[mig$year==year,]
  grd$mcount <- lengths(st_intersects(grd, temp))
  #grid$year <- temp$year
  #grid$city <- temp$city
  #grid$country <- temp$country
}



#aa <- sapply(st_intersects(grd, mig), function(i){print(i)})
#mig[unlist(aa[32603]), ]


#which(sapply(aa, length) >0)
#aa[32603]
#mapview(list(grd[grd$id == 32603,], mig[22213, ]),        layer.name = c("grid", "mig"))



df <- sapply(st_intersects(grd, mig), function(i){if(length(i)>0){ st_drop_geometry( mig[i,]) %>% group_by(year, month) %>% tally() }else{NA} })

for(id in 1:length(df)){ if(length(df[[id]]) > 1){df[[id]] <- data.frame(id = id, df[[id]]) }}
df_final <- do.call(rbind, df)

#tt <- geodata::gadm(c('ETH',"KEN", 'UGA', 'TZA'),path = tempdir(), level=1, version='latest')