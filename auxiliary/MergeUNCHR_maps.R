g <- gc(reset = T); 
rm(list = ls()) # Empty garbage collector
#root <- 'D:/OneDrive - CGIAR/Data/UNCHR/conflict'
root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data'
out <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/UNCHR'
countries_lst <- c('BDI','DJI', 'ERI', 'ETH', 'KEN','RWA','SDN','SOM', 'SSD', 'TZA', 'UGA')
region <- 'East_and_Horn_of_Africa'
shps_list <- list()
library(terra)
for(i in 1:length(countries_lst)){
  shp <- vect(paste0(root,'/',countries_lst[i],'/_results/cluster_results/conflict/conflict_regular_clust.shp'))
  labels <- read.csv(paste0(root,'/',countries_lst[i],'/_results/cluster_results/conflict/conflict_cluster_text_description.csv'))
  shp <- merge(shp, labels[c('clust','label','short_label')], by='clust')
  shps_list[[i]] <- shp
}
shp <- do.call(rbind, shps_list)
labs <- c("High", "Moderate", "Limited")
shp$cluster <- ""
for(i in 1:length(labs)){
  shp$cluster[shp$short_label==labs[i]] <- i
}
plot(shp, "label", col=rainbow(3), lwd=3, type="classes")


#files <- list.files(root, recursive = T, full.names = T, pattern = '*.shp')
#shps_list <- lapply(files, vect)

wrld_shp <- terra::vect(paste0(root, "/_global/world_shapefile/all_country/all_countries.shp"))
regionShp <- wrld_shp[wrld_shp$ISO3 %in% countries_lst,]
regionShp <- sf::st_as_sf(regionShp)
bdy <- sf::st_as_sf(shp)
library(tmap)
library(mapview)
tmap_mode("view")
map <- tm_shape(regionShp)+
  tm_borders()+
  tm_shape(bdy)+
  tm_fill(col= "cluster", palette="-YlOrRd", title='Conflict', labels = labs) +
  tm_compass(type = "8star", position = c("right", "top")) +#c('red','yellow','orange')
  tm_scale_bar(breaks = c(0, 50, 100), text.size = 1, width=1,
               position = c("left", "bottom"))+
  tm_layout(legend.outside=F, 
            legend.text.size = 1.1,
            legend.title.size= 1.3,
            legend.frame=F,
            legend.just = c("right", "bottom"), 
            #legend.width= 1,
            #legend.height= -0.2 
  )#+ 
#tm_add_legend(type = "text", text = "labels")

#tm_format("World")
map

tmap_save(map,  dpi= 600,  height=8, width=10, units="in",
          filename=paste0(out,'/', region,'.png'))

to_save <- shp[,  c("label", "short_label", "cluster")]
writeVector(to_save, filename = paste0(out,'/', region,'.shp'))
