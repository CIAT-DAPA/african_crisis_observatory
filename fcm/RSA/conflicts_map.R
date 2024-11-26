#'Author:Brenda Chepngetich
#'This script makes a conflict cluster map for Karamoja cluster
#'

library(sf)
library(terra)
library(tmap)
library(geodata)
library(dplyr)

wd <- "C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team"
setwd(wd)
#karamoja and country shapefiles
Karamoja <- st_read("C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/igad_cluster_1_1/igad_cluster_1_1.shp")
Karamoja_buffer <- st_read("C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/igad_cluster_1_1/IGAD_cluster_1_buffer_4km.shp")
uga <- st_read('C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Data/Admin/Ilemi/admin_digitized/UGANDA.shp')%>% st_set_crs(st_crs(Karamoja_buffer)) %>% dplyr::select(c('ADM0_EN', 'geometry'))
ken <- st_read('C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Brenda/PGIS_IILRI/DATA/Kenya_0.shp')%>% st_set_crs(st_crs(Karamoja_buffer)) %>% dplyr::select(c('COUNTY_NAM', 'geometry'))
colnames(ken)[colnames(ken) == "COUNTY_NAM"] <- 'ADM0_EN'
ssd <- st_read('C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Data/Admin/Ilemi/admin_digitized/SOUTH_SUDAN.shp')%>% st_set_crs(st_crs(Karamoja_buffer)) %>% dplyr::select(c('ADM0_EN', 'geometry'))
eth <- st_read('C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Data/Admin/Ilemi/admin_digitized/ETHIOPIA.shp')%>% st_set_crs(st_crs(Karamoja_buffer)) %>% dplyr::select(c('ADM0_EN', 'geometry'))
ilemi <- st_read('C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Data/Admin/Ilemi/admin_digitized/ILEMI_TRIANGLE.shp')%>% st_set_crs(st_crs(Karamoja_buffer))
gha <- st_as_sf(rbind(ssd, uga, ken, eth))

#conflicts data
ken <- geojsonio::geojson_sf(paste0(wd,"/Data/CSO/KEN/clim_conflict_ips_overlays (ACCLED-1997-2022).geojson"))%>% st_set_crs(st_crs(Karamoja))
ken_ <- ken[st_within(ken,Karamoja_buffer, sparse = FALSE),]
eth <- geojsonio::geojson_sf(paste0(wd,"/Data/CSO/ETH/clim_conflict_ips_overlays.geojson"))%>% st_set_crs(st_crs(Karamoja))
eth_ <- eth[st_within(eth,Karamoja_buffer, sparse = FALSE),]
uga <- geojsonio::geojson_sf(paste0(wd,"/Data/CSO/UGA/clim_conflict_ips_overlays.geojson"))%>% st_set_crs(st_crs(Karamoja))
uga_ <- uga[st_within(uga,Karamoja_buffer, sparse = FALSE),]
ssd <- geojsonio::geojson_sf(paste0(wd,"/Data/CSO/UGA/clim_conflict_ips_overlays.geojson"))%>% st_set_crs(st_crs(Karamoja))
plot(ken_setl)
#place names
ken_setl <- st_read('C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/korir/LULC/hotosm_ken_populated_places_points_shp.shp')%>% st_set_crs(st_crs(Karamoja)) %>% st_intersection(., Karamoja)%>% filter(., place == 'town')
uga_setl <- st_read('C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/korir/LULC/hotosm_uga_populated_places_points_shp.shp')%>% st_set_crs(st_crs(Karamoja)) %>% st_intersection(., Karamoja) %>% filter(., place == 'town')
ssd_setl <- st_read('C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/korir/LULC/SSD_PopulatedAreas_Dataset/ssd_pppls_ocha_20221216.shp' )%>% st_set_crs(st_crs(Karamoja))%>% st_intersection(., Karamoja) 



unique(uga_$conflict_clust_label)
#map
custom_palette <- c("Limited conflict" = "orange2","Moderate conflict" ="yellow", "High conflict" = "red4")
map <- tm_shape(Karamoja)+
  tm_borders(col = 'black', lwd = 1)+
  tm_shape(ken_)+
  tm_fill(col="conflict_clust_label",
          palette = custom_palette,
          title = "Conflict Cluster")+
  tm_shape(uga_)+
  tm_fill(col="conflict_clust_label", palette=custom_palette, legend.show = FALSE)+
  tm_shape(gha) + tm_borders(col = 'black', lwd = 1)+
  tm_shape(ilemi)+ tm_borders(lty = 'dotted', col = 'darkred', lwd = 2)+
  tm_shape(uga_setl)+tm_text('name', remove.overlap = T, print.tiny = T)+
  tm_shape(ken_setl)+tm_text('name', remove.overlap = T, print.tiny = T)+
  tmap::tm_legend(
    position = c(0.01, 0.02),
    legend.title.size = 1.2,
    legend.title.fontface = "bold",
    legend.text.size = 1.0
  ) +
  tmap::tm_compass(type = "arrow",
                   position = c("right", "top"),
                   size = 4) +
  tmap::tm_scale_bar(
    breaks = c(0, 50, 100),
    position = c(0.50, 0.001),
    text.size = 0.7
  ) +
  tmap::tm_graticules(
    n.x = 6,
    n.y = 6,
    lines = FALSE,
    #alpha = 0.1
    labels.rot = c(0, 90)
  ) +
  tm_layout(inner.margins = c(0.02, 0.02, 0.02, 0.02))
map  
tmap_save(map,paste0(wd,'/Brenda/conflicts mapping/karamoja_conf.png'), device = png, dpi = 600)
