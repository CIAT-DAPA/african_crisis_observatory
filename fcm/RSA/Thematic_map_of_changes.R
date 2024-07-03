#* Land Use Land Cover change analysis to inform resource sharing agreement in Karamoja
#* 
#* 
#* Author:: Victor Korir
####################################################################################
library(sf)
library(terra)
library(tidyverse)
library(raster)
library(tmap)
library(OpenLand)
library(geodata)
library(RColorBrewer)


Kenya_wards <- geodata::gadm('KENYA', level=3, path = tempdir(), version = 'latest') 
uganda_wards <- geodata::gadm('Uganda', level=3, path = tempdir(), version = 'latest') 
cluster_wards <- rbind(st_as_sf(Kenya_wards), st_as_sf(uganda_wards))
Turkana_wards <- st_as_sf(Kenya_wards%>% subset(Kenya_wards$NAME_1 == 'Turkana',))
Kenya <-  geodata::gadm('KENYA', level=0, path = tempdir(), version = 'latest') 
Uganda <- geodata::gadm('Uganda', level=0, path = tempdir(), version = 'latest') 
Countries <- st_as_sf(rbind(Kenya, Uganda))

AOI <- st_read("D:/OneDrive - CGIAR/SA_Team/korir/LULC/igad_cluster_1_1/igad_cluster_1_1.shp")
AOI_buffer <- st_buffer(AOI, 0.44996400287491034)
plot(AOI)
cluster_wards <- cluster_wards[st_within(cluster_wards, AOI, sparse = F),]

#ESA_CCI Land Cover
CCI_1995 <- terra::rast("D:/OneDrive - CGIAR/SA_Team/korir/LULC/ESA_CCI/ESACCI-LC-L4-LCCS-Map-300m-P1Y-1995-v2.0.7cds.nc")
CCI_2004 <- terra::rast("D:/OneDrive - CGIAR/SA_Team/korir/LULC/ESA_CCI/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2004-v2.0.7cds.nc")
CCI_2013 <- terra::rast("D:/OneDrive - CGIAR/SA_Team/korir/LULC/ESA_CCI/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2013-v2.0.7cds.nc")
CCI_2007 <- terra::rast("D:/OneDrive - CGIAR/SA_Team/korir/LULC/ESA_CCI/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2007-v2.0.7cds.nc")
CCI_2012 <- terra::rast("D:/OneDrive - CGIAR/SA_Team/korir/LULC/ESA_CCI/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2012-v2.0.7cds.nc")
CCI_2017 <- terra::rast("D:/OneDrive - CGIAR/SA_Team/korir/LULC/ESA_CCI/C3S-LC-L4-LCCS-Map-300m-P1Y-2017-v2.1.1.nc")
CCI_2022 <- terra::rast("D:/OneDrive - CGIAR/SA_Team/korir/LULC/ESA_CCI/C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc")



#Cropping and masking to the AOI
CCI_1995 <- CCI_1995 %>% terra::crop(., st_bbox(AOI)) %>% terra::mask(., sf::st_as_sf(AOI))
CCI_2004 <- CCI_2004 %>% terra::crop(., st_bbox(AOI)) %>% terra::mask(., sf::st_as_sf(AOI))
CCI_2013 <- CCI_2013 %>% terra::crop(., st_bbox(AOI)) %>% terra::mask(., sf::st_as_sf(AOI))
CCI_2007 <- CCI_2007 %>% terra::crop(., st_bbox(AOI)) %>% terra::mask(., sf::st_as_sf(AOI))
CCI_2012 <- CCI_2012 %>% terra::crop(., st_bbox(AOI)) %>% terra::mask(., sf::st_as_sf(AOI))
CCI_2017 <- CCI_2017 %>% terra::crop(., st_bbox(AOI)) %>% terra::mask(., sf::st_as_sf(AOI))
CCI_2022 <- CCI_2022 %>% terra::crop(., st_bbox(AOI)) %>% terra::mask(., sf::st_as_sf(AOI))


#Defining reclass matrix based on the IPCC classification guidelines
#*1-cropland,  2-Forest, 3- grassland,
#* 4- Settlements, 5- wetland, 6- others
#* 

class_names <- c("1" = "Cropland", "2" = "Forest", "3" = "Grassland", "4" = "Settlements", "5" = "Wetland", "6" = "Others")
Reclass_mat <- matrix(c(10, 1, 11, 1, 20, 1, 30, 1, 40, 1, 60, 2, 61, 2, 50, 2, 
                        62, 2,100, 2,110, 3, 122, 3, 120, 3, 130, 3, 150, 3, 151,
                        3, 152,3, 153,3, 160,2,170,2, 180, 5, 190, 4, 200, 6,201, 6, 202, 6, 210, 6),
                      ncol = 2,
                      byrow = TRUE)
#Reclassifying

LULC1995_reclass <- terra::classify(CCI_1995, Reclass_mat)
LULC2004_reclass <- terra::classify(CCI_2004, Reclass_mat)
LULC2013_reclass <- terra::classify(CCI_2013, Reclass_mat)
LULC2022_reclass <- terra::classify(CCI_2022, Reclass_mat)


#Masking grassland
Grassland_95 <- LULC1995_reclass$lccs_class == 3
Grassland_04 <- LULC2004_reclass$lccs_class == 3
Grassland_13 <- LULC2013_reclass$lccs_class == 3
Grassland_22 <- LULC2022_reclass$lccs_class == 3

#COmputing the changes
change_95_04 <- Grassland_04 - Grassland_95
change_04_13 <- Grassland_13 - Grassland_04
change_22_13 <- Grassland_22 - Grassland_13
#change_22_95 <- Grassland_22 - Grassland_95

ward_changes1 <- extract(c(change_95_04,change_04_13,change_22_13), cluster_wards, fun = sum, na.rm = T)
colnames(ward_changes1)<- c('ID', 'Change_1995_2004', 'Change_2004_2013', 'Change_2013_2022')
#ward_changes2 <- extract(change_04_13, cluster_wards, fun = sum, na.rm = T)
#ward_changes3 <- extract(change_22_13, cluster_wards, fun = sum, na.rm = T)
#ward_changes4 <- extract(change_22_95, cluster_wards, fun = sum, na.rm = T)

cluster_wards <- cbind(cluster_wards, ward_changes1)


cluster_wards <- pivot_longer(cluster_wards, cols = c( 'Change_1995_2004', 'Change_2004_2013', 'Change_2013_2022'), names_to = "Period", values_to = "Change")
cluster_wards$Change <- cluster_wards$Change*9

#Thematic maps
them_map <- tm_shape(cluster_wards) +
  tm_polygons("Change", style = 'pretty',palette ='YlOrBr', title ='Grassland Change(Ha)')+
  tm_facets('Period', nrow = 2)+
  tm_shape(Countries)+tm_borders(col = 'black', lwd = 2)+
  tmap::tm_legend(
    position = c(0.01, 0.01),
    legend.title.size = 1.2,
    legend.title.fontface = "bold",
    legend.text.size = 0.8
  ) +
  tmap::tm_compass(type = "8star",
                   position = c("right", "top"),
                   size = 1) +
  tmap::tm_scale_bar(
    breaks = c(seq(0, 40, 10)),
    position = c(0.76, 0.001),
    text.size = 0.6
  ) +
  
  tmap::tm_graticules(
    n.x = 6,
    n.y = 6,
    lines = FALSE,
    #alpha = 0.1
    labels.rot = c(0, 90)
  ) +
  tmap::tm_layout(inner.margins = c(0.02, 0.02, 0.02, 0.02))

tmap_save(them_map,'Grassland Change.png',width=1920, height=1080, asp=0, dpi = 300)
