#* LULC Transition Matrix
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


#Kenya_wards <- geodata::gadm('KENYA', level=2, path = tempdir(), version = 'latest') 
#Turkana_wards <- Kenya_wards%>% subset(Kenya_wards$NAME_1 == 'Turkana',)
AOI <- st_read("D:/OneDrive - CGIAR/SA_Team/korir/LULC/igad_cluster_1_1/igad_cluster_1_1.shp")
AOI_buffer <- st_buffer(AOI, 0.44996400287491034)
plot(AOI)


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


kara1_class <- c(LULC1995_reclass$lccs_class, LULC2004_reclass$lccs_class, LULC2013_reclass$lccs_class, 
                 LULC2022_reclass$lccs_class)

# Ward level transitions
kara_class_df <- as.data.frame(kara1_class, na.rm = TRUE)
colnames(kara_class_df) <- c('lc1995', 'lc2004', 'lc2013', 'lc2022')

create_transition_matrix <- function(raster_a, raster_b) {
  transition_matrix <- table(raster_a, raster_b)
  transition_matrix_percent <- prop.table(transition_matrix, 1) * 100
  return(as.data.frame.matrix(transition_matrix_percent))
}

transition_matrix_95_04 <- create_transition_matrix(kara_class_df$lc1995, kara_class_df$lc2004)
transition_matrix_04_13 <- create_transition_matrix(kara_class_df$lc2004, kara_class_df$lc2013)
transition_matrix_13_22 <- create_transition_matrix(kara_class_df$lc2013, kara_class_df$lc2022)
transition_matrix_95_22 <- create_transition_matrix(kara_class_df$lc1995, kara_class_df$lc2022)
print(transition_matrix_95_04)

melt_transition_matrix <- function(transition_matrix, name) {
  transition_matrix <- as.data.frame(transition_matrix)
  transition_matrix$From <- rownames(transition_matrix)
  transition_matrix <- transition_matrix %>%
    gather(key = "To", value = "Percentage", -From) %>%
    mutate(From = class_names[as.character(From)],
           To = class_names[as.character(To)],
           Transition = name)
  return(transition_matrix)
}

melted_95_04 <- melt_transition_matrix(transition_matrix_95_04, "1995 to 2004")
melted_04_13 <- melt_transition_matrix(transition_matrix_04_13, "2004 to 2013")
melted_13_22 <- melt_transition_matrix(transition_matrix_13_22, "2013 to 2022")

melted_all <- bind_rows(melted_95_04, melted_04_13, melted_13_22)

ggplot(melted_all, aes(x = From, y = To, fill = Percentage)) +
  geom_tile(color = "black", size = 0.5) +  # Add black borders to the tiles
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
  facet_wrap(~ Transition, scales = "free", ncol = 1) +
  labs(title = "Transition Matrices", x = "From Class", y = "To Class", fill = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
