library(sf)
library(tmap)
library(scatterpie)
library(dplyr)
library(ggplot2)
rm(list=ls(all=TRUE))

reasons <- read.csv("C:/Users/bchepngetich/OneDrive - CGIAR/IOM - CGIAR Climate Security Coordination/Data/OneDrive_1_7-23-2024/FMS_migration_reasons.csv")
exclude <- c("","Not Specified")
reasons <- reasons[!reasons$region_of_departure %in% exclude, ]
shp <- st_read(paste0("C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Brenda/IOM/data/eth_admbnda_admins_csa_bofedb_2024.shp/eth_admbnda_adm1_csa_bofedb_2024.shp"))
shp <- shp[ ,"admin1Name"]
names(shp)[names(shp) == "admin1Name"] <- "region_of_departure"
#rename regions
reasons$region_of_departure[reasons$region_of_departure == "Hareri"] <- "Harari"
reasons$region_of_departure[reasons$region_of_departure == "Gambella"] <- "Gambela"
reasons$region_of_departure[reasons$region_of_departure == "South West Ethiopia Peoples"] <- "South West Ethiopia"
reasons$region_of_departure[reasons$region_of_departure == "Beneshangul Gumuz"] <- "Benishangul-Gumuz"
##
snnp <- reasons[reasons$region_of_departure == "SNNP", ]
Central_zones <- c("Hadiya","Guraghe","Siltie","Kembata Tembaro","Halaba")
South_zones <- c("South Omo","Gofa","Gamo","Amaro","Burji","Wolayita","Gedeo","Konso","Derashe")
reasons$region_of_departure[reasons$region_of_departure == "SNNP" & reasons$departure_admin2 %in% South_zones] <- "South Ethiopia"
reasons$region_of_departure[reasons$region_of_departure == "SNNP" & reasons$departure_admin2 %in% Central_zones] <- "Central Ethiopia"
reasons$region_of_departure[reasons$departure_admin2 == 'Dawuro'] <- 'South West Ethiopia'
reasons$region_of_departure[reasons$departure_admin2 == 'Konta Special'] <- 'South West Ethiopia'
reasons$region_of_departure[reasons$departure_admin2 == 'Alle'] <- 'South Ethiopia'
reasons$region_of_departure[reasons$departure_admin2 == 'Sidama'] <- 'Sidama'
reasons <- reasons[!reasons$region_of_departure == "SNNP",]
##

push <- reasons[, c("region_of_departure","push_reason_economic","push_reason_education","push_reason_family","push_reason_acc_serv","push_reason_natural_disaster",
                      "push_reason_environmental","push_reason_conf","push_reason_violence","push_reason_other","push_reason_agropast","destination_EHoA")]
push_within <- push[push$destination_EHoA == 1,]
push_outside <- push[push$destination_EHoA == 0,]
push_within_ <- push_within %>% group_by(region_of_departure) %>% summarize(
  Agro_pastoral = sum(push_reason_agropast),
  Conflict = sum(push_reason_conf, na.rm = TRUE),
  Violence = sum(push_reason_violence),
  Environmental = sum(push_reason_environmental, na.rm = TRUE),
  Natural_disaster = sum(push_reason_natural_disaster),
  Economic = sum(push_reason_economic),
  Education = sum(push_reason_education, na.rm = TRUE),
  Family = sum(push_reason_family),
  Service_access = sum(push_reason_acc_serv),
  Other = sum(push_reason_other)
)
within_pie <- merge(shp,push_within_, by = "region_of_departure")
within_centroids <- st_centroid(within_pie)
within_centroids <- cbind(within_centroids, st_coordinates(within_centroids))
within_centroids_df <- as.data.frame(within_centroids)

push_outside_ <- push_outside %>% group_by(region_of_departure) %>% summarize(
  Agro_pastoral = sum(push_reason_agropast),
  Conflict = sum(push_reason_conf, na.rm = TRUE),
  Violence = sum(push_reason_violence),
  Environmental = sum(push_reason_environmental, na.rm = TRUE),
  Natural_disaster = sum(push_reason_natural_disaster),
  Economic = sum(push_reason_economic),
  Education = sum(push_reason_education, na.rm = TRUE),
  Family = sum(push_reason_family),
  Service_access = sum(push_reason_acc_serv),
  Other = sum(push_reason_other)
)
push_outside_ <- push_outside_[!push_outside_$region_of_departure == "SNNP",] #to remove later
outside_pie <- merge(shp,push_outside_, by = "region_of_departure")
outside_centroids <- st_centroid(outside_pie)
outside_centroids <- cbind(outside_centroids, st_coordinates(outside_centroids))
outside_centroids_df <- as.data.frame(outside_centroids)


#plot
plot <- ggplot() +
  geom_sf(data = within_pie, fill = "gray80", color = "black")+
  geom_scatterpie(data = within_centroids_df, aes(x = X, y = Y, r = 0.5, group=region_of_departure),  # Add pie charts using centroid coordinates
                  cols = c("Agro_pastoral","Conflict","Violence","Environmental", "Natural_disaster","Economic", "Education", "Family","Service_access","Other"),
                  color = "black",
                  size=0.3)+
  geom_text(data=within_centroids,
            aes(x=X + 0.7, y=Y + 0.4,label=region_of_departure),
            size=3,
            color="black",
            fontface="bold"
  )+
  coord_sf()+
  theme_minimal() +
  # theme_classic()+
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_blank(),        # Remove x-axis text (labels)
    axis.text.y = element_blank(),        # Remove y-axis text (labels)
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position.inside = c(0.95, 0.05),  # Position the legend in the bottom right
    legend.justification = c("right", "bottom"),
    panel.border = element_rect(color = "gray50", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black")
    
  )+
  scale_fill_manual(values = c("Agro_pastoral" = "orange","Conflict"="red4","Violence"="red3","Environmental"="green", 
                               "Natural_disaster"="red1", "Economic" = "lightgreen", "Education" = "white",
                               "Family" = "purple","Service_access" = "thistle", "Other" = "lavender"))+  # Custom colors
  labs(title = "",fill="Push Reason")
plot
output <- "C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Brenda/IOM/Final results/push_within_pie.png"
ggsave(output, plot = plot, dpi= 600,  height=8.3, width=11.7, units="in")

#pull agro pastoralists
pull <- reasons[,c("region_of_departure","pull_reason_agropastoral","destination_EHoA")]
pull_within <- pull[pull$destination_EHoA == 1,]
pull_outside <- pull[pull$destination_EHoA == 0,]
pull_within_ <- pull_within %>% group_by(region_of_departure) %>% summarize(Agro_pastoral = sum(pull_reason_agropastoral, na.rm=TRUE))
pull_outside_ <- pull_outside %>% group_by(region_of_departure) %>% summarize(Agro_pastoral = sum(pull_reason_agropastoral, na.rm=TRUE))
agropastoral_within <- merge(shp, pull_within_, by="region_of_departure")
agropastoral_outside <- merge(shp, pull_outside_, by="region_of_departure")

#map
map <- tm_shape(agropastoral_within)+
  tm_polygons(col ="Agro_pastoral", title = "Pull Agro-pastoral", style = "cont", palette="YlOrBr", legend.show = T, border.col = "black", lwd = 1) +
  tm_text("region_of_departure", size =0.6, col="black", auto.placement = TRUE)+
  tm_compass(type = "8star", size=2,position = c("right", "bottom")) +
  tm_scale_bar(breaks = c(0, 50, 100), text.size = 1.5, 
               position = c("right", "bottom"))+
  tm_layout(legend.outside=F,
            legend.text.size = 0.6,
            legend.text.color = "black",
            legend.title.size= 0.9,
            legend.title.color = "black",
            legend.title.fontface = 2,
            legend.frame=F,
            legend.width=0.7,
            asp=1.4,
            inner.margins = c(0,0.05,0,0.05)
 
  )
map
tmap_save(map,  dpi= 600,  height=8.3, width=11.7, units="in",scale=1.6,
          filename="C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Brenda/IOM/Final results/agropastoral_within.png")
