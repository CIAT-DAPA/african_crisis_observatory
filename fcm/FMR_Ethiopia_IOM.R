# **********
#Author:Brenda Chepngetich
# Analyze IOM FMR migration data to visualize trends, analyse co-occurrence of climate-conflict clusters with migration, analyze socioeconomic profiles of migrants.
# This is done using admin level 1 (regions)
# **************
#load packages
library(tidyverse)
library(sf)
library(tmap)
library(readxl)
library(ggplot2)
library(networkD3)
library(htmlwidgets)
library(webshot)
library(tidyr)
library(scatterpie)
library(data.table)

#set working directory
wd <- "C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Brenda/IOM"
setwd(wd)

#read fmr data
file_path <- paste0(wd,"/data/FMR_ethiopia.csv")
fmr <- fread(file_path)


#data cleaning
fmr <- fmr[ ,c("survey_date","total_number_of_persons","departure_country",
               "departure_admin_1","departure_admin_2","departure_admin_3",
               "route","next_destination_country","final_destination_country","destination_country")]
next_count <- sum(!fmr$next_destination_country %in% exclude & !is.na(fmr$next_destination_country)) #327581
final_count <- sum(!fmr$final_destination_country %in% exclude & !is.na(fmr$final_destination_country)) #26047
dest_count <- sum(!fmr$destination_country %in% exclude & !is.na(fmr$destination_country)) #327570
names(fmr)[names(fmr) == "departure_admin_1"] <- "Region"
names(fmr)[names(fmr) == "departure_admin_2"] <- "Zones"
names(fmr)[names(fmr) == "departure_admin_3"] <- "Woredas"
names(fmr)[names(fmr) == "next_destination_country"] <- "Destination"
fmr$survey_date <- as.Date(fmr$survey_date,format = "%m/%d/%Y")
fmr_region <- fmr[!is.na(fmr$Region), ]
fmr_region <- fmr_region[!is.na(fmr_region$Destination), ]
fmr_region <- fmr_region[!fmr_region$Destination == ""]
exclude <- c("Not Specified","","unknown")
fmr_region <- subset(fmr_region, !Region %in% exclude)
fmr_region <- subset(fmr_region, !Destination %in% exclude)
fmr_region <- fmr_region[!(fmr_region$Region == "SNNP" & fmr_region$Zones == "Not Specified"),]
fmr_region <- fmr_region[!(fmr_region$Region == "SNNP" & fmr_region$Zones == ""),]
Central_zones <- c("Hadiya","Guraghe","Siltie","Kembata Tembaro","Halaba","Yem Special")
South_zones <- c("South Omo","Gofa","Gamo","Amaro","Burji","Wolayita","Gedeo","Konso","Derashe","Alle")
fmr_region$Region[fmr_region$Region == "SNNP" & fmr_region$Zones %in% South_zones] <- "South Ethiopia"
fmr_region$Region[fmr_region$Region == "SNNP" & fmr_region$Zones %in% Central_zones] <- "Central Ethiopia"
fmr_region$Region[fmr_region$Region == 'South West Ethiopia Peoples'] <- 'South West Ethiopia'
unique(fmr_region$Region)
fmr_region$Region[fmr_region$Region == 'Gambella'] <- 'Gambela'
Regions <- c("Oromia","Dire Dawa","Tigray","Gambela","Addis Ababa","Afar","South Ethiopia",
             "Central Ethiopia","Somali","Amhara","South West Ethiopia","Sidama","Beneshangul Gumuz",
             "Hareri")
fmr_region <- subset(fmr_region, Region %in% Regions)
fmr_region$Region[fmr_region$Region == 'Hareri'] <- 'Harari'
fmr_region$Region[fmr_region$Region == 'Beneshangul Gumuz'] <- 'Benishangul-Gumuz' 
fmr_region$Destination[fmr_region$Destination == 'ETHIOPIA'] <- 'Ethiopia'
#define HOA countries
HOA <- list("Somalia","Djibouti","Kenya","Uganda","Ethiopia","South Sudan","United Republic of Tanzania", "Tanzania", "Rwanda","Eritrea","Burundi","Ilemi triangle")
unique(fmr_region$Destination)
fmr_region$HOA <- ifelse(fmr_region$Destination %in% HOA, "within", "outside")
fmr_region$year <- format(fmr_region$survey_date, "%Y")
fmr_within <- subset(fmr_region, HOA == "within")
fmr_within_year <- aggregate(total_number_of_persons ~ Region + year, data = fmr_within, FUN = sum)
fmr_outside <- subset(fmr_region, HOA == "outside")
fmr_outside_year <- aggregate(total_number_of_persons ~ Region + year, data = fmr_outside, FUN = sum)

#line plots
trends <- ggplot(fmr_within_year, aes(x = year, y = total_number_of_persons,color = Region, group = Region)) +
  geom_line(linewidth = 0.8) +
  geom_smooth(se = FALSE, method = "loess", linetype = "dashed", size = 0.7) + #Adds a smoothed line (trend line) to the plot to show the overall trend in the data.
  geom_point(size = 0.8) +
  labs(
    title = "Migration Trends from Ethiopia to countries within HOA (2018-2024)",
    x = "Time",
    y = "Number of Migrants"
  ) +
  theme_bw()+
  theme(legend.text = element_text(size=9),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

trends
#save the plot
file <- paste0(wd,"/FMR/trends_within.png")
ggsave(file, plot = trends, width = 11, height = 8, dpi = 600)

#facets
shp <- st_read(paste0(wd, "/data/eth_admbnda_admins_csa_bofedb_2024.shp/eth_admbnda_adm1_csa_bofedb_2024.shp"))
shp <- shp[ ,"admin1Name"]
names(shp)[names(shp) == "admin1Name"] <- "Region"
shp_outside <- sf::st_as_sf(merge(fmr_outside_year, shp, by="Region"))
shp_within <- sf::st_as_sf(merge(fmr_within_year, shp, by="Region"))

map <- tm_shape(shp_within) +
  tm_polygons(col ="total_number_of_persons", title = "Migrants", style = "cont", palette="YlOrBr", legend.show = T, border.col = "black", lwd = 1) +
  tm_text("Region", size =0.42, col="black", auto.placement = TRUE)+
  tm_facets(by= "year", ncol = 3) +
  tm_compass(type = "8star", size = 2,position = c("right", "top")) +
  tm_scale_bar(breaks = c(0, 50, 100), text.size = 1, 
               position = c("right", "bottom"))+
  tm_layout(main.title = "Ethiopia migration within HOA",
            legend.outside = T,
            asp = 1.4,
            # inner.margins = c(0,0,0,0)
  )
map
tmap_save(map,  dpi= 600,  height=8.3, width=11.7, units="in",scale = 1.6,
          filename=paste0(wd,"/FMR/within_facets.png"))
#sankey plots
within_sankey <- fmr_within
names(within_sankey)[names(within_sankey) == "Region"] <- "Origin"
within_sankey <- aggregate(total_number_of_persons ~ Origin + Destination, data=within_sankey, FUN=sum)
within_nodes <- data.frame(name = unique(c(as.character(within_sankey$Origin), 
                                    as.character(within_sankey$Destination))))
View(within_nodes)

# Generate links dataframe using indices
within_links <- data.frame(
  source = match(within_sankey$Origin, within_nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  target = match(within_sankey$Destination, within_nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  value = within_sankey$total_number_of_persons
)
within_sankey_ <- sankeyNetwork(
  Links = within_links,
  Nodes = within_nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",  # optional
  fontSize = 16,  # optional
  nodeWidth = 30  # optional
)
file_path <- paste0(wd,"/FMR/outside_sankey_plot.html")
saveWidget(outside_sankey_, file_path)
webshot::install_phantomjs()
webshot(file_path, paste0(wd,"/FMR/outside_sankey_plot.png"), zoom = 600/72)

outside_sankey <- fmr_outside
names(outside_sankey)[names(outside_sankey) == "Region"] <- "Origin"
outside_sankey <- aggregate(total_number_of_persons ~ Origin + Destination, data=outside_sankey, FUN=sum)
outside_sankey$region <- ""
unique(outside_sankey$Destination)
outside_sankey$region <- ifelse(outside_sankey$Destination %in% MENA, "MENA", 
                                ifelse(outside_sankey$Destination %in% Sub_saharan, "Sub Saharan Africa",
                                       ifelse(outside_sankey$Destination %in% Europe, "Europe",
                                              ifelse(outside_sankey$Destination %in% Asia, "Asia",
                                                     ifelse(outside_sankey$Destination %in% America, "America",
                                                            ifelse(outside_sankey$Destination %in% pacific, "Pacific",
                                                                   ifelse(outside_sankey$Destination %in% Caribbean, "Caribbean",
                                                                          "Unknown")))))))
outside_sankey <- aggregate(total_number_of_persons ~ Origin + region, data=outside_sankey, FUN=sum)
outside_nodes <- data.frame(name = unique(c(as.character(outside_sankey$Origin), 
                                           as.character(outside_sankey$region))))
outside_links <- data.frame(
  source = match(outside_sankey$Origin, outside_nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  target = match(outside_sankey$region, outside_nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  value = outside_sankey$total_number_of_persons
)
outside_sankey_ <- sankeyNetwork(
  Links = outside_links,
  Nodes = outside_nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",  # optional
  fontSize = 16,  # optional
  nodeWidth = 30  # optional
)

#overlay climate_conflict hotspots with migration
#get total migrants over the years
fmr_within_migrants <- merge(shp,aggregate(total_number_of_persons ~ Region, data = fmr_within_year, FUN = sum))
fmr_outside_migrants <- merge(shp,aggregate(total_number_of_persons ~ Region, data = fmr_outside_year, FUN = sum))
#load climate-conflict hotspots
clim_conflict <- sf::st_read(paste0(wd, "/data/clim_conflict_ips_overlays.geojson"))

#rename columns
x <- clim_conflict$intersect_conf_clim
unique(clim_conflict$intersect_conf_clim)
clim_conflict$intersect_conf_clim[x == "High conflict-[Low levels of precipitation/High levels of Heat stress]"] <- "High conflict + High drought stress"
clim_conflict$intersect_conf_clim[x == "Limited conflict-[Low levels of precipitation/High levels of Heat stress]"] <- "Limited conflict + High drought stress"
clim_conflict$intersect_conf_clim[x == "Moderate conflict-[Low levels of precipitation/High levels of Heat stress]"] <- "Moderate conflict + High drought stress"
clim_conflict$intersect_conf_clim[x == "High conflict-[High levels of precipitation/Low levels of Heat stress]"] <- "High conflict + Low drought stress"
clim_conflict$intersect_conf_clim[x == "High conflict-[Moderate levels of precipitation/moderate levels of Heat stress]"] <- "High conflict + Moderate drought stress"

#extract the needed data
required <- c("High conflict + High drought stress","Limited conflict + High drought stress"
              ,"Moderate conflict + High drought stress","High conflict + Low drought stress"
              ,"High conflict + Moderate drought stress")
clusters <- subset(clim_conflict, intersect_conf_clim %in% required)
c <- clusters$intersect_conf_clim
clusters$clust[c=="High conflict + High drought stress"] <- 1
clusters$clust[c=="High conflict + Moderate drought stress"] <- 2
clusters$clust[c=="High conflict + Low drought stress"] <- 3
clusters$clust[c=="Limited conflict + High drought stress"] <- 4
clusters$clust[c=="Moderate conflict + High drought stress"] <- 5
clusters$clust <- as.factor(clusters$clust)

label <- c("High conflict + High drought", "High conflict + Moderate drought",
           "High conflict + Low drought", "Limited conflict + High drought",
           "Moderate conflict + High drought")
#mapping
tmap_mode("plot")
map <- tm_shape(fmr_within_migrants) +
  tm_polygons(col ="total_number_of_persons", title = "Migrants", style = "cont", palette="Blues", 
              legend.show = T, border.col = "black", lwd = 1) +
  tm_shape(clusters) +
  tm_fill(col= "clust", palette="-YlOrRd", title="Conflict-Climate Intersection",
          legend.show = T, labels= label)+
  tm_shape(fmr_outside_migrants) +
  tm_text("Region", size=0.6,col="black")+
  tm_compass(type = "8star", size=4,position = c("right", "bottom")) +
  tm_scale_bar(breaks = c(0, 50, 100), text.size = 1.5, 
               position = c("right", "bottom"))+
  tm_layout(legend.outside=F,
            main.title = "Climate-Conflict Migration co-occurence within HOA",
            main.title.position = "center",
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
tmap_save(map,  dpi= 600,  height=8.3, width=11.7, units="in",scale = 1.6,
          filename=paste0(wd,"/FMR/within_Overlay.png"))



