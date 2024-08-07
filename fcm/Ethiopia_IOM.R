#Author:Brenda Chepngetich
#install and load required packages
install.packages("networkD3")
install.packages("htmlwidgets")
install.packages("webshot")
library(geodata)
library(tidyverse)
library(sf)
library(tmap)
library(readxl)
library(ggplot2)
library(networkD3)
library(htmlwidgets)
library(webshot)

#set working directory
wd <- "D:/OneDrive - CGIAR/SA_Team/Brenda/IOM"
setwd(wd)

#read fms data
file <- "D:/OneDrive - CGIAR/IOM - CGIAR Climate Security Coordination/Data/OneDrive_1_7-23-2024/FMS.xlsx"
fms <- read_excel(file)
View(fms)
#data cleaning
fms <- subset(fms, departure_country == "Ethiopia")
names(fms)[names(fms) == "departure_admin1"] <- "Region"
names(fms)[names(fms) == "final_destination_country"] <- "Destination"
fms$survey_date <- as.Date(fms$survey_date)
exclude <- c("Not Specified", "Not specified", "NA", "Other", "Unknown","I don’t know","transit","Other Country","0")
fms <- subset(fms,!(Destination %in% exclude))
fms <- subset(fms, Region != "Not Specified")
#correct region names
fms$Region[fms$Region == 'Hareri'] <- 'Harari'
fms$Region[fms$Region == 'Beneshangul Gumuz'] <- 'Benishangul-Gumuz' 
fms$Region[fms$Region == 'Gambella'] <- 'Gambela'
#remove na values
fms_ <- fms[!is.na(fms$Destination), ]
View(fms_merged_regions)
#merged regions
fms_merged_regions <- fms_
fms_merged_regions$Region[fms_merged_regions$Region == "South Ethiopia"] <- "SNNP"
fms_merged_regions$Region[fms_merged_regions$Region == "Central Ethiopia"] <- "SNNP"
fms_merged_regions$Region[fms_merged_regions$Region == "South West Ethiopia Peoples"] <- "South West Ethiopia"

#within HOA route
HOA <- list("Somalia","Djibouti","Kenya","Uganda","Ethiopia","South Sudan","United Republic of Tanzania", "Tanzania", "Rwanda","Eritrea","Burundi")
fms_merged_regions$HOA <- ifelse(fms_merged_regions$Destination %in% HOA, "within", "outside")

#visualize migration trends per region
fms_region <- fms_merged_regions[, c("survey_date","Region","Destination","HOA")]
fms_within <- subset(fms_region, HOA == "within")
fms_outside <- subset(fms_region, HOA == "outside")
#fms within
fms_within$migrants <- 1
fms_within <- aggregate(migrants ~ survey_date + Region, data = fms_within, FUN = sum)
fms_within$year <- format(fms_within$survey_date, "%Y")
within_by_year <- aggregate(migrants ~ Region + year, data = fms_within, FUN = sum)
#fms outside
fms_outside$migrants <- 1
fms_outside <- aggregate(migrants ~ survey_date + Region, data = fms_outside, FUN = sum)
fms_outside$year <- format(fms_outside$survey_date, "%Y")
outside_by_year <- aggregate(migrants ~ Region + year, data = fms_outside, FUN = sum)

#plotting line plot
region_trends <- ggplot(within_by_year, aes(x = year, y = migrants,color = Region, group = Region)) +
  geom_line(size = 0.8) +
  geom_smooth(se = TRUE, method = "loess", linetype = "dashed", size = 0.7) +
  geom_point(size = 0.8) +
  labs(
    title = "Migration Trends from Ethiopia to countries within HOA (2018-2024)",
    x = "Time",
    y = "Number of Migrants"
  ) +
  theme_bw()
region_trends
#save the plot
file <- paste0(wd,"/Results/within/trends_by_year.png")
ggsave(file, plot = region_trends, width = 11, height = 8, dpi = 300)

#map plotting
# ETH_adm1 <- st_as_sf(geodata::gadm(country = 'ETH', level = 1, version = 'latest' , path = tempdir()))
# ETH_adm1 <- ETH_adm1[,"NAME_1"]
shp <- st_read(paste0(wd, "/data/admin areas OCHA/eth_admbnda_adm1_csa_bofedb_2021.shp"))
shp <- shp[ ,"ADM1_EN"]
shp$Region[shp$Region == "Gambella"] <- "Gambela"
shp$Region[shp$Region == "Benishangul Gumz"] <- "Benishangul-Gumuz"
names(shp)[names(shp) == "ADM1_EN"] <- "Region"
plot(shp)
#to do: add admin 2 data

# #rename regions to match those in the dataframe
# ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Addis Abeba'] <- 'Addis Ababa'
# ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Benshangul-Gumaz'] <- 'Beneshangul Gumuz'
# ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Gambela Peoples'] <- 'Gambella'
# ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Southern Nations, Nationalities'] <- 'SNNP'
# ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Harari People'] <- 'Harari'
# #rename
# names(ETH_adm1)[names(ETH_adm1) == "NAME_1"] <- "Region"

#merge shp with data
# ETH <- merge(fms_by_year,ETH_adm1,by="Region")
# ETH <- sf::st_as_sf(ETH)
Eth_outside <- sf::st_as_sf(merge(outside_by_year, shp, by="Region"))
Eth_within <- sf::st_as_sf(merge(within_by_year, shp, by="Region"))
#map plotting to create facets
map <- tm_shape(Eth_within) +
  tm_polygons(col ="migrants", title = "Migrants", style = "cont", palette="YlOrBr", legend.show = T, border.col = "black", lwd = 1) +
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
tmap_save(map,  dpi= 300,  height=8.3, width=11.7, units="in",scale = 1.6,
          filename=paste0(wd,"/Results/within/facets.png"))


#plot sankey diagram
within_sankey <- subset(fms_region, HOA == "within")
names(within_sankey)[names(within_sankey) == "Region"] <- "Origin"
within_sankey$migrants <- 1
within_sankey <- aggregate(migrants ~ Origin + Destination, data=within_sankey, FUN=sum)

outside_sankey <- subset(fms_region, HOA == "outside")
names(outside_sankey)[names(outside_sankey) == "Region"] <- "Origin"
outside_sankey$migrants <- 1
outside_sankey <- aggregate(migrants ~ Origin + Destination, data=outside_sankey, FUN=sum)
#define regions outside HOA
Northern_Africa <- c("Algeria", "Egypt", "Libya", "Morocco", "Sudan", "Tunisia", "Western Sahara")
Western_Africa <- c("Benin", "Burkina Faso", "Cape Verde", "Côte d'Ivoire", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Liberia", "Mali", "Mauritania", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo")
Central_Africa <- c("Angola", "Cameroon", "Central African Republic", "Chad", "Congo", "Democratic Republic of the Congo", "Equatorial Guinea", "Gabon", "São Tomé and Príncipe")
Southern_Africa <- c("Botswana", "Eswatini", "Lesotho", "Namibia", "South Africa")
Eastern_Africa <- c("Burundi", "Comoros", "Djibouti", "Eritrea", "Ethiopia", "Kenya", "Madagascar", "Malawi", "Mauritius", "Mozambique", "Rwanda", "Seychelles", "Somalia", "South Sudan", "Tanzania", "Uganda", "Zambia", "Zimbabwe")
outside_sankey$region <- ifelse(outside_sankey$Destination %in% Northern_Africa, "Northern Africa", 
                         ifelse(outside_sankey$Destination %in% Western_Africa, "Western Africa",
                         ifelse(outside_sankey$Destination %in% Southern_Africa, "Southern Africa",
                         ifelse(outside_sankey$Destination %in% Eastern_Africa, "Eastern Africa",
                         ifelse(outside_sankey$Destination %in% Central_Africa, "Central Africa",
                                "")))))

# Generate unique nodes
nodes <- data.frame(name = unique(c(as.character(outside_sankey$Origin), 
                                    as.character(outside_sankey$Destination))))
View(links)

# Generate links dataframe using indices
links <- data.frame(
  source = match(outside_sankey$Origin, nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  target = match(outside_sankey$Destination, nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  value = outside_sankey$migrants
)
sankey_ <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",  # optional
  fontSize = 12,  # optional
  nodeWidth = 30  # optional
)
saveWidget(sankey_, paste0(wd,"/Results/outside/sankey_plot.html"))
#initialize webshot
webshot::install_phantomjs()

# Save as PNG
webshot(paste0(wd,"/Results/within/sankey_plot.html"), paste0(wd,"/Results/within/sankey_plot.png"))

#overlay climate_conflict hotspots with migration
#get total migrants over the years
total_migrants <- aggregate(migrants ~ Region, data = fms_by_year, FUN = sum)
ETH_total <- merge(ETH_adm1, total_migrants, by = "Region") 
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
View(clusters)
unique(clusters$clust)
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
map <- tm_shape(ETH_total) +
  tm_polygons(col ="migrants", title = "Migrants", style = "cont", palette="YlOrBr", 
              legend.show = T, border.col = "black", lwd = 1) +
  tm_shape(clusters) +
  tm_fill(col= "clust", palette="-YlOrRd", title="Conflict-Climate Intersection",
          legend.show = T, labels= label)+
  tm_shape(ETH_total) +
  tm_text("Region", size=0.6,col="black")+
  tm_compass(type = "8star", size=4,position = c("right", "bottom")) +
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
tmap_save(map,  dpi= 300,  height=8.3, width=11.7, units="in",scale = 1.6,
          filename=paste0(wd,"/Overlay2.png"))
