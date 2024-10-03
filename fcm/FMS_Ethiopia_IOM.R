# ****************************
#Author:Brenda Chepngetich
# Analyze IOM FMS migration data to visualize trends using line plots and facet maps
# Analyse co-occurrence of climate-conflict clusters with migration
# Analyze socioeconomic profiles of migrants
# Use admin level 1 (regions)
# *******************************************
#install and load required packages
install.packages("networkD3")
install.packages("htmlwidgets")
install.packages("webshot")
install.packages("scatterpie")
install.packages("ggplot2")

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

#set working directory
wd <- "D:/OneDrive - CGIAR/SA_Team/Brenda/IOM"
setwd(wd)

#read fms data
file <- "D:/OneDrive - CGIAR/IOM - CGIAR Climate Security Coordination/Data/OneDrive_1_7-23-2024/FMS.xlsx"
fms <- read_excel(file)

#data cleaning
fms <- subset(fms, departure_country == "Ethiopia")
#TOTALS VALUES = 110,066
adm1_count <- sum(!fms$departure_admin1 %in% exclude & !is.na(fms$departure_admin1)) #107,688
adm2_count <- sum(!fms$departure_admin2 %in% exclude & !is.na(fms$departure_admin2)) #92,682
adm3_count <- sum(!fms$departure_admin3 %in% exclude & !is.na(fms$departure_admin3)) #80,488
names(fms)[names(fms) == "departure_admin1"] <- "Region"
names(fms)[names(fms) == "departure_admin2"] <- "Zones"
names(fms)[names(fms) == "departure_admin3"] <- "Woredas"
names(fms)[names(fms) == "final_destination_country"] <- "Destination"
fms$survey_date <- as.Date(fms$survey_date)
exclude <- c("Not Specified", "Not specified", "NA", "Other", "Unknown","I don’t know","transit","Other Country","0")
fms <- subset(fms,!(Destination %in% exclude))
fms <- subset(fms, Region != "Not Specified")
#correct region names
fms$Region[fms$Region == 'Hareri'] <- 'Harari'
fms$Region[fms$Region == 'Beneshangul Gumuz'] <- 'Benishangul-Gumuz' 
fms$Region[fms$Region == 'Gambella'] <- 'Gambela'
snnp <- fms[fms$Region == "SNNP", ]
Central_zones <- c("Hadiya","Guraghe","Siltie","Kembata Tembaro","Halaba")
South_zones <- c("South Omo","Gofa","Gamo","Amaro","Burji","Wolayita","Gedeo","Konso","Derashe")
fms$Region[fms$Region == "SNNP" & fms$Zones %in% South_zones] <- "South Ethiopia"
fms$Region[fms$Region == "SNNP" & fms$Zones %in% Central_zones] <- "Central Ethiopia"
fms$Region[fms$Zones == 'Dawuro'] <- 'South West Ethiopia'
fms$Region[fms$Zones == 'Konta Special'] <- 'South West Ethiopia'
fms$Region[fms$Zones == 'Sidama'] <- 'Sidama'
fms$Region[fms$Zones == 'Alle'] <- 'South Ethiopia'
fms <- subset(fms, !fms$Region == "SNNP")
fms$Region[fms$Region == 'South West Ethiopia Peoples'] <- 'South West Ethiopia'
#remove na values
fms <- fms[!is.na(fms$Destination), ]
fms <- fms[, c("survey_date","Region","Zones","Woredas","Destination","current_employment_status",
               "highest_education_level","recent_employment_status_before_journey","main_profession_of_recent_or_current_job","other_profession_of_recent_or_current_job")]
fms$Destination[fms$Destination == "Tanzania"] <- "United Republic of Tanzania"

#Define HOA route
HOA <- list("Somalia","Djibouti","Kenya","Uganda","Ethiopia","South Sudan","United Republic of Tanzania", "Tanzania", "Rwanda","Eritrea","Burundi","Ilemi triangle")

#visualize migration trends per region
fms_region <- fms 
fms_region$HOA <- ifelse(fms_region$Destination %in% HOA, "within", "outside")
fms_region$migrants <- 1
fms_region$year <- format(fms_region$survey_date, "%Y")
fms_within <- subset(fms_region, HOA == "within")
fms_outside <- subset(fms_region, HOA == "outside")

#fms within
# within_by_date <- aggregate(migrants ~ survey_date + Region, data = fms_within, FUN = sum)
within_by_year <- aggregate(migrants ~ Region + year, data = fms_within, FUN = sum)
#fms outside
# outside_by_date <- aggregate(migrants ~ survey_date + Region, data = fms_outside, FUN = sum)
outside_by_year <- aggregate(migrants ~ Region + year, data = fms_outside, FUN = sum)
# here*********
#plotting line plot
region_trends <- ggplot(within_by_year, aes(x = year, y = migrants,color = Region, group = Region)) +
  geom_line(linewidth = 0.8) +
  geom_smooth(se = FALSE, method = "loess", linetype = "dashed", size = 0.7) + #Adds a smoothed line (trend line) to the plot to show the overall trend in the data.
  geom_point(size = 0.8) +
  labs(
    title = "",
    x = "Time",
    y = "Number of Migrants"
  ) +
  theme_bw()+
  theme(legend.text = element_text(size=9),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
region_trends
#save the plot
file <- paste0(wd,"/Final results/trends_within.png")
ggsave(file, plot = region_trends, width = 11, height = 8, dpi = 600)

#map plotting
admin1 <- st_read(paste0(wd, "/data/eth_admbnda_admins_csa_bofedb_2024.shp/eth_admbnda_adm1_csa_bofedb_2024.shp"))
admin1 <- admin1[ ,"admin1Name"]
names(admin1)[names(admin1) == "admin1Name"] <- "Region"

#merge shp with data
Eth_outside <- sf::st_as_sf(merge(outside_by_year, admin1, by="Region"))
Eth_within <- sf::st_as_sf(merge(within_by_year, admin1, by="Region"))

#map plotting to create facets
map <- tm_shape(Eth_within) +
  tm_polygons(col ="migrants", title = "Migrants", style = "cont", palette="YlOrBr", legend.show = T, border.col = "black", lwd = 1) +
  tm_text("Region", size =0.42, col="black", auto.placement = TRUE)+
  tm_facets(by= "year", ncol = 3) +
  tm_compass(type = "8star", size = 2,position = c("right", "top")) +
  tm_scale_bar(breaks = c(0, 50, 100), text.size = 1, 
               position = c("right", "bottom"))+
  tm_layout(legend.outside = T,
            asp = 1.4,
            # inner.margins = c(0,0,0,0)
            )
map
tmap_save(map,  dpi= 600,  height=8.3, width=11.7, units="in",scale = 1.6,
          filename=paste0(wd,"/Final results/within_facets.png"))


#plot sankey diagram
within_sankey <- fms_within
names(within_sankey)[names(within_sankey) == "Region"] <- "Origin"
within_sankey <- aggregate(migrants ~ Origin + Destination, data=within_sankey, FUN=sum)
outside_sankey <- fms_outside
names(outside_sankey)[names(outside_sankey) == "Region"] <- "Origin"
outside_sankey <- aggregate(migrants ~ Origin + Destination, data=outside_sankey, FUN=sum)
outside_sankey$Destination[outside_sankey$Destination == "Saudi arabie"] <- "Kingdom of Saudi Arabia"
outside_sankey$Destination[outside_sankey$Destination == "SWZ"] <- "Swaziland"
outside_sankey$Destination[outside_sankey$Destination == "sweden"] <- "Sweden"
#define regions outside HOA
MENA <- list("Algeria","Morocco","Libya","Tunisia","Egypt","Jordan","Syria", "Saudi Arabia",
             "Iraq","United Arab Emirates","Iran","Yemen","Oman","Israel","Kingdom of Saudi Arabia","Gaza Strip","Hala'ib triangle","Ma'tan al-Sarra")
Northern_Africa <- c("Algeria", "Egypt", "Libya", "Morocco", "Sudan", "Tunisia", "Western Sahara")
Western_Africa <- c("Benin", "Burkina Faso", "Cape Verde", "Côte d'Ivoire", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Liberia", "Mali", "Mauritania", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo")
Central_Africa <- c("Angola", "Cameroon", "Central African Republic", "Chad", "Congo", "Democratic Republic of the Congo", "Equatorial Guinea", "Gabon", "São Tomé and Príncipe")
Southern_Africa <- c("Botswana", "Eswatini", "Lesotho", "Namibia", "South Africa","Swaziland")
Eastern_Africa <- c("Burundi", "Comoros", "Djibouti", "Eritrea", "Ethiopia", "Kenya", "Madagascar", "Malawi", "Mauritius", "Mozambique", "Rwanda", "Seychelles", "Somalia", "South Sudan", "Tanzania", "Uganda", "Zambia", "Zimbabwe")
Sub_saharan <- c(Western_Africa,Central_Africa,Southern_Africa, Eastern_Africa, "Sudan","Western Sahara")
Eastern_Europe = c("Belarus", "Bulgaria", "Czech Republic", "Hungary", "Moldova", "Poland", "Romania", "Russia", "Slovakia", "Ukraine")
Northern_Europe = c("Denmark", "Estonia", "Finland", "Iceland", "Ireland", "Latvia", "Lithuania", "Norway", "Bouvet Island", "Sweden", "United Kingdom", "London")
Southern_Europe = c("Albania", "Andorra", "Bosnia and Herzegovina", "Croatia", "Greece", "Italy", "Kosovo", "Malta", "Montenegro", "North Macedonia", "Portugal", "San Marino", "Serbia", "Slovenia", "Spain", "Vatican City")
Western_Europe = c("Austria", "Belgium", "France", "French Southern and Antarctic Territories", "Germany", "Liechtenstein", "Luxembourg", "Monaco", "Netherlands", "Switzerland")
Europe <- c(Eastern_Europe,Northern_Europe,Southern_Europe,Western_Europe,"UNK","Cayman Islands","Europe","U.K. of Great Britain and Northern Ireland")
Central_Asia = c("Kazakhstan", "Kyrgyzstan", "Tajikistan", "Turkmenistan", "Uzbekistan")
East_Asia = c("China", "Japan", "Mongolia", "North Korea", "South Korea", "Taiwan", "Republic of Korea")
South_Asia = c("Afghanistan", "Bangladesh", "Bhutan", "India", "Maldives", "Nepal", "Pakistan", "Sri Lanka")
Southeast_Asia = c("Brunei", "Cambodia", "East Timor", "Indonesia", "Laos", "Malaysia", "Myanmar", "Philippines", "Singapore", "Thailand", "Vietnam")
Western_Asia = c("Armenia", "Azerbaijan", "Bahrain", "Cyprus", "Georgia", "Kuwait", "Lebanon","Qatar", "Turkey")
Asia <- c(Central_Asia,East_Asia,South_Asia,Southeast_Asia,Western_Asia,"Aksai Chin","Arunachal Pradesh","China/India","Jammu and Kashmir","Kuril islands","Paracel Islands","Scarborough Reef","Senkaku Islands")
North_America = c("Canada", "Mexico", "United States","United States of America")
Central_America = c("Belize", "Costa Rica", "El Salvador", "Guatemala", "Honduras", "Nicaragua", "Panama")
South_America = c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Ecuador", "Guyana", "Paraguay", "Peru", "Suriname", "Uruguay", "Venezuela")
America <- c(North_America,Central_America,South_America,"American Samoa")
Australia_and_New_Zealand = c("Australia", "New Zealand")
Melanesia = c("Fiji", "Papua New Guinea", "Solomon Islands", "Vanuatu")
Micronesia = c("Kiribati", "Marshall Islands", "Micronesia (Federated States of)", "Nauru", "Palau")
Polynesia = c("Samoa", "Tonga", "Tuvalu")
pacific <- c(Australia_and_New_Zealand,Melanesia,Micronesia,Polynesia)
Caribbean = c("Antigua and Barbuda", "Bahamas", "Barbados", "Cuba", "Dominica", "Dominican Republic", "Grenada", "Haiti", "Jamaica", "Saint Kitts and Nevis", "Saint Lucia", "Saint Vincent and the Grenadines", "Trinidad and Tobago")


outside_sankey$region <- ifelse(outside_sankey$Destination %in% MENA, "MENA", 
                         ifelse(outside_sankey$Destination %in% Sub_saharan, "Sub Saharan Africa",
                         ifelse(outside_sankey$Destination %in% Europe, "Europe",
                         ifelse(outside_sankey$Destination %in% Asia, "Asia",
                         ifelse(outside_sankey$Destination %in% America, "America",
                         ifelse(outside_sankey$Destination %in% pacific, "Pacific",
                         ifelse(outside_sankey$Destination %in% Caribbean, "Caribbean",
                                "Unknown")))))))




outside_sankey_ <- aggregate(migrants ~ Origin + region, data=outside_sankey, FUN=sum)
# Generate unique nodes
nodes <- data.frame(name = unique(c(as.character(within_sankey$Origin), 
                                    as.character(within_sankey$Destination))))
View(links)

# Generate links dataframe using indices
links <- data.frame(
  source = match(within_sankey$Origin, nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  target = match(within_sankey$Destination, nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  value = within_sankey$migrants
)
sankey_ <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",  # optional
  fontSize = 16,  # optional
  nodeWidth = 30  # optional
)
file_path <- paste0(wd,"/Final results/within_sankey_plot.html")
saveWidget(sankey_, file_path)
#initialize webshot
webshot::install_phantomjs()

# Save as PNG
webshot(file_path, paste0(wd,"/Final results/within_sankey_plot.png"), zoom = 600/72)

#overlay climate_conflict hotspots with migration
#get total migrants over the years
within_migrants <- merge(admin1,aggregate(migrants ~ Region, data = within_by_year, FUN = sum))
outside_migrants <- merge(admin1,aggregate(migrants ~ Region, data = outside_by_year, FUN = sum))
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
map <- tm_shape(within_migrants) +
  tm_polygons(col ="migrants", title = "Migrants", style = "cont", palette="Blues", 
              legend.show = T, border.col = "black", lwd = 1) +
  tm_shape(clusters) +
  tm_fill(col= "clust", palette="-YlOrRd", title="Conflict-Climate Intersection",
          legend.show = T, labels= label)+
  tm_shape(within_migrants) +
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
tmap_save(map,  dpi= 600,  height=8.3, width=11.7, units="in",scale = 1.6,
          filename=paste0(wd,"/Final results/within_Overlay.png"))
#socioeconomic profiles
#define professions
agriculture <- c("Agriculture, Fishery, and/or Forestry workers",
                 "Skilled agricultural, forestry and fishery worker (e.g. gardeners, farmers, fishers, gatherers)",
                 "Unskilled farmer", "Traditional farmer")
pastoralist <- c("Pastoralist","Pastoralist or livestock herder","Livestock herder",
                 "Livestock rearing","Livestock keeper","Livestock","Rearing reduced","Rearing induce")
elementary <- c("Elementary occupation (e.g. cleaners, mining/ construction labourers, street vendors, refuse workers)",
                "Aide conducteur","Motorbike operator","Construction worker's")
profession <- c("Professional (e.g. doctors, nurses, teachers, accountants)","Managers, Professionals , Office work (ex: public servant, NGO / UN worker)",
                "Mining","Health Professional","Education Professional")
services <- c("Waitress","Waiter at a restaurant","Waiter","vente","Vendeur de son boutique",
              "Vegetable seller","Dishwasher","Selling of milk","Casual labour","Beauty salon",
              "Shop keeper","Car wash","Domestic worker","Garage", 
              "Services and sales worker (e.g. cooks, hairdressers, protective services",
              "Service and sales workers (ex: make tea, serve food, sell at market)")
skilled_manual <- c("Carpenter","i. Skilled manual (craft, transport)")
unskilled_manual <- c("Unskilled Manual")
#within
within_profession <- fms_within[,c("survey_date","Region","main_profession_of_recent_or_current_job","other_profession_of_recent_or_current_job")]
within_profession <- within_profession[!is.na(within_profession$main_profession_of_recent_or_current_job), ]
within_profession <- subset(within_profession, !main_profession_of_recent_or_current_job == "Don’t know/ No answer")
within_profession$profession <- ""
within_profession$main_profession_of_recent_or_current_job <- ifelse(
  within_profession$main_profession_of_recent_or_current_job == "Other" & !is.na(within_profession$other_profession_of_recent_or_current_job),
  within_profession$other_profession_of_recent_or_current_job,
  within_profession$main_profession_of_recent_or_current_job
)

within_profession$main_profession_of_recent_or_current_job <- ifelse(
  within_profession$main_profession_of_recent_or_current_job == "Other, specify" & !is.na(within_profession$other_profession_of_recent_or_current_job),
  within_profession$other_profession_of_recent_or_current_job,
  within_profession$main_profession_of_recent_or_current_job
)
within_profession$profession <- ifelse(within_profession$main_profession_of_recent_or_current_job %in% pastoralist, "Pastoralist", 
                                        ifelse(within_profession$main_profession_of_recent_or_current_job %in% agriculture, "Agriculture, Fishery&Forestry",
                                               ifelse(within_profession$main_profession_of_recent_or_current_job %in% elementary, "Elementary Occupation",
                                                      ifelse(within_profession$main_profession_of_recent_or_current_job %in% profession, "Professional",
                                                             ifelse(within_profession$main_profession_of_recent_or_current_job %in% services, "Services&Sales",
                                                                    ifelse(within_profession$main_profession_of_recent_or_current_job %in% skilled_manual, "Skilled Manual",
                                                                           ifelse(within_profession$main_profession_of_recent_or_current_job %in% unskilled_manual, "Unskilled Manual",
                                                                                  "Others")))))))
within_profession_ <- within_profession[,c("Region","profession")]
within_profession_$number <- 1
within_profession_ <- aggregate(number ~ Region+profession, data=within_profession_, FUN=sum)
#convert to long format
within_prof_wide <- pivot_wider(
  data = within_profession_,
  names_from = profession,
  values_from = number
)
#replace NA value with 0
within_prof_wide[is.na(within_prof_wide)] <- 0
#merge with shp
within_pie <- merge(admin1,within_prof_wide, by="Region")
centroids <- st_centroid(within_pie)
c <- cbind(centroids, st_coordinates(centroids))
labels <- cbind(centroids, st_coordinates(centroids))
names(c)[names(c) == "Agriculture..Fishery.Forestry"] <- "Agriculture,Fishery&Forestry"
names(c)[names(c) == "Elementary.Occupation"] <- "Elementary Occupation"
names(c)[names(c) == "Services.Sales"] <- "Services&Sales"
names(c)[names(c) == "Skilled.Manual"] <- "Skilled Manual"
names(c)[names(c) == "Unskilled.Manual"] <- "Unskilled Manual"
c <- as.data.frame(c)

# Plot the map with pie charts
within_plot <- ggplot() +
  geom_sf(data = within_pie, fill = "gray80", color = "black")+
  geom_scatterpie(data = c, aes(x = X, y = Y, r = 0.5, group=Region),  # Add pie charts using centroid coordinates
  cols = c("Pastoralist", "Agriculture,Fishery&Forestry", "Professional","Elementary Occupation", "Skilled Manual","Unskilled Manual","Services&Sales","Others"),
  color = "black",
  size=0.3)+
  geom_text(data=labels,
    aes(x=X + 0.7, y=Y + 0.4,label=Region),
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
  scale_fill_manual(values = c("Agriculture,Fishery&Forestry" = "green", "Others" = "white", "Pastoralist" = "orange","Elementary Occupation" = "thistle",
                               "Professional"="purple","Services&Sales"="lavender","Skilled Manual"="black", "Unskilled Manual"="blue"))+  # Custom colors
  labs(title = "",fill="Profession")
within_plot
output <- paste0(wd,"/Final results/within_pie.png")
ggsave(output, plot = within_plot, dpi= 600,  height=8.3, width=11.7, units="in")

#plot outside HOA pie charts
outside_profession <- fms_outside[,c("survey_date","Region","main_profession_of_recent_or_current_job","other_profession_of_recent_or_current_job")]
outside_profession <- outside_profession[!is.na(outside_profession$main_profession_of_recent_or_current_job), ]
outside_profession <- subset(outside_profession, !main_profession_of_recent_or_current_job == "Don’t know/ No answer")
outside_profession$profession <- ""
outside_profession$main_profession_of_recent_or_current_job <- ifelse(
  outside_profession$main_profession_of_recent_or_current_job == "Other" & !is.na(outside_profession$other_profession_of_recent_or_current_job),
  outside_profession$other_profession_of_recent_or_current_job,
  outside_profession$main_profession_of_recent_or_current_job
)

outside_profession$main_profession_of_recent_or_current_job <- ifelse(
  outside_profession$main_profession_of_recent_or_current_job == "Other, specify" & !is.na(outside_profession$other_profession_of_recent_or_current_job),
  outside_profession$other_profession_of_recent_or_current_job,
  outside_profession$main_profession_of_recent_or_current_job
)
outside_profession$profession <- ifelse(outside_profession$main_profession_of_recent_or_current_job %in% pastoralist, "Pastoralist", 
                                       ifelse(outside_profession$main_profession_of_recent_or_current_job %in% agriculture, "Agriculture, Fishery&Forestry",
                                              ifelse(outside_profession$main_profession_of_recent_or_current_job %in% elementary, "Elementary Occupation",
                                              ifelse(outside_profession$main_profession_of_recent_or_current_job %in% profession, "Professional",
                                                     ifelse(outside_profession$main_profession_of_recent_or_current_job %in% services, "Services&Sales",
                                                            ifelse(outside_profession$main_profession_of_recent_or_current_job %in% skilled_manual, "Skilled Manual",
                                                                   ifelse(outside_profession$main_profession_of_recent_or_current_job %in% unskilled_manual, "Unskilled Manual",
                                              "Others")))))))

outside_profession_ <- outside_profession[,c("Region","profession")]
outside_profession_$number <- 1
outside_profession_ <- aggregate(number ~ Region+profession, data=outside_profession_, FUN=sum)
#convert to long format
outside_prof_wide <- pivot_wider(
  data = outside_profession_,
  names_from = profession,
  values_from = number
)
#replace NA value with 0
outside_prof_wide[is.na(outside_prof_wide)] <- 0
#merge with shp
outside_pie <- merge(admin1,outside_prof_wide, by="Region")
out_centroids <- st_centroid(outside_pie)
out_c <- cbind(out_centroids, st_coordinates(out_centroids))
labels_ <- cbind(out_centroids, st_coordinates(out_centroids))
names(out_c)[names(out_c) == "Agriculture..Fishery.Forestry"] <- "Agriculture,Fishery&Forestry"
names(out_c)[names(out_c) == "Elementary.Occupation"] <- "Elementary Occupation"
names(out_c)[names(out_c) == "Services.Sales"] <- "Services&Sales"
names(out_c)[names(out_c) == "Skilled.Manual"] <- "Skilled Manual"
names(out_c)[names(out_c) == "Unskilled.Manual"] <- "Unskilled Manual"
out_c <- as.data.frame(out_c)

# Plot the map with pie charts
outside_plot <- ggplot() +
  geom_sf(data = outside_pie, fill = "gray80", color = "black")+
  geom_scatterpie(data = out_c, aes(x = X, y = Y, r = 0.5, group=Region),  # Add pie charts using centroid coordinates
                  cols = c("Pastoralist", "Agriculture,Fishery&Forestry", "Professional","Elementary Occupation", "Skilled Manual","Unskilled Manual","Services&Sales","Others"),
                  color = "black",
                  size=0.3)+
  geom_text(data=labels_,
            aes(x=X + 0.7, y=Y + 0.4,label=Region),
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
  scale_fill_manual(values = c("Agriculture,Fishery&Forestry" = "green", "Others" = "white", "Pastoralist" = "orange","Elementary Occupation" = "thistle",
                               "Professional"="purple","Services&Sales"="lavender","Skilled Manual"="black", "Unskilled Manual"="blue"))+  # Custom colors
  labs(title = "",fill="Profession")
outside_plot
output_file <- paste0(wd,"/Final results/outside_pie.png")
ggsave(output_file, plot = outside_plot, dpi= 600,  height=8.3, width=11.7, units="in")
