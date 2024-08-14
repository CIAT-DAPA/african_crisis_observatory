#Author:Brenda Chepngetich
"Analyze IOM FMS migration data to visualize trends, analyse co-occurrence of climate-conflict clusters with migration, analyze socioeconomic profiles of migrants.
This is done using admin level 1 (regions)"
#install and load required packages
install.packages("networkD3")
install.packages("htmlwidgets")
install.packages("webshot")
install.packages("scatterpie")
library(geodata)
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
View(fms)
#data cleaning
fms <- subset(fms, departure_country == "Ethiopia")
unique(fms$departure_admin3) #TOTALS VALUES = 110,066
adm1_count <- sum(!fms$departure_admin1 %in% exclude & !is.na(fms$departure_admin1)) #107,688
adm2_count <- sum(!fms$departure_admin2 %in% exclude & !is.na(fms$departure_admin2)) #92,682
adm3_count <- sum(!fms$departure_admin3 %in% exclude & !is.na(fms$departure_admin3)) #80,488
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
View(snnp)
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
within_by_date <- aggregate(migrants ~ survey_date + Region, data = fms_within, FUN = sum)
within_by_year <- aggregate(migrants ~ Region + year, data = fms_within, FUN = sum)
#fms outside
outside_by_date <- aggregate(migrants ~ survey_date + Region, data = fms_outside, FUN = sum)
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

admin1 <- st_read(paste0(wd, "/data/eth_admbnda_admins_csa_bofedb_2024.shp/eth_admbnda_adm1_csa_bofedb_2024.shp"))
admin1 <- admin1[ ,"admin1Name"]
names(admin1)[names(admin1) == "admin1Name"] <- "Region"
plot(admin1)
#to do: add admin 2 data

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
  tm_layout(main.title = "Ethiopia migration within HOA",
            legend.outside = T,
            asp = 1.4,
            # inner.margins = c(0,0,0,0)
            )
map
tmap_save(map,  dpi= 300,  height=8.3, width=11.7, units="in",scale = 1.6,
          filename=paste0(wd,"/Results/within/facets.png"))


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
Northern_Europe = c("Denmark", "Estonia", "Finland", "Iceland", "Ireland", "Latvia", "Lithuania", "Norway", "Sweden", "United Kingdom")
Southern_Europe = c("Albania", "Andorra", "Bosnia and Herzegovina", "Croatia", "Greece", "Italy", "Kosovo", "Malta", "Montenegro", "North Macedonia", "Portugal", "San Marino", "Serbia", "Slovenia", "Spain", "Vatican City")
Western_Europe = c("Austria", "Belgium", "France", "Germany", "Liechtenstein", "Luxembourg", "Monaco", "Netherlands", "Switzerland")
Europe <- c(Eastern_Europe,Northern_Europe,Southern_Europe,Western_Europe,"UNK","Cayman Islands","Europe","U.K. of Great Britain and Northern Ireland")
Central_Asia = c("Kazakhstan", "Kyrgyzstan", "Tajikistan", "Turkmenistan", "Uzbekistan")
East_Asia = c("China", "Japan", "Mongolia", "North Korea", "South Korea", "Taiwan")
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
nodes <- data.frame(name = unique(c(as.character(outside_sankey_$Origin), 
                                    as.character(outside_sankey_$region))))
View(links)

# Generate links dataframe using indices
links <- data.frame(
  source = match(outside_sankey_$Origin, nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  target = match(outside_sankey_$region, nodes$name) - 1,  # match returns positions, subtract 1 for zero-indexing
  value = outside_sankey_$migrants
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
file_path <- paste0(wd,"/Results/outside/sankey_plot.html")
saveWidget(sankey_, file_path)
#initialize webshot
webshot::install_phantomjs()

# Save as PNG
webshot(file_path, paste0(wd,"/Results/outside/sankey_plot.png"))

#overlay climate_conflict hotspots with migration
#get total migrants over the years
within_migrants <- merge(shp,aggregate(migrants ~ Region, data = within_by_year, FUN = sum))
outside_migrants <- merge(shp,aggregate(migrants ~ Region, data = outside_by_year, FUN = sum))
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
View(clim_conflict)
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
            main.title = "Migration flows within HOA",
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
tmap_save(map,  dpi= 300,  height=8.3, width=11.7, units="in",scale = 1.6,
          filename=paste0(wd,"/Results/within/Overlay.png"))
#socioeconomic profiles
#within
within_profession <- fms_within[,c("survey_date","Region","main_profession_of_recent_or_current_job","other_profession_of_recent_or_current_job")]
View(within_profession)
within_profession <- within_profession[!is.na(within_profession$main_profession_of_recent_or_current_job), ]
within_profession <- subset(within_profession, !main_profession_of_recent_or_current_job == "Don’t know/ No answer")
within_profession$profession <- ""
within_profession$main_profession_of_recent_or_current_job <- ifelse(
  within_profession$main_profession_of_recent_or_current_job == "Other" & !is.na(within_profession$other_profession_of_recent_or_current_job),
  within_profession$other_profession_of_recent_or_current_job,
  within_profession$main_profession_of_recent_or_current_job
)

within_profession$profession <- ifelse(
  within_profession$main_profession_of_recent_or_current_job == "Other, specify" & !is.na(within_profession$other_profession_of_recent_or_current_job),
  within_profession$other_profession_of_recent_or_current_job,
  within_profession$profession
)

# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Elementary occupation (e.g. cleaners, mining/ construction labourers, street vendors, refuse workers)"] <- "Elementary Occupation"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Plant and machine operator, assembler (e.g. truck/ bus drivers, mining/ rubber machine operators)"] <- "Assembling,Plant&Machine Operation"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Technician and associate professional (e.g. sales and purchasing agents, religious associate professionals)"] <- "Technician&Associates"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Services and sales worker (e.g. cooks, hairdressers, protective services)"] <- "Services&Sales"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Service and sales workers (ex: make tea, serve food, sell at market)"] <- "Services&Sales"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Craft and related trades worker (e.g. metal workers, repairers, woodworkers, electronic installers)"] <- "Trade&Craft"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Education Professional"] <- "Professional"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Manager (e.g. directors, senior officials)"] <- "Management"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Clerical support worker (e.g. general secretaries, customer service clerks)"] <- "Clerical Work"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Professional (e.g. doctors, nurses, teachers, accountants)"] <- "Professional"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Health Professional"] <- "Professional"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Armed forces occupation"] <- "Armed Forces"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Government Civil Servants and Administrators"] <- "Civil Servants"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Managers, Professionals , Office work (ex: public servant, NGO / UN worker)"] <- "Professional"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
#                                "Unskilled Manual"] <- "Unskilled Manual"
# within_profession$profession[within_profession$main_profession_of_recent_or_current_job == 
                               # "i. Skilled manual (craft, transport)"] <- "Skilled Manual"
# elementary <- c("Aide conducteur","Motorbike operator","Construction worker's")
# services <- c("Waitress","Waiter at a restaurant","Waiter","vente","Vendeur de son boutique",
#               "Vegetable seller","Dishwasher","Selling of milk","Casual labour","Beauty salon",
#               "Shop keeper","Car wash","Domestic worker","Garage")
# profession <- c("Mining")
# business <- c("Small Business woman","Business Woman","Small Business Owner","Business woman")
agriculture <- c("Agriculture, Fishery, and/or Forestry workers",
                 "Skilled agricultural, forestry and fishery worker (e.g. gardeners, farmers, fishers, gatherers)",
                 "Unskilled farmer", "Traditional farmer")
# skilled_manual <- c("Carpenter")
pastoralist <- c("Pastoralist","Pastoralist or livestock herder","Livestock herder",
                 "Livestock rearing","Livestock keeper","Livestock","Rearing reduced","Rearing induce")
# unwanted <- c("Employed","Nongovernmental employed","UN Worker","Doker","Housewife")


within_profession$profession <- ifelse(within_profession$main_profession_of_recent_or_current_job %in% pastoralist, "Pastoralist", 
                                ifelse(within_profession$main_profession_of_recent_or_current_job %in% agriculture, "Agriculture, Fishery&Forestry",
                                       "Others"))
within_profession_ <- within_profession[,c("Region","profession")]
unique(within_profession$profession)
within_profession_$number <- 1
within_profession_ <- aggregate(number ~ Region+profession, data=within_profession_, FUN=sum)
percentages <- within_total$number / sum(within_total$number) * 100
within_prof_wide <- pivot_wider(
  data = within_profession_,
  names_from = profession,
  values_from = number
)
names(within_pie)
View(within_pie)
within_pie <- merge(admin1,within_prof_wide, by="Region")
centroids <- st_centroid(within_pie)
plot(centroids)
c <- cbind(centroids, st_coordinates(centroids))
c$Pastoralist <- as.numeric(c$Pastoralist)
c$X <- round(c$X)
View(c)
c$Pastoralist[is.na(c$Pastoralist)] <- 0
c$Others[is.na(c$Others)] <- 0
# Plot the map with pie charts

ggplot() +
  geom_sf(data = within_pie, fill = "gray80", color = "white") +  # Plot the world map
  geom_scatterpie(data = c, aes(x = X, y = Y, r = 2),  # Add pie charts using centroid coordinates
                  cols = c("Pastoralist", "Agriculture..Fishery.Forestry", "Others"),
                  color = "black") +  # Border color for the pie charts
  # theme_minimal() +  # Apply a minimal theme
  labs(title = "Distribution of Migrants occupation in Selected African Countries") +  # Add title
  scale_fill_manual(values = c("Agriculture..Fishery.Forestry" = "green", "Others" = "blue", "Pastoralist" = "orange"))  # Custom colors

# Create a tmap plot with pie charts
tmap_mode("view")  # Set to interactive view mode (use "plot" for static mode)

map <- tm_shape(within_pie) +
       tm_borders() +  
       tm_symbols(size = 0.5,  # Size of the symbols (pie charts)
             shape = 21,  # Circle shape
             pies = c("Pastoralist", "Agriculture, Fishery&Forestry", "Others"),  # Data for pie charts
             col = c("red", "green", "blue"),  # Colors for the categories
             border.col = "black",  # Border color of pie charts
             border.lwd = 0.5) +  # Border line width
      tm_layout(title = "Migrants Occupation")
#total profession pie chart
png(paste0(wd,"/Results/within/total_pie.png"), width = 2500, height = 2000)
#plotting area for i row and 1 column
par(mfrow = c(1, 1))
pie(
  within_total$number, 
  labels = percentages, 
  main = "Migrants within HOA Occupation",
  col = rainbow(length(within_total$profession)),
  radius = 1,
  font=2
)
legend("topright", 
       legend = within_total$profession, 
       fill = rainbow(length(within_total$profession)), 
       cex = 1,  # Size of the text in the legend
       bty = "n")
dev.off()
# Loop through each region and create a pie chart
#plot pie charts
png(paste0(wd,"/Results/within/region_pie.png"), width = 2500, height = 2000)
# Set up the plotting area for 5 rows and 3 columns for regions
par(mfrow = c(5, 3), mar = c(2, 2, 2, 2) + 0.1)
regions <- unique(within_profession_region$Region)

for (region in regions) {
  region_data <- subset(within_profession_region, Region == region)
  
  # Create a pie chart
  pie(
    region_data$number, 
    labels = NA, 
    main = region,
    col = rainbow(length(region_data$profession)),
    radius=1,
    cex = 2
  )
  legend("topright", 
         legend = within_total$profession, 
         fill = rainbow(length(within_total$profession)), 
         cex = 1,  # Size of the text in the legend
         bty = "n")
}
dev.off()
