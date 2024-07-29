#Author:Brenda Chepngetich
#load required packages
library(geodata)
library(tidyverse)
library(sf)
library(tmap)
library(readxl)

#read data
file <- "D:/OneDrive - CGIAR/IOM - CGIAR Climate Security Coordination/Data/OneDrive_1_7-23-2024/FMS.xlsx"
fms <- read_excel(file)
#filter to have only those originating from Ethiopia
fms <- subset(fms, departure_country == "Ethiopia")
#rename admin 1 to "region"
names(fms)[names(fms) == "departure_admin1"] <- "Region"
View(fms_region)
#analysis to visualize migration trends per region
fms_region <- fms[, c("survey_date","Region")]
fms_region$survey_date <- as.Date(fms_region$survey_date)
str(fms_region)#check data types
#remove rows whose region names have not been specified
fms_region <- subset(fms_region, Region!="Not Specified")
#correct region names
fms_region$Region[fms_region$Region == 'Hareri'] <- 'Harari'
#count migrants per region
fms_region$migrants <- 1
fms_region <- aggregate(migrants ~ survey_date + Region, data = fms_region, FUN = sum)
View(fms_region)
#plotting 
region_trends <- ggplot(fms_by_year, aes(x = year, y = migrants,color = Region, group = Region)) +
  geom_line(size = 0.8) +
  geom_smooth(se = TRUE, method = "loess", linetype = "dashed", size = 0.7) +
  geom_point(size = 0.8) +
  labs(
    title = "Migration Trends in Ethiopia (2018-2024)",
    x = "Time",
    y = "Number of Migrants"
  ) +
  theme_bw()
region_trends
#save the plot
graph <- "D:/OneDrive - CGIAR/SA_Team/Brenda/IOM/migration_trends_by_year.png"
ggsave(graph, plot = region_trends, width = 11, height = 8, dpi = 300)

#map plotting
ETH_adm1 <- st_as_sf(geodata::gadm(country = 'ETH', level = 1, version = 'latest' , path = tempdir()))
ETH_adm1 <- ETH_adm1[,"NAME_1"]
View(ETH)
unique(fms_region$Region)
#rename regions to match those in the dataframe
ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Addis Abeba'] <- 'Addis Ababa'
ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Benshangul-Gumaz'] <- 'Beneshangul Gumuz'
ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Gambela Peoples'] <- 'Gambella'
ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Southern Nations, Nationalities'] <- 'SNNP'
ETH_adm1$NAME_1[ETH_adm1$NAME_1 == 'Harari People'] <- 'Harari'
#rename
names(ETH_adm1)[names(ETH_adm1) == "NAME_1"] <- "Region"
#group migration data by year
fms_region$year <- format(fms_region$survey_date, "%Y")
fms_by_year <- aggregate(migrants ~ Region + year, data = fms_region, FUN = sum)

# fms_by_year_ <- reshape(fms_by_year,
#                      idvar = c("Region"),
#                      timevar = "year",
#                      direction = "wide")
# ETH_Migrants <- merge(ETH_adm1,fms_by_year_, by="Region")
View(ETH_Migrants)
plot(ETH_Migrants)
ETH <- merge(fms_by_year,ETH_adm1,by="Region")
unique(ETH$Region)
#map plotting


