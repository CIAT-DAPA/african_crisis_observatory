
#* TODO: Script to Geocode and visualize village/town names
#* Requires: Google maps API
#* Reference: https://guides.library.duke.edu/r-geospatial/geocode
#* 
#* Author: KORIR, Victor K.
################################################################################
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,ggmap,tidyverse,lubridate, readxl))

# Enter the API key here
register_google(key = "API KEY", write = TRUE)

#read IOM_FMS data
FMS_data <- readxl::read_xlsx("D:/OneDrive - CGIAR/SA_Team/korir/FCM/Somalia/New IOM Data/FMS_dataset_combined_clean_dept_and_transit.xlsx")

#Remove NAs and any character than cannot be Geocoded in village of departure
FMS_data <- FMS_data[-c(which(is.na(FMS_data$`City (Cleaned)`))),]
FMS_data <- FMS_data[,c("Date of Survey","What is the main reason for your Journey?", "Have you been forcibly displaced in your own country?"  ,"District (Cleaned)" , "City (Cleaned)" )]
names <- FMS_data$`City (Cleaned)`
names <- names[-c(which(names=='Not Specified'))]
names <- as.data.frame(names)

######### GEOCODING METHOD 1 #################################################
geo <- mutate_geocode(names, location = names, output = 'latlona')
write.csv(geo, file = 'Z:/1.Data/Palmira/IOM/Migration/geo_villages.csv')


######### GEOCODING METHOD 2 #################################################
library(tidygeocoder)
df <- tidygeocoder::geocode(names, address = names, method = "osm")
write.csv(df, file = 'Z:/1.Data/Palmira/IOM/Migration/tidygeo_villages.csv' )


########### Coalescing coordinates from the two approaches ################
gg_Geocoded <- read.csv('Z:/1.Data/Palmira/IOM/Migration/geo_villages.csv')
tidy_Geocoded <- read.csv('Z:/1.Data/Palmira/IOM/Migration/tidygeo_villages.csv')

Geocoded <- gg_Geocoded %>% mutate(lon = coalesce(tidy_Geocoded$long)) %>%
  mutate(lat = coalesce(tidy_Geocoded$lat))


######### Merging with the initial data #####################################
#Filtering for displaced persons only
FMS_data <- FMS_data[which(FMS_data$`Have you been forcibly displaced in your own country?` == 'Yes'),]


colnames(Geocoded) <- c("X","City (Cleaned)","lon","lat","address")
Geocoded <- subset(Geocoded,select = -X)
Geocoded <- distinct(Geocoded)

FMS_data <- FMS_data[-c(which(FMS_data$`City (Cleaned)` == 'Not Specified')),]
FMS_data <- FMS_data %>% mutate(ID= row_number())
FMS_geo <- FMS_data %>% 
  left_join(Geocoded, by = "City (Cleaned)") %>% 
  distinct(ID, .keep_all = TRUE)
#Removing missing coordinates
FMS_geo <- FMS_geo[-c(which(is.na(FMS_geo$lat))),]

FMS_geo$`Date of Survey` <- lubridate::year(as.Date(FMS_geo$`Date of Survey`))

Fgeo_2018 <- FMS_geo[which(FMS_geo$`Date of Survey`==2018),]
Fgeo_2019 <- FMS_geo[which(FMS_geo$`Date of Survey`==2019),]
Fgeo_2020 <- FMS_geo[which(FMS_geo$`Date of Survey`==2020),]
Fgeo_2021 <- FMS_geo[which(FMS_geo$`Date of Survey`==2021),]
Fgeo_2022 <- FMS_geo[which(FMS_geo$`Date of Survey`==2022),]
Fgeo_2023 <- FMS_geo[which(FMS_geo$`Date of Survey`==2023),]

write.csv(Fgeo_2018, file = 'Z:/1.Data/Palmira/IOM/Migration/geo_2018.csv')
write.csv(Fgeo_2019, file = 'Z:/1.Data/Palmira/IOM/Migration/geo_2019.csv')
write.csv(Fgeo_2020, file = 'Z:/1.Data/Palmira/IOM/Migration/geo_2020.csv')
write.csv(Fgeo_2021, file = 'Z:/1.Data/Palmira/IOM/Migration/geo_2021.csv')
write.csv(Fgeo_2022, file = 'Z:/1.Data/Palmira/IOM/Migration/geo_2022.csv')
write.csv(Fgeo_2023, file = 'Z:/1.Data/Palmira/IOM/Migration/geo_2023.csv')


