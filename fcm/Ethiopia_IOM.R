#Author:Brenda Chepngetich
#load required packages
library(geodata)
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
#count migrants per region
fms_region$migrants <- 1
fms_region <- aggregate(migrants ~ survey_date + Region, data = fms_region, FUN = sum)
fms_ethiopia$migrants <- 1
#create a new dataframe containing departure, number of immigrants and date
fms_ <- fms_ethiopia[,c("survey_date","departure_admin1","migrants")]
View(fms_)
unique(fms_$departure_admin1)
fms_ <- aggregate(migrants~survey_date + departure_admin1, data = fms_, FUN=sum)

#rename admin 1 column to region
colnames(fms_)[2] <- "Region"

#plotting 
ggplot(fms_, aes(x = survey_date, y = migrants,color = Region, group = Region)) +
  geom_line(size = 0.8) +
  geom_smooth(se = TRUE, method = "loess", linetype = "dashed", size = 0.7) +
  geom_point(size = 0.8) +
  labs(
    title = "Migration Trends in Ethiopia Regions (2018-2024)",
    x = "Time",
    y = "Number of Migrants"
  ) +
  theme_bw()
