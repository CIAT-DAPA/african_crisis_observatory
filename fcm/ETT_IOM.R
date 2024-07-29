library(geodata)
library(sf)
library(tidyverse)
library(tmap)
#Read ETT exported as CSV
ETT <- read.csv("D:/OneDrive - CGIAR/SA_Team/korir/FCM/Somalia/ETT/IOM_DTM_ETT_GedoDroughts_Tracker_R64_publish_0.csv", strip.white = TRUE)


#Cleaning
ETT <- ETT[-1,]

#Somalia admin 2 shapefile
SOM_adm2 <- st_as_sf(geodata::gadm(country = 'SOM', level = 2, version = 'latest' , path = tempdir()))
SOM_adm2_Gedo <- SOM_adm2[which(SOM_adm2$NAME_1 == 'Gedo'),]
SOM_adm2_Bay <- SOM_adm2[which(SOM_adm2$NAME_1 == 'Bay'),]

Gedo_Bay <- rbind(SOM_adm2_Bay, SOM_adm2_Gedo)

#renaming districts to match the districts and the ETT
ETT$District[ETT$District == 'Doolow'] <- 'Dolow'
ETT$District[ETT$District == 'Luuq'] <- 'Luuk'
ETT$District[ETT$District == 'Belet Xaawo'] <- 'Beled Xaawo'
ETT$District[ETT$District == 'Garbahaarey'] <- 'Garbahaaray'
ETT$District[ETT$District == 'Baardheere'] <- 'Baar-Dheere'



#further cleaning to remove unnecessary columns
ETT <- ETT %>% select(-contains('change'))
ETT <- ETT %>% select(1:10, ends_with('IDPS'))
ETT <- na.omit(ETT)
ETT <- ETT %>%
  mutate_at(vars(11:75), ~as.numeric(str_replace_all(., "[^0-9.-]", "")))



#Aggregate at district level levels
ETT_sum <- ETT %>% group_by(District) %>% summarise(across(13:74,sum))


# Specify the start date, end date, and interval
start_date <- as.Date("2022-03-09")
end_date <- as.Date("2023-05-24")
interval <- "1 week"

# Create a sequence of dates to name the columns
date_sequence <- seq(start_date, end_date, by = interval)

colnames(ETT_sum)[2:63]<- as.character(date_sequence)

baar_long <- tidyr::gather(ETT_sum, key = "Column", value = "Sum", -District)

baar_long$Column <- as.Date(baar_long$Column)

selected_data_long$Index <- seq_along(selected_data_long$Column)
selected_data_long$Date_c <-date_sequence

# Plot using ggplot2
ggplot(baar_long, aes(x = Column, y = Sum, color = District, group = District)) +
 geom_line(size = 0.8) +
  geom_smooth(se = T, method = "loess", linetype = "dashed", size = 0.7) +
  geom_point( size = 0.8,) +
  labs(title = "Emergency Trend Tracking - Gedo Region",
       x = "Time",
       y = "Number of IDPS") +
  theme_bw()

#***To pick from here tom:
#*Co-plot the trends in every disctrict
#*Map the IDPS

# Assuming your_table is your original data frame with 64 columns

row_sum_first_10 <- rowSums(ETT_sum[, 2:12], na.rm = T)
row_sum_second_10 <- rowSums(ETT_sum[, 13:23], na.rm = T)
row_sum_3_10 <- rowSums(ETT_sum[, 24:34], na.rm = T)
row_sum_4_10 <- rowSums(ETT_sum[, 35:45], na.rm = T)
row_sum_5_10 <- rowSums(ETT_sum[, 46:56], na.rm = T)
row_sum_6_10 <- rowSums(ETT_sum[, 57:63], na.rm = T)

Monthly_sums <- data.frame(ETT_sum$District, row_sum_first_10, row_sum_second_10, 
                           row_sum_3_10, row_sum_4_10, row_sum_5_10, row_sum_6_10)

Monthly_sums <- pivot_longer(Monthly_sums, !ETT_sum.District, names_to = 'week', values_to = 'IDPS' )


shp_idps <- left_join(SOM_adm2, Monthly_sums, by= c('NAME_2'= 'ETT_sum.District'))
shp_idps <- shp_idps[-c(which(is.na(shp_idps$IDPS))),]

tmap_mode("plot")
map <- tm_shape(shp_idps) +
  tm_fill(col= "IDPS", style = 'jenks', palette= "YlOrBr", title="IDPs",
          legend.show = T, popup.vars="NAME_1") + #"-YlOrRd"
  tm_text("NAME_2", col='black', size = 1.1)+
  tm_facets(by = "week", ncol = 3, free.coords = TRUE) +
  tm_shape(SOM_adm2)+
  #tm_text("NAME_2", size = 1.0, col='black', remove.overlap = TRUE)+ 
  tm_borders(col = "black")+
  tm_compass(type = "4star", position = c("left", "top")) +
  tm_scale_bar(breaks = c(0, 75, 150),lwd = 2, text.size = 2, 
               position = c("right", "bottom"))+
  tm_mouse_coordinates()+
  tm_layout(legend.outside=F, 
            legend.text.size = 0.8,
            legend.title.size= 1.2,
            legend.frame=F,
            legend.just = c("left", "top"), 
            legend.position  = c("right", "bottom"),
            legend.width = 0.65,
            inner.margins = c(0.1, 0.1, 0.06, 0.1)
  )
map
#*Write about the three datasets
#*Load and the map the different componets of the baseline data
#*Explore Upssala data


#Load climate data

#Load conflict data

#Pearson correlation
register_google