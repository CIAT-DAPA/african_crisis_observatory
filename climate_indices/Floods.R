list.of.packages <- c("tidyverse","terra","gtools", "sf", 'furrr', 'future')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)

#temp and precipation
prec_baseline <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/Historical_Climate_data_1991_2020/prec_1991_2020"
prec_baseline <- list.files(prec_baseline, pattern = "\\.tif$", full.names = TRUE)
prec_baseline <- terra::rast(prec_baseline) #resolution = 5.5km
tmax_baseline <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/Historical_Climate_data_1991_2020/tmax_1991_2020"
tmax_baseline <- list.files(tmax_baseline, pattern = "\\.tif$", full.names = TRUE)
tmax_baseline <- terra::rast(tmax_baseline)
tmin_baseline <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/Historical_Climate_data_1991_2020/tmin_1991_2020"
tmin_baseline <- list.files(tmin_baseline, pattern = "\\.tif$", full.names = TRUE)
tmin_baseline <- terra::rast(tmin_baseline)

#solar radiation
era5Dir <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/AgERA5/solar_radiation_flux'
sr_baseline <- list.files(era5Dir, pattern = "\\.tif$", full.names = TRUE)
sr_baseline <- terra:rast(sr_baseline)