#' ############################################################################################
#' Compute environmental indices using worldclim rainfall and temperature data
#' Author: Brenda Chepngetich
#' ############################################################################################
rm(list=ls(all=TRUE))
list.of.packages <- c("tidyverse", "terra", "gtools", "sf", "furrr", "future", "caTools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)
tmax_files <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/Historical_Climate_data_1991_2020/tmax_1991_2020"
tmin_files <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/Historical_Climate_data_1991_2020/tmin_1991_2020"
prec_files <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/Historical_Climate_data_1991_2020/prec_1991_2020"

tmax <- list.files(tmax_files, pattern = "*.tif$", full.names = TRUE)
tmax <- terra::rast(tmax)
tmax_avg <- app(tmax, mean)
tmin <- list.files(tmin_files, pattern = "*.tif$", full.names = TRUE)
tmin <- terra::rast(tmin)
prec <- list.files(prec_files, pattern = "*.tif$", full.names = TRUE)
prec <- terra::rast(prec)
summary(tmax)


era5Dir <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/AgERA5'
srd_pth <- paste0(era5Dir,'/solar_radiation_flux')
