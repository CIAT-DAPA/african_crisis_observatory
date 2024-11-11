#' ############################################################################################
#' Compute environmental indices using worldclim rainfall and temperature data
#' Author: Brenda Chepngetich
#' ############################################################################################
library(terra)

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
tmin <- list.files(tmin_files, pattern = "*.tif$", full.names = TRUE)
tmin <- terra::rast(tmin)
avg_temp <- (tmax + tmin) / 2
prec <- list.files(prec_files, pattern = "*.tif$", full.names = TRUE)
prec <- terra::rast(prec)

#compute monthly averages for solar radiation
era5Dir <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/AgERA5'
srd_pth <- paste0(era5Dir,'/solar_radiation_flux')
srd_fls <- gtools::mixedsort(list.files(srd_pth, pattern = '*.nc$', full.names = T))
srd_dts <- strsplit(x = srd_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
srd_dts <- strsplit(x = srd_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
srd_dts <- as.Date(srd_dts, "%Y%m%d")
yrs <- lubridate::year(srd_dts)
yrs <- yrs[yrs >= "1991" & yrs <= "2020"]
srd_fls <- srd_fls[lubridate::year(srd_dts) %in% yrs]
srd_dts <- srd_dts[lubridate::year(srd_dts) %in% yrs]
monthly_averages <- list()

for (year in 1991:2020){
  for (month in 1:12){
    month_files <- srd_fls[format(srd_dts, "%Y") == as.character(year) & format(srd_dts, "%m") == sprintf("%02d", month)]
    if (length(month_files) > 0) {
      rasters <- rast(month_files)
      monthly_mean <- mean(rasters, na.rm = TRUE)
      monthly_averages[[paste(year, month, sep = "_")]] <- monthly_mean
    }
  }
  
srd_monthly <- terra::rast(monthly_averages) 
  
  
  
}





