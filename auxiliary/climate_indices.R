
#* TODO: Script to count the number of dry days in every month of each season
#* Requires: Specify the ISO of country of interest in line 45
#* Requires: Specify the  range of years to compute the indices in line  47
#* 
#* 
#* 
#* Author: Victor Korir
################################################################################
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

calc_ndd <- function(yr, mn){
  outfile <- paste0(out_dir, '/',iso, '/', 'season',index,'/','NDD','/NDD-',yr,'-',mn,'.tif')
  cat(outfile, "\n")
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    # Files
    
    fls <- paste0(pr_pth,'/chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    print(fls)
    print('Reading files complete---')
    fls <- fls[file.exists(fls)]
    print(fls)
    # Read precipitation data
    prc <- terra::rast(fls)
    prc <- prc %>% terra::crop(terra::ext(shp_fl)) %>% terra::mask(shp_fl)
    prc[prc == -9999] <- NA
    print('Precipitation data loading complete')
    # Calculate number of dry days
    terra::app(x   = prc,
               fun = function(x){ ndd = sum(x < 1, na.rm = T); return(ndd) },
               filename = outfile)
  }
}
pr_pth <- '//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/Chirps'

#Specify country ISO here#################################################################################################################
iso <- 'NER'
#specify the range of years to compute the indeces #########################################################################################
years <- seq(2015, 2017, by = 1)

#source('calc_NDD.R')
shp_fl <-st_read(paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/_shps/", iso, ".shp"))
ref <- "//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif"
out_dir <- paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/Hazards")

seasons <- switch(iso, "KEN" = list(season_type_1 = 1:6, season_type_2 = 7:12),
                  "SEN" = list(season_type_1 = 6:12),
                  "ETH" = list(season_type_1 = 4:10),
                  "GTM" = list(season_type_1 = 5:11),
                  "MLI" = list(season_type_1 = 6:12),
                  "NGA" = list(season_type_1 = 6:12),
                  "PHL" = list(season_type_1 = 6:12),
                  "SDN" = list(season_type_1 = 6:10),
                  "UGA" = list(season_type_1 = 3:10),
                  "SOM" = list(season_type_1 = 3:10),
                  "ZMB" = list(season_type_1 = 1:6, season_type_2 = 7:12),
                  "ZWE" = list(season_type_1 = 1:6, season_type_2 = 7:12),
                  "NER" = list(season_type_1 = 5:10, season_type_2 = c(11,12,1,2,3,4)),
                  "SSD" = list(season_type_1 = 4:11, season_type_2 = c(12,1,2,3) ),
                  "BFA" = list(seanon_type_1 = 4:10)
)



for (yr in years) {
  print(yr)
  for (season in seasons){
    index <- which(sapply(seasons, function(x) identical(x, season)))
    if (season[1] > season[length(season)]){
      start <- as.Date(paste0(yr,'-', season[1], '-', 01))
      lastday <- lubridate::days_in_month(as.Date(paste0(yr+1,'-',season[length(season)],'-01')))
      end   <- as.Date(paste0(yr+1,'-', season[length(season)], '-', lastday))
      season_mnths <- seq(from = start, to = end, by = 'month')
      for (i in 1:length(season_mnths)){
        calc_ndd(lubridate::year(season_mnths[i]), lubridate::month(season_mnths[i]))
      }
      print(start)
      print(end)
    } else{
      start <- as.Date(paste0(yr,'-', season[1], '-', 01))
      lastday <- lubridate::days_in_month(as.Date(paste0(yr,'-',season[length(season)],'-01')))
      end   <- as.Date(paste0(yr,'-', season[length(season)], '-', lastday))
      season_mnths <- seq(from = start, to = end, by = 'month')
      for (i in 1:length(season_mnths)){
        calc_ndd(lubridate::year(season_mnths[i]), lubridate::month(season_mnths[i]))
      }
      print(start)
      print(end)
      
    }
    
  }
  
}




