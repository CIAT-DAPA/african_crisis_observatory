# Produce daily climate tables per country - master code
# African Crisis Observatory
# By: H. Achicanoy
# Alliance Bioversity-CIAT, 2021

options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(doSNOW,terra,plyr,tidyverse,parallel,countrycode,raster,RSAGA,ncdf4,sf,future,furrr,lubridate,glue,cowsay,fst))
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/get_daily_climate/01_country_5km_srtm.R') # Get 5 km elevation raster
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/get_daily_climate/02_CPC_downscaling_process.R') # Get 5 km min and max temperature rasters
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/get_daily_climate/03_extract_daily_climate_data.R') # Final table with summarized information at the lowest administrative level

# Root directory
root <- '//dapadfs.cgiarad.org/workspace_cluster_9/Sustainable_Food_System/Grazia'

# Selected countries
geodata <- data.frame(iso = c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA'),
                      country = c('Sudan','Zimbabwe','Senegal','Mali','Nigeria','Kenya','Uganda'))

1:nrow(geodata) %>%
  purrr::map(.f = function(i){
    # Step 1: produce 5 km elevation raster
    get_5km_srtm(iso = geodata$iso[i], country = geodata$country[i])
    # Step 2: CPC temperature downscaling
    cpc_downscaling(iso = geodata$iso[i], country = geodata$country[i])
    # Step 3: produce table
    get_table_wrap(iso = geodata$iso[i], country = geodata$country[i], ncores = 5)
    return('Done\n')
  })
