options(warn = -1, scipen = 999)
suppressMessages(if(!require("pacman")) install.packages("pacman"))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend, vroom))

source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Main functions

# Example: cwd_path <- 'D:/OneDrive - CGIAR/climate_water_deficit'
cwd_path <- 'Your_own_path/climate_water_deficit'

# Summarize Climate Water Deficit time series
iso <- 'ZWE' # ISO 3 country code
shp <- lowest_gadm(iso = iso, out = NULL) # Download from GADM the lowest administrative level for the country
mxL <- grep(pattern = '^NAME_', x = names(shp@data), value = T) %>% .[length(.)] %>% readr::parse_number(.)
r <- raster::stack(list.files(path = cwd_path, pattern = '.nc$', full.names = T))
names(r) <- seq(from = as.Date('1990-01-01'), to = as.Date('2020-12-01'), by = 'month')
r <- r * 0.1
df <- cbind(shp@data[,grep(pattern = '^NAME_', x = names(shp@data), value = T)],
            exactextractr::exact_extract(x = r, y = shp, c('mean','median','min','max','stdev')))

vroom::vroom_write(x = df, file = paste0(cwd_path,'/',iso,'_cwd.csv'), delim = ',')
