options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

# Summarize climate water deficit and SPI time series

iso <- 'ZWE'
shp <- raster::shapefile(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
mxL <- grep(pattern = '^NAME_', x = names(shp@data), value = T) %>% .[length(.)] %>% readr::parse_number(.)
r <- raster::stack(list.files(path = "C:/Users/haachicanoy/Downloads/SPI_ImageCollection", pattern = '.tif$', full.names = T))
# r <- raster::stack(paste0("C:/Users/haachicanoy/Downloads/def_zwe",1:3,".nc"))
names(r) <- seq(from = as.Date('1990-01-01'), to = as.Date('2020-12-01'), by = 'month')
# r[] <- r[] * 0.1 # Just for climate water deficit
df <- cbind(shp@data[,grep(pattern = '^NAME_', x = names(shp@data), value = T)],
            exactextractr::exact_extract(x = r, y = shp, c('mean','median','min','max','stdev')))

names(df)
tst <- df %>%
  dplyr::select(NAME_0,NAME_1,NAME_2,
                mean.X1990.01.01:mean.X2020.12.01)
plot(x = 1:(ncol(tst)-3), y = as.numeric(tst[1, grep(pattern = 'mean', x = names(tst))]), ty = 'l', ylim = range(as.matrix(tst[,grep(pattern = 'mean', x = names(tst))])))
for(i in 2:nrow(tst)){
  lines(x = 1:(ncol(tst)-3), y = as.numeric(tst[i, grep(pattern = 'mean', x = names(tst))]))
}

vroom::vroom_write(x = df, file = paste0('D:/',iso,'_spi.csv'), delim = ',')
# vroom::vroom_write(x = df, file = paste0('D:/',iso,'_climate_water_deficit.csv'), delim = ',')
