# Produce summarized metrics per admin level - master code
# African Crisis Observatory
# By: H. Achicanoy
# Alliance Bioversity-CIAT, 2021

options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, terra, exactextractr, vroom))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

iso <- 'KEN'
shp <- raster::shapefile(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
mxL <- grep(pattern = '^NAME_', x = names(shp@data), value = T) %>% .[length(.)] %>% readr::parse_number(.)
vrs <- list.files(path = paste0(root,'/data/',iso), pattern = '*.tif$', recursive = T, full.names = T)
vrs <- vrs[-grep(pattern = '_results', x = vrs)]

znl <- vrs %>%
  purrr::map(.f = function(vr){
    r <- raster::raster(x = vr)
    df <- cbind(shp@data[,grep(pattern = '^NAME_', x = names(shp@data), value = T)],
                exactextractr::exact_extract(x = r, y = shp, c('mean','median','min','max','stdev')))
    df$variable <- names(r)
    return(df)
  })
znl <- znl %>%
  dplyr::bind_rows()
znl <- znl %>%
  tidyr::pivot_longer(cols = mean:stdev, names_to = "metric", values_to = "value") %>%
  tidyr::pivot_wider(names_from = 'variable', values_from = 'value')

vroom::vroom_write(x = znl, file = paste0('D:/',iso,'_stats.csv'), delim = ',')
