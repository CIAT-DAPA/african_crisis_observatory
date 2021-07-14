options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA')

create_folders <- function(iso = 'SDN'){
  
  out <- paste0(root,'/data/',iso,'/')
  vrs <- c('accessibility','climate_water_deficit','land_cover','nightlights','population_density','precipitation_annual','tmax_annual')
  
  paste0(out,vrs) %>%
    purrr::map(.f = function(pth){
      dir.create(path = pth, F, T)
    })
  
  return('Done\n')
}
isos %>%
  purrr::map(.f = function(iso){
    create_folders(iso = iso)
  })
