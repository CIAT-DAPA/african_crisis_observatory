# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Create folders structure
# Author: Harold Achicanoy
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA','GTM')

create_folders <- function(iso = 'SDN'){
  
  out <- paste0(root,'/data/',iso,'/')
  vrs <- c('accessibility','agricultural_area','child_growth_failure','climate_water_deficit','climatic_indexes','conflict','deforestation','education','evapotranspiration','flooding','irrigation','land_cover','livestock','mask','migration','net_primary_product','nightlights','pasture_area','population_density','precipitation_annual','sanitation','soil','tmax_annual','wealth_index')
  
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
