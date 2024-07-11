# ------------------------------------------ #
# Monthly water balance: precipitation - evapotranspiration
# By:  Harold Achicanoy
# Adopted by Benson Kenduiywo
# ABC
# Feb. 2024
# ------------------------------------------ #

# R options and packages loading
rm(list=ls(all=TRUE))


g <- gc(reset = T); 
options(warn = -1, scipen = 999)
#install.packages("remotes")
#remotes::install_github("rspatial/luna")
list.of.packages <- c("tidyverse","terra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)


# Root directory
root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/intermediate/spei'

# Historical setup
yrs <- 1981:2023
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
# Input directories
evp_pth <- paste0(root,'/monthly_evapotranspiration')
prc_pth <- paste0(root,'/monthly_precipitation')
# Output directory

out_pth <- paste0(root,'/monthly_balance')
dir.create(out_pth, recursive = T)

calc_balance <- function(yr, mn){
  evp <- terra::rast(paste0(evp_pth,'/ET-',yr,'-',mn,'.tif'))
  prc <- terra::rast(paste0(prc_pth,'/prec-',yr,'-',mn,'.tif'))
  bal <- prc - evp
  terra::writeRaster(x = bal, filename = paste0(out_pth,'/bal-',yr,'-',mn,'.tif'))
}

# loop for each year and month
1:nrow(stp) |>
  purrr::map(.f = function(i){calc_balance(yr = stp$yrs[i], mn = stp$mns[i]); gc(T)})
