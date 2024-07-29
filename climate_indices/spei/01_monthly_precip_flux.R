# ------------------------------------------ #
# Monthly precipitation
# By: Harold Achicanoy
# Adopted by Benson Kenduiywo
# WUR & ABC
# Feb. 2024
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
list.of.packages <- c("dplyr","terra", "future", "furrr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)

# Root directory
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Calculate monthly precipitation function
calc_mprec <- function(yr, mn){
  
  # Define output file
  outfile <- paste0(out_dir,'/monthly_precipitation/prec-',yr,'-',mn,'.tif')
  cat('>>> Processing year: ',yr,', month: ',mn,'\n')
  if(!file.exists(outfile)){
    # Create directory
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    # Daily precipitation files
    pr_fls <- paste0(pr_pth,'/Precipitation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    # Load daily precipitation
    prc <- terra::rast(pr_fls)
    prc <- sum(prc)
    # Write monthly precipitation file
    terra::writeRaster(prc, outfile, overwrite = T)
    gc(T)
    return(cat('Done ...\n'))
    
  }
}

pr_pth <- paste0(root,'/1.Data/AgERA5/precipitation_flux')
out_dir <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/intermediate/spei' #paste0(root,'/agroclimExtremes/agex_raw_data')

# Define setup
yrs <- 1981:2023
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) |> base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp |>
  dplyr::arrange(yrs, mns) |>
  base::as.data.frame()

# loop for each year and month
# 1:nrow(stp) |>
#   purrr::map(.f = function(i){calc_mprec(yr = stp$yrs[i], mn = stp$mns[i]); gc(T)})

ncores <- 20
future::plan(cluster, workers = ncores, gc = T)
1:nrow(stp) |>
  furrr::future_map(.f = function(i){calc_mprec(yr = stp$yrs[i], mn = stp$mns[i]); gc(T)})
future:::ClusterRegistry('stop')
