# ------------------------------------------ #
# FAO-Penman-Monteith Evapotranspiration method
# By: Cesar Saavedra & Harold Achicanoy
# Adopted by Benson Kenduiywo
# ABC
# Feb. 2024
# ------------------------------------------ #

# R options and packages loading
rm(list=ls(all=TRUE))
options(warn = -1, scipen = 999)

g <- gc(reset = T); 
#install.packages("remotes")
#remotes::install_github("rspatial/luna")
list.of.packages <- c("tidyverse","terra","gtools", 'lubridate', 'geodata')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)


# Root directory
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Calculate FAO-Penman-Monteith ET function
calc_ET0 <- function(yr, mn){
  # Define output file
  outfile <- paste0(out_dir,'/monthly_evapotranspiration/ET-',yr,'-',mn,'.tif')
  cat('>>> Processing year: ',yr,', month: ',mn,'\n')
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T) # Create directory
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    # Maximum temperature files
    tx_fls <- paste0(tx_pth,'/Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    # Minimum temperature files
    tm_fls <- paste0(tm_pth,'/Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    # Average temperature files
    tav_fls <- paste0(tav_pth,'/Temperature-Air-2m-Mean-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
    tav_fls <- tav_fls[file.exists(tav_fls)]
    # Dew point temperature files
    dp_fls <- paste0(dp_pth,'/Dew-Point-Temperature-2m-Mean_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
    dp_fls <- dp_fls[file.exists(dp_fls)]
    # Solar radiation files
    sr_fls <- paste0(sr_pth,'/Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1','.nc')
    sr_fls <- sr_fls[file.exists(sr_fls)]
    # Wind speed files
    ws_fls <- paste0(ws_pth,'/Wind-Speed-10m-Mean_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1','.nc')
    ws_fls <- ws_fls[file.exists(ws_fls)]
    
    # Load maximum temperature
    tmx <- terra::rast(tx_fls)
    tmx <- tmx - 273.15 # Kelvin to Celsius
    # Load minimum temperature
    tmn <- terra::rast(tm_fls)
    tmn <- tmn - 273.15 # Kelvin to Celsius
    # Load average temperature
    tmean <- terra::rast(tav_fls)
    tmean <- tmean - 273.15 # Kelvin to Celsius
    # Load dew point temperature
    dp <- terra::rast(dp_fls)
    dp <- dp - 273.15 # Kelvin to Celsius
    # Load solar radiation
    sr <- terra::rast(sr_fls)
    sr <- sr/1000000 # Solar radiation from j to kj
    sr <- 1/2.45*(sr) # Solar radiation from kj to mm day-1
    # Load wind speed
    ws <- terra::rast(ws_fls)
    ws_2m <- ws*((4.87)/(log((67.8*10)-5.42)))
    # Load elevation data
    elev <- geodata::elevation_global(res = 10,tempdir())
    elev <- terra::resample(x = elev, y = tmx[[1]], threads = T)
    lst <- list()
    for(i in 1:length(dts)){
      print(i)
      lst[[i]] <- elev
    }
    elv <- terra::rast(lst)
    # calculate ET constants 
    p <- (101.3)*((293-0.0065*elv)/(293))^5.26   # Presion
    cps <- ((1.013*10^-3)*p)/(0.622*2.45)        # Constante psicrometrica
    etmx <- 0.6108*exp((17.27*tmx)/(tmx+237.3))  # Presion media de vapor de la saturacion
    etmn <- 0.6108*exp((17.27*tmn)/(tmn+237.3))  # Presion media de vapor de la saturacion
    es <- (etmx+etmn)/2                          # Presion media de vapor de la saturacion
    Dlt <- (4098*(0.6108*(exp((17.27*tmean)/(tmean+237.3)))))/((tmean+237.3)^2) # Delta/Pendiente de la curva de presion de saturacion de vapor
    ea <- 0.6108*(exp((17.27*dp)/(dp+237.3)))    # Presion real de vapor derivada del punto de rocio
    # Evapotranspiration function
    ET_0_idx <- function(Dlt, sr, cps, tmean, ws_2m, es){
      et_0 <- ((0.408*Dlt*sr)+(cps*(900/(tmean+273))*ws_2m*es))/(Dlt+(cps*(1+(0.34*ws_2m))))
      return(et_0)
    }
    # Calculate evapotranspiration 
    ET_0 <- terra::lapp(x = terra::sds(Dlt, sr, cps, tmean, ws_2m, es), fun = ET_0_idx) %>% sum()
    terra::writeRaster(ET_0, outfile, overwrite = T)
    
  } else {
    cat('File exist. Creating next file.\n')
  }
  cat("Done", "\n", ":)", "\n")
}

# Historical setup
yrs <- 1981:2023
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
tx_pth <- paste0(root,'/1.Data/AgERA5/2m_temperature-24_hour_maximum') # Maximum temperature
tm_pth <- paste0(root,'/1.Data/AgERA5/2m_temperature-24_hour_minimum') # Minimun temperature
tav_pth <- paste0(root,'/1.Data/AgERA5/2m_temperature-24_hour_mean')   # Mean temperature
sr_pth <- paste0(root,'/1.Data/AgERA5/solar_radiation_flux')           # Solar radiation
ws_pth <- paste0(root,'/1.Data/AgERA5/10m_wind_speed')                 # Wind speed
dp_pth <- paste0(root,'/1.Data/AgERA5/2m_dewpoint_temperature')        # Dewpoint temperature

out_dir <- out_dir <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/intermediate/spei'#paste0(root,'/agroclimExtremes/agex_raw_data')
# loop for each year and month
1:nrow(stp) |>
  purrr::map(.f = function(i){calc_ET0(yr = stp$yrs[i], mn = stp$mns[i]); gc(T)})
