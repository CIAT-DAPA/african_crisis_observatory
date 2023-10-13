suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend, vegan, VGAM))

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/AWCPTF.R')


iso <- 'MOZ'


#' Compute a PDF from  Pareto distribution.
#' @param x (data.frame): dataframe for country relative wealth index
#' @param theta 
#' @param gini (numeric): GINI index for selected country
#' @param gdppc
ICDF <- function(x, theta, gini, gdppc){
  
  alpha = (1+gini)/(2*gini)
  thr <- (1-(1/alpha))*gdppc
  
  dens_pareto <-  VGAM::qpareto(p = x , 
                                scale = thr, 
                                shape = alpha)  #(x/ (theta+y))^alpha inverse cdf  #(alpha*theta*x^(alpha - 1))/((theta+x)^(alpha+1)) - density
  
  sd_log <- sqrt(2)*VGAM::probitlink(theta = (gini+1)/2)
  mean_log <-  log(gdppc)-((sd_log)^2/2)
  
  dens_logN <- qlnorm(p = x, 
                      meanlog = mean_log, 
                      sdlog = sd_log)
  
  return(dens_pareto^(0.32)*dens_logN^(1-0.32))
  
}


#' Function to download or just load country shapefile
#' @param iso (character): COuntry ISO-3 code
#' @return Country shapefile
#' .../_shps/{iso}.shp
get_iso_shapefile <- function(iso){
  if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
    dir.create(path = dirname(paste0(root,'/data/',iso,'/_shps/',iso,'.shp')), recursive = TRUE)
    shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  } else {
    shp <- raster::shapefile(x = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  }
  return(shp)
}
shp <- get_iso_shapefile(iso = iso)


#' Load global variables
fld <- terra::rast(paste0(root,'/data/_global/flooding/fl_frequency.tif'))
mgr <- terra::rast(paste0(root,'/data/_global/migration/1990-2000/30arcsec-net-migration-1990-2000.tif'))
irr <- terra::rast(paste0(root,'/data/_global/irrigation/irrigation.tif'))
ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
crp <- terra::rast(paste0(root,'/data/_global/cropland/Cropland2000_5m.tif')) # Read global cropland raster
pas <- terra::rast(paste0(root,'/data/_global/cropland/Pasture2000_5m.tif'))  # Read global pasture area raster
smm <- crp + pas # Compute cropland and pasture area raster



#' Function for get Child growth failure variables from global raster file
#' @param iso (character): country ISO-3 code
#' @returns Create 90th percentile, Median, Stdv and Trend raster files for stunting, Underweight, Wasting
#' Output files are:
#' .../child_growth_failure/medn_stunting.tif, 
#' .../child_growth_failure/p90_stunting.tif, 
#' .../child_growth_failure/trnd_stunting.tif
#' .../child_growth_failure/medn_underweight.tif, 
#' .../child_growth_failure/p90_underweight.tif, 
#' .../child_growth_failure/trnd_underweight.tif
#' .../child_growth_failure/medn_wasting.tif, 
#' .../child_growth_failure/medn_wasting.tif, 
#' .../child_growth_failure/medn_wasting.tif
get_chld_grwht_fail <- function(iso = 'SDN'){
  

  
  # Reference raster at 1 km
  ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
  shpr <- terra::rasterize(x = shp, y = ref, field = 'key')
  
  # ----------------------------------------------------------------------- #
  # Child growth failure variables: mean percentage (time series 2000-2019)
  # ----------------------------------------------------------------------- #
  
  # Stunting
  stntg <- list.files(path = paste0(root,'/data/_global/prevalence_stunting'), pattern = '*.TIF$', full.names = T) %>% terra::rast(.)
  names(stntg) <- paste0('yr',2000:2019)
  
  # Underweight
  uwght <- list.files(path = paste0(root,'/data/_global/prevalence_underweight'), pattern = '*.TIF$', full.names = T) %>% terra::rast(.)
  names(uwght) <- paste0('yr',2000:2019)
  
  # Wasting
  wstng <- list.files(path = paste0(root,'/data/_global/prevalence_wasting'), pattern = '*.TIF$', full.names = T) %>% terra::rast(.)
  names(wstng) <- paste0('yr',2000:2019)
  
  # ------------------------------------ #
  # 90th percentile
  # ------------------------------------ #
  p90stntg <- quantile(stntg, probs = 0.9)
  p90Avguwght <- quantile(uwght, probs = 0.9)
  p90wstng <- quantile(wstng, probs = 0.9)
  
  
  # ------------------------------------ #
  # Median
  # ------------------------------------ #
  
  mstntg <- median(stntg)
  muwght <- median(uwght)
  mwstng <- median(wstng)
  
  # ------------------------------------ #
  # Coefficient of variation
  # ------------------------------------ #
  
  vstntg <- terra::stdev(stntg)/mean(stntg)
  vuwght <- terra::stdev(uwght)/mean(uwght)
  vwstng <- terra::stdev(wstng)/mean(wstng)
  
  # Crop Median
  mstntg_crp <- terra::crop(x = mstntg, terra::ext(shp))
  mstntg_crp <- terra::resample(x = mstntg_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/child_growth_failure/medn_stunting.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = mstntg_crp, out, overwrite=TRUE)
  
  muwght_crp <- terra::crop(x = muwght, terra::ext(shp))
  muwght_crp <- terra::resample(x = muwght_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/child_growth_failure/medn_underweight.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = muwght_crp, out, overwrite=TRUE)
  
  mwstng_crp <- terra::crop(x = mwstng, terra::ext(shp))
  mwstng_crp <- terra::resample(x = mwstng_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/child_growth_failure/medn_wasting.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = mwstng_crp, out, overwrite=TRUE)
  
  # 90th percentile
  p90stntg_crp <- terra::crop(x = p90stntg, terra::ext(shp))
  p90stntg_crp <- terra::resample(x = p90stntg_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/child_growth_failure/p90_stunting.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = p90stntg_crp, out, overwrite=TRUE)
  
  p90Avguwght_crp <- terra::crop(x = p90Avguwght, terra::ext(shp))
  p90Avguwght_crp <- terra::resample(x = p90Avguwght_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/child_growth_failure/p90_underweight.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = p90Avguwght_crp, out, overwrite=TRUE)
  
  p90wstng_crp <- terra::crop(x = p90wstng, terra::ext(shp))
  p90wstng_crp <- terra::resample(x = p90wstng_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/child_growth_failure/p90_wasting.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = p90wstng_crp, out, overwrite=TRUE) 
  
  
  # Crop Coefficient of variation
  vstntg_crp <- terra::crop(x = vstntg, terra::ext(shp))
  vstntg_crp <- terra::resample(x = vstntg_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/child_growth_failure/cvar_stunting.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = vstntg_crp, out, overwrite=TRUE) 
  
  vuwght_crp <- terra::crop(x = vuwght, terra::ext(shp))
  vuwght_crp <- terra::resample(x = vuwght_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/child_growth_failure/cvar_underweight.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = vuwght_crp, out, overwrite=TRUE) 
  
  vwstng_crp <- terra::crop(x = vwstng, terra::ext(shp))
  vwstng_crp <- terra::resample(x = vwstng_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/child_growth_failure/cvar_wasting.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = vwstng_crp, out, overwrite=TRUE) 
  
  # Crop Trend
  out <- paste0(root,'/data/',iso,'/child_growth_failure/trnd_stunting.tif')
  if(!file.exists(out)){ 
    stntg_crp <- terra::crop(x = stntg, terra::ext(shp))
    tstntg_crp <- terra::app(x = stntg_crp, fun = function(x){
      x <- as.numeric(na.omit(x))
      if(length(x) > 1){
        y <- trend::sens.slope(x)$estimates
      } else {
        y <- NA
      }
      return(y)
    })
    tstntg_crp <- terra::resample(x = tstntg_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
    
    dir.create(dirname(out), showWarnings = F, recursive = T)
    terra::writeRaster(x = tstntg_crp, out, overwrite=TRUE) }
  
  out <- paste0(root,'/data/',iso,'/child_growth_failure/trnd_underweight.tif')
  if(!file.exists(out)){
    uwght_crp <- terra::crop(x = uwght, terra::ext(shp))
    tuwght_crp <- terra::app(x = uwght_crp, fun = function(x){
      x <- as.numeric(na.omit(x))
      if(length(x) > 1){
        y <- trend::sens.slope(x)$estimates
      } else {
        y <- NA
      }
      return(y)
    })
    tuwght_crp <- terra::resample(x = tuwght_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
    
    dir.create(dirname(out), showWarnings = F, recursive = T)
    terra::writeRaster(x = tuwght_crp, out) }
  
  out <- paste0(root,'/data/',iso,'/child_growth_failure/trnd_wasting.tif')
  #if(!file.exists(out)){ 
  wstng_crp <- terra::crop(x = wstng, terra::ext(shp))
  twstng_crp <- terra::app(x = wstng_crp, fun = function(x){
    x <- as.numeric(na.omit(x))
    if(length(x) > 1){
      y <- trend::sens.slope(x)$estimates
    } else {
      y <- NA
    }
    return(y)
  })
  twstng_crp <- terra::resample(x = twstng_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = twstng_crp, out, overwrite=TRUE) #}
  
}
#'run
get_chld_grwht_fail(iso = iso)

#'Function for get flooding variable  
#'@param iso (character): COuntry ISO-3 code
#'@return Raster file of flooding stored in path
#' .../flooding/flood.tif
get_flooding_var <- function(iso){
  shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
  ref <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
  shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
  
  crd <- shpr %>% terra::as.data.frame(xy = T, cells = T, na.rm = T)
  crd <- cbind(crd, terra::extract(x = fld, y = crd[,c('x','y')]))
  names(crd)[ncol(crd)] <- 'Flood'
  crd <- crd[,c('x','y','Flood')]
  
  tryCatch(expr = {rst <- terra::rast(x = crd, type = 'xyz', crs = terra::crs(shp))},
           error = function(e){cat('Terra format failed\n')},
           finally = {
             rst <- raster::rasterFromXYZ(xyz = crd, crs = raster::crs(shp))
             rst <- terra::rast(rst)
           })
  
  out <- paste0(root,'/data/',iso,'/flooding/flood.tif')
  
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = rst, out, overwrite=T) }
}
#'run
get_flooding_var(iso = iso)
#'Function for get recent migration variable  
#'@param iso (character): COuntry ISO-3 code
#'@return Raster file of recent migration stored in path
#' .../migration/rcnt_migration.tif        
get_recent_migration <- function(iso){
  shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
  ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
  shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
  
  smm <- mgr %>%
    terra::crop(shpr) %>%
    terra::resample(shpr) %>%
    terra::mask(shpr)
  
  out <- paste0(root,'/data/',iso,'/migration/rcnt_migration.tif')
  
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = smm, out) }
}
#'run
get_recent_migration(iso = iso)
#'Function to download gender population data 
#'@param iso (character): COuntry ISO-3 code
#'@param root (character): path to CSO root folder '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'
#'@param method (character): method for downloading data one of 'wget', 'curl', 'libcurl' or 'auto'
#'@return Rasters file of recent female and Male population
#' .../gender_population/{iso}_female_population.tif
#'.../gender_population/{iso}_male_population.tif
get_gender_pop_var <- function(iso, root, method = 'auto'){
  
  root_dir <- paste0(root, "/data")
  ####################################################
  ######### female population #######################
  ##################################################
  
  file_names <- paste0(tolower(iso),"_f_", c(0,1, seq(5,80, by = 5)), "_2020_constrained.tif" )
  cat(" >>> Downloading data for Female population \n")
  for(i in file_names){
    
    download.file(paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_Constrained/2020/",toupper(iso),"//",i), 
                  paste0(root_dir, "/", iso, "/gender_population/", i),
                  method = method
    )
    Sys.sleep(30)
  }
  
  stk <- list.files(paste0(root_dir, "/", iso, "/gender_population"), pattern = "_f_", full.names = T) %>%
    purrr::map(., raster)
  
  sum_r <- sum(raster::stack(stk), na.rm = T)
  sum_r[sum_r[] <= 1] <- NA
  writeRaster(sum_r,paste0(root_dir, "/", iso, "/gender_population/", iso, "_female_population.tif"), overwirte = T)
  
  
  #####################################
  ##### Male population ##############
  ###################################
  
  file_names_male <- paste0(tolower(iso),"_m_", c(0,1, seq(5,80, by = 5)), "_2020_constrained.tif" )
  cat(" >>> Downloading data for Male population \n")
  
  for(i in file_names_male){
    
    download.file(paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_Constrained/2020/",toupper(iso),"//",i), 
                  
                  paste0(root_dir, "/", iso, "/gender_population/", i),
                  method = method
    )
    Sys.sleep(30)
  }
  
  
  
  stk <- list.files(paste0(root_dir, "/", iso, "/gender_population"), pattern = "_m_", full.names = T) %>%
    purrr::map(., raster)
  
  sum_r <- sum(raster::stack(stk), na.rm = T)
  sum_r[sum_r[] <= 1] <- NA
  writeRaster(sum_r,paste0(root_dir, "/", iso, "/gender_population/",iso , "_male_population.tif"), overwirte = T)
  
  
}
#'run
get_gender_pop_var(root = root, iso = iso, method = 'wget')

#' Function to get Sanitation variables shuch as
#' Piped water and sanitation facilities (median, Average and trend)
#' @param iso (character): Country ISO-3 code
#' @returns Raster files for sanitation facilities
#' .../sanitation/medn_piped_water.tif
#' .../sanitation/avg_piped_water.tif
#' .../sanitation/tren_piped_water.tif
#' .../sanitation/medn_sanitation_facilities.tif
#' .../sanitation/avg_sanitation_facilities.tif
#' .../sanitation/tren_sanitation_facilities.tif
get_sanit_vars <- function(iso = 'PHL'){
  
  # Load the country lowest administrative level shapefile
  if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
    dir.create(path = dirname(paste0(root,'/data/',iso,'/_shps/',iso,'.shp')), recursive = TRUE)
    shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  } else {
    shp <- raster::shapefile(x = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  }
  
  # Reference raster at 1 km
  ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
  shpr <- terra::rasterize(x = shp, y = ref, field = 'key')
  
  # ----------------------------------------------------------------- #
  # Sanitation facilities: mean of percentage (time series 2000-2017)
  # ----------------------------------------------------------------- #
  # Piped water
  pp_wtr <- list.files(path = paste0(root,'/data/_global/piped_water'), pattern = '*.TIF$', full.names = T) %>% terra::rast(.)
  names(pp_wtr) <- paste0('yr',2000:2017)
  # Sanitation facilities
  st_flt <- list.files(path = paste0(root,'/data/_global/sanitation_facilities'), pattern = '*.TIF$', full.names = T) %>% terra::rast(.)
  names(st_flt) <- paste0('yr',2000:2017)
  
  # ------------------------------------ #
  # average
  # ------------------------------------ #
  avgpp_wtr <- mean(pp_wtr)
  avgst_flt <- mean(st_flt)
  
  
  
  # ------------------------------------ #
  # Median
  # ------------------------------------ #
  
  mpp_wtr <- median(pp_wtr)
  mst_flt <- median(st_flt)
  
  # ------------------------------------ #
  # Coefficient of variation
  # ------------------------------------ #
  
  vpp_wtr <- terra::stdev(pp_wtr)/mean(pp_wtr)
  vst_flt <- terra::stdev(st_flt)/mean(st_flt)
  
  # Crop Median
  mpp_wtr_crp <- terra::crop(x = mpp_wtr, terra::ext(shp))
  mpp_wtr_crp <- terra::resample(x = mpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/medn_piped_water.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mpp_wtr_crp, out) }
  
  mst_flt_crp <- terra::crop(x = mst_flt, terra::ext(shp))
  mst_flt_crp <- terra::resample(x = mst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/medn_sanitation_facilities.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mst_flt_crp, out) }
  
  # Crop average
  avgpp_wtr_crp <- terra::crop(x = avgpp_wtr, terra::ext(shp))
  avgpp_wtr_crp <- terra::resample(x = avgpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/avg_piped_water.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mpp_wtr_crp, out) }
  
  avgst_flt_crp <- terra::crop(x = avgst_flt, terra::ext(shp))
  avgst_flt_crp <- terra::resample(x = avgst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/avg_sanitation_facilities.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mst_flt_crp, out) }
  
  
  # Crop Coefficient of variation
  vpp_wtr_crp <- terra::crop(x = vpp_wtr, terra::ext(shp))
  vpp_wtr_crp <- terra::resample(x = vpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/cvar_piped_water.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vpp_wtr_crp, out) }
  
  vst_flt_crp <- terra::crop(x = vst_flt, terra::ext(shp))
  vst_flt_crp <- terra::resample(x = vst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/sanitation/cvar_sanitation_facilities.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vst_flt_crp, out) }
  
  # Crop Trend
  out <- paste0(root,'/data/',iso,'/sanitation/trnd_piped_water.tif')
  if(!file.exists(out)){ 
    pp_wtr_crp <- terra::crop(x = pp_wtr, terra::ext(shp))
    tpp_wtr_crp <- terra::app(x = pp_wtr_crp, fun = function(x){
      x <- as.numeric(na.omit(x))
      if(length(x) > 1){
        y <- trend::sens.slope(x)$estimates
      } else {
        y <- NA
      }
      return(y)
    })
    tpp_wtr_crp <- terra::resample(x = tpp_wtr_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
    
    dir.create(dirname(out), showWarnings = F, recursive = T)
    
    terra::writeRaster(x = tpp_wtr_crp, out) }
  
  
  
  
  out <- paste0(root,'/data/',iso,'/sanitation/trnd_sanitation_facilities.tif')
  if(!file.exists(out)){
    st_flt_crp <- terra::crop(x = st_flt, terra::ext(shp))
    tst_flt_crp <- terra::app(x = st_flt_crp, fun = function(x){
      x <- as.numeric(na.omit(x))
      if(length(x) > 1){
        y <- trend::sens.slope(x)$estimates
      } else {
        y <- NA
      }
      return(y)
    })
    tst_flt_crp <- terra::resample(x = tst_flt_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
    
    dir.create(dirname(out), showWarnings = F, recursive = T)
    terra::writeRaster(x = tst_flt_crp, out, overwrite = TRUE) }
  
}
#'run
get_sanit_vars(iso = iso)


#' Get irrigation variable 
#' @param iso (character): Country ISO-3 code
#' @return Raster file for irrigated areas
#' .../irrigation/irrigation.tif
get_irrigation_var <- function(iso){
  shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
  ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
  shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
  
  irr <- irr %>%
    terra::crop(shpr) %>%
    terra::resample(shpr) %>%
    terra::mask(shpr)
  
  out <- paste0(root,'/data/',iso,'/irrigation/irrigation.tif')
  
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = irr, out, overwrite = T) }
}
#'run
get_irrigation_var(iso = iso)

#'Get  diversity livestock variable
#' @param iso (Character): Country ISO-3 code
#' @return Raster file for livestock diversity variable
#' .../livestock/lvst_diver.tif
get_livestock_var <- function(iso){
  
  pth <- paste0(root,'/data/_global/livestock')
  lvst <- list.files(path = pth, pattern = '*_2010_Da.tif$', full.names = T, recursive = T) %>%
    terra::rast() %>%
    terra::crop(shpr) %>%
    terra::resample(shpr) %>%
    terra::mask(shpr)
  lvst <- round(lvst)
  
  lvst_crp <- terra::app(x = lvst, fun = function(x){
    x <- as.numeric(na.omit(x))
    if(length(x) > 1){
      y <- vegan::diversity(x = x, index = "shannon", MARGIN = 1)
    } else {
      y <- NA
    }
    return(y)
  })
  out <- paste0(root,'/data/',iso,'/livestock/lvst_diver.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = lvst_crp, out, overwrite=TRUE) }
  
  return("ok")
}
#'run
get_livestock_var(iso = iso)


#'Get livestock TLU Variable
#' @param iso (character): Country ISO-3 code
#' @return Raster file for Livestock TLU variable
#' .../livestock/lvst_tlu.tif
get_livestock_TLU_var <- function(iso){
  
  pth <- paste0(root,'/data/_global/livestock')
  grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern', SIMPLIFY = T)
  lvst <- list.files(path = pth, pattern = '*_2010_Da.tif$', full.names = T, recursive = T) %>%
    grep2(pattern = c('cattle','sheep','goats','chicken','pig','horse'), x = ., value = T) %>%
    terra::rast() %>%
    terra::crop(shpr) %>%
    terra::resample(shpr) %>%
    terra::mask(shpr)
  
  lvst_crp <- terra::app(x = lvst, fun = function(x){
    x <- as.numeric(x)
    x[which(is.na(x))] <- 0
    y <- 0.7*x[1] + 0.1*x[2] + 0.1*x[3] + 0.01*x[4] + 0.2*x[5] + 0.8*x[6]
    y <- ifelse(y == 0, NA, y)
    return(y)
  })
  out <- paste0(root,'/data/',iso,'/livestock/lvst_tlu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = lvst_crp, out, overwrite = TRUE) }
  return("ok")
  # TLU: tropical livestock units
  # TLU = Cattle * 0.7 + Sheep * 0.1 + Goat * 0.1 + Chicken * 0.01 + Pig * 0.2 + Horse * 0.8
  # Percentile: 0.1
  
  # + vacas - recursos
  # - vacas + pobreza
}
#'run
get_livestock_TLU_var(iso = iso)

#' Get Migration variables
#' @param iso (character): Country ISO-3 code
#' @returns  Raster files for migration (median, Percentil 90 P90, Coeficient of variation and Trend)
get_migr_vars <- function(iso = 'SDN'){
  
  # Load the country lowest administrative level shapefile
  if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
    dir.create(path = dirname(paste0(root,'/data/',iso,'/_shps/',iso,'.shp')), recursive = TRUE)
    shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  } else {
    shp <- raster::shapefile(x = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  }
  
  # Reference raster at 1 km
  ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
  shpr <- terra::rasterize(x = shp, y = ref, field = 'key')
  
  # ----------------------------------------------------------------- #
  # Migration: mean of percentage (time series 1970-2000)
  # ----------------------------------------------------------------- #
  # Migration
  mgrt <- list.files(path = paste0(root,'/data/_global/migration'), pattern = '*.tif$', full.names = T, recursive = T) %>% terra::rast(.)
  names(mgrt) <- paste0('yr',c(1975,1985,1995))
  mgrt <- terra::crop(x = mgrt, terra::ext(shp))
  
  # ------------------------------------ #
  # Median
  # ------------------------------------ #
  
  mmgrt <- median(mgrt)
  mmgrt_crp <- terra::resample(x = mmgrt, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/migration/medn_migration.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mmgrt_crp, out) }
  
  # ------------------------------------- #
  # 90th percentile
  # ------------------------------------- #
  
  p95mgrt <- terra::quantile(mgrt, probs = 0.90)
  vmgrt_crp <- terra::resample(x = p95mgrt, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/migration/p90_migration.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vmgrt_crp, out) }
  
  
  # ------------------------------------ #
  # Coefficient of variation
  # ------------------------------------ #
  
  vmgrt <- terra::stdev(mgrt)/mean(mgrt)
  vmgrt_crp <- terra::resample(x = vmgrt, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/migration/cvar_migration.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vmgrt_crp, out) }
  
  # Crop Trend
  out <- paste0(root,'/data/',iso,'/migration/trnd_migration.tif')
  if(!file.exists(out)){ tmgrt <- terra::app(x = mgrt, fun = function(x){
    x <- as.numeric(na.omit(x))
    if(length(x) > 1){
      y <- trend::sens.slope(x)$estimates
    } else {
      y <- NA
    }
    return(y)
  })
  tmgrt_crp <- terra::resample(x = tmgrt, y = ref) %>% terra::mask(x = ., mask = shpr)
  
  dir.create(dirname(out), showWarnings = F, recursive = T)
  terra::writeRaster(x = tmgrt_crp, out, overwrite = TRUE)
  }
  
}
#'run
get_migr_vars(iso = iso)

#'Get Pasture ares
#' @param iso (character): Country ISO-3 code
#' @return Raster File for pasture area variable
#' .../pasture_area/pasture.tif
get_pasture_vars <- function(iso){
  
  smm <- pas %>%
    terra::crop(shpr) %>%
    terra::resample(shpr) %>%
    terra::mask(shpr)
  
  out <- paste0(root,'/data/',iso,'/pasture_area/pasture.tif')
  
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = smm, out) }
}
#'run
get_pasture_vars(iso = iso)

#'Function to get Soil variables
#' check user has access to this folder '//192.168.20.97/data_cluster17/GLOBAL/Biofisico/SoilGrids250m'
#' @param iso (character): Country ISO-3 code
#' @returns Raster files for Soil clay content and Soil saturation
#' .../climatic_indexes/temp/soilcp.tif
#' ../climatic_indexes/temp/soilsat.tif
get_soil_var <- function(iso = iso){
  
  cat("Processing: ", iso, "/n")
  shp_fl <- paste0('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/', iso,'/_shps/', iso, '.shp')
  
  out_dir_soil <- paste0('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/', iso,'/climatic_indexes/temp/')
  
  if(!dir.exists(out_dir_soil)){dir.create(out_dir_soil)}
  
  soil_cp  <- paste0('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/', iso,'/climatic_indexes/temp/soilcp.tif')
  soil_sat <- paste0('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/', iso,'/climatic_indexes/temp/soilsat.tif')
  
  #' Input parameters:
  #   shp_fl: shapefile with the regions of interest
  #   root_depth: root depth in cm (it's assumed to be constant over all coordinates)
  #   outfiles: output file paths
  # Output:
  #   Raster files of soil capacity and soil saturation values
  
  get_soil <- function(shp_fl = shp_fl, root_depth = 60, outfiles = c(soil_cp,
                                                                      soil_sat)){
    
    if(sum(!file.exists(outfiles)) != 0){
      # Load packages
      if(!require(pacman)){install.packages('pacman'); library(pacman)} else {suppressMessages(library(pacman))}
      suppressMessages(pacman::p_load(terra, tidyverse))
      
      # Load CHIRPS template
      #//catalogue/BaseLineDataCluster01/observed/gridded_products/chirps/daily/chirps-v2.0.2020.01.01.tif
      tmp <- terra::rast("//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif")
      ## ROI: regions of interest
      shp <- terra::vect(shp_fl)
      r <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      r[r == -9999] <- NA
      r[!is.na(r)]  <- 1
      crd <- r %>% terra::as.data.frame(xy = T, na.rm = T)
      names(crd)[3] <- 'vals'
      crd$id <- 1:nrow(crd)
      crd$vals <- NULL
      crd <- crd[,c('id','x','y')]
      
      # Soil data repository. ISRIC soil data 250 m
      soils_root <- '//192.168.20.97/data_cluster17/GLOBAL/Biofisico/SoilGrids250m'
      # Soil organic carbon content
      orc <- terra::rast(list.files(paste0(soils_root,'/Chemical soil properties/Soil organic carbon content'), pattern = '.tif$', full.names = T) %>% sort())
      # Cation exchange capacity
      cec <- terra::rast(list.files(paste0(soils_root,'/Chemical soil properties/Cation exchange capacity (CEC)'), pattern = '.tif$', full.names = T) %>% sort())
      # Soil ph in H2O
      phx <- terra::rast(list.files(paste0(soils_root,'/Chemical soil properties/Soil ph in H2O'), pattern = '.tif$', full.names = T) %>% sort())
      # Sand content
      snd <- terra::rast(list.files(paste0(soils_root,'/Physical soil properties/Sand content'), pattern = '.tif$', full.names = T) %>% sort())
      # Silt content
      slt <- terra::rast(list.files(paste0(soils_root,'/Physical soil properties/Silt content'), pattern = '.tif$', full.names = T) %>% sort())
      # Clay content
      cly <- terra::rast(list.files(paste0(soils_root,'/Physical soil properties/Clay content (0-2 micro meter) mass fraction'), pattern = 'sl[1-7]_250m_ll.tif$', full.names = T) %>% sort())
      # Bulk density
      bld <- terra::rast(list.files(paste0(soils_root,'/Physical soil properties/Bulk density (fine earth)'), pattern = '.tif$', full.names = T) %>% sort())
      
      # Put all layers together and resampling them to the proper resolution 5 km
      soil <- terra::rast(list(orc,cec,phx,snd,slt,cly,bld))
      soil <- soil %>%
        terra::crop(., terra::ext(r)) %>%
        terra::resample(., r) %>%
        terra::mask(., mask = r)
      
      # Obtain soil data for the corresponding coordinates
      soil_data <- cbind(crd, terra::extract(soil, crd[,c('x','y')]))
      soil_data$ID <- NULL
      
      # Arrange the soil data at different depth levels
      soil_data2 <- soil_data %>%
        tidyr::pivot_longer(names_to = 'var', values_to = 'val', -(1:3)) %>%
        tidyr::separate(col = 'var', sep = '_M_', into = c('var','depth')) %>%
        tidyr::pivot_wider(names_from = 'var', values_from = 'val') %>%
        dplyr::arrange(id)
      soil_data2$depth <- gsub('_250m_ll','',soil_data2$depth)
      
      # Get Available soil water capacity per depth level
      soil_data2 <- cbind(soil_data2,AWCPTF(SNDPPT = soil_data2$SNDPPT,
                                            SLTPPT = soil_data2$SLTPPT,
                                            CLYPPT = soil_data2$CLYPPT,
                                            ORCDRC = soil_data2$ORCDRC,
                                            BLD = soil_data2$BLDFIE,
                                            CEC = soil_data2$CECSOL,
                                            PHIHOX = soil_data2$PHIHOX/10,
                                            h1=-10, h2=-20, h3=-33))
      
      #now calculate the ASW in mm for each soil horizon
      soil_data2$tetaFC <- soil_data2$WWP + soil_data2$AWCh3 #volumetric water content at field capacity (fraction)
      soil_data2$AWSat <- soil_data2$tetaS - soil_data2$tetaFC
      
      soil_data2$depth[soil_data2$depth == "sl1"] <- 0
      soil_data2$depth[soil_data2$depth == "sl2"] <- 5
      soil_data2$depth[soil_data2$depth == "sl3"] <- 15
      soil_data2$depth[soil_data2$depth == "sl4"] <- 30
      soil_data2$depth[soil_data2$depth == "sl5"] <- 60
      soil_data2$depth[soil_data2$depth == "sl6"] <- 100
      soil_data2$depth[soil_data2$depth == "sl7"] <- 200
      soil_data2$depth <- as.numeric(soil_data2$depth)
      
      soilcap_calc <- function(x, y, rdepth=60, minval, maxval) {
        if (length(x) != length(y)) {stop("length of x and y must be the same")}
        rdepth <- max(c(rdepth,minval)) #cross check
        rdepth <- min(c(rdepth,maxval)) #cross-check
        wc_df <- data.frame(depth=y,wc=x)
        if (!rdepth %in% wc_df$depth) {
          wc_df1 <- wc_df[which(wc_df$depth < rdepth),]
          wc_df2 <- wc_df[which(wc_df$depth > rdepth),]
          y1 <- wc_df1$wc[nrow(wc_df1)]; y2 <- wc_df2$wc[1]
          x1 <- wc_df1$depth[nrow(wc_df1)]; x2 <- wc_df2$depth[1]
          ya <- (rdepth-x1) / (x2-x1) * (y2-y1) + y1
          wc_df <- rbind(wc_df1,data.frame(depth=rdepth,wc=ya),wc_df2)
        }
        wc_df <- wc_df[which(wc_df$depth <= rdepth),]
        wc_df$soilthick <- wc_df$depth - c(0,wc_df$depth[1:(nrow(wc_df)-1)])
        wc_df$soilcap <- wc_df$soilthick * wc_df$wc
        soilcp <- sum(wc_df$soilcap) * 10 #in mm
        return(soilcp)
      }
      
      soil_data4 <- soil_data2 %>%
        dplyr::group_by(id) %>%
        dplyr::group_split(id) %>%
        purrr::map(.f = function(px){
          scp  <- soilcap_calc(x=px$AWCh3, y=px$depth, rdepth = root_depth, minval=45, maxval=100)
          ssat <- soilcap_calc(x=px$AWSat, y=px$depth, rdepth = root_depth, minval=45, maxval=100)
          df <- data.frame(id = unique(px$id),
                           x  = unique(px$x),
                           y  = unique(px$y),
                           scp  = scp,
                           ssat = ssat)
          return(df)
        }) %>%
        dplyr::bind_rows()
      
      scp  <- terra::rast(x = as.matrix(soil_data4[,c('x','y','scp')]), type = 'xyz', crs = terra::crs(r))
      ssat <- terra::rast(x = as.matrix(soil_data4[,c('x','y','ssat')]), type = 'xyz', crs = terra::crs(r))
      
      dir.create(path = dirname(outfiles[1]), FALSE, TRUE)
      
      terra::writeRaster(x = scp, filename = outfiles[1], overwrite = T)
      terra::writeRaster(x = ssat, filename = outfiles[2], overwrite = T)
      
    } else {
      cat('Soil capacity and soil saturation variables are already calculated.\n')
    }
    return(cat('Get soil data: finished successfully!\n'))
  }
  
  
  
  
  
  if(file.exists(soil_cp) & file.exists(soil_sat)){
    print("file already exists")
  }else{
    get_soil(shp_fl = shp_fl, root_depth = 60, outfiles = c(soil_cp,
                                                            soil_sat))
    
  }
  
  return(NULL)  
}
#'run
get_soil_var(iso = iso)

#'Function to get Years of Education variables for Male, Female and the Difference
#' @param iso (character): Country ISO-3 code
#' @returns Raster files for education variables (median, coefficient of variation and Trend)
#' .../education/medn_male_edu.tif
#' .../education/cvar_male_edu.tif
#' .../education/trnd_male_edu.tif
#' .../education/medn_female_edu.tif
#' .../education/cvar_female_edu.tif
#' .../education/trnd_female_edu.tif
#' .../education/medn_difference_edu.tif
#' .../education/cvar_difference_edu.tif
#' .../education/trnd_difference_edu.tif
get_edu_vars <- function(iso = 'PHL'){
  
  # Load the country lowest administrative level shapefile
  if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
    dir.create(path = dirname(paste0(root,'/data/',iso,'/_shps/',iso,'.shp')), recursive = TRUE)
    shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  } else {
    shp <- raster::shapefile(x = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
    shp <- as(shp, 'SpatVector')
  }
  
  # Reference raster at 1 km
  ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
  shpr <- terra::rasterize(x = shp, y = ref, field = 'key')
  
  # ------------------------------------------------------------- #
  # Education: mean of years of education (time series 2000-2017)
  # ------------------------------------------------------------- #
  # By male
  mle <- terra::rast(paste0(root,'/data/_global/education/IHME_LMIC_EDU_2000_2017_MEAN_15_49_MALE_MEAN_Y2019M12D24.TIF'))
  # By female
  fle <- terra::rast(paste0(root,'/data/_global/education/IHME_LMIC_EDU_2000_2017_MEAN_15_49_FEMALE_MEAN_Y2019M12D24.TIF'))
  
  # Difference of years of education
  if(!file.exists(paste0(root,'/data/_global/education/difference.TIF'))){
    dff <- mle - fle
    terra::writeRaster(dff, paste0(root,'/data/_global/education/difference.TIF'))
  } else {
    dff <- terra::rast(paste0(root,'/data/_global/education/difference.TIF'))
  }
  
  # Assign years names to each layer
  names(mle) <- names(fle) <- names(dff) <- paste0('yr',2000:2017)
  
  # ------------------------------------ #
  # Median
  # ------------------------------------ #
  
  mmle <- median(mle)
  mfle <- median(fle)
  mdff <- median(dff)
  
  # ------------------------------------ #
  # Coefficient of variation
  # ------------------------------------ #
  
  vmle <- terra::stdev(mle)/mean(mle)
  vfle <- terra::stdev(fle)/mean(fle)
  vdff <- terra::stdev(dff)/mean(dff)
  
  # Crop Median
  mmle_crp <- terra::crop(x = mmle, terra::ext(shp))
  mmle_crp <- terra::resample(x = mmle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/education/medn_male_edu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mmle_crp, out) }
  
  mfle_crp <- terra::crop(x = mfle, terra::ext(shp))
  mfle_crp <- terra::resample(x = mfle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/education/medn_female_edu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mfle_crp, out) }
  
  mdff_crp <- terra::crop(x = mdff, terra::ext(shp))
  mdff_crp <- terra::resample(x = mdff_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/education/medn_difference_edu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = mdff_crp, out) }
  
  # Crop Coefficient of variation
  vmle_crp <- terra::crop(x = vmle, terra::ext(shp))
  vmle_crp <- terra::resample(x = vmle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/education/cvar_male_edu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vmle_crp, out) }
  
  vfle_crp <- terra::crop(x = vfle, terra::ext(shp))
  vfle_crp <- terra::resample(x = vfle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/education/cvar_female_edu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vfle_crp, out) }
  
  vdff_crp <- terra::crop(x = vdff, terra::ext(shp))
  vdff_crp <- terra::resample(x = vdff_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/education/cvar_difference_edu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = vdff_crp, out) }
  
  # Crop Trend
  mle_crp <- terra::crop(x = mle, terra::ext(shp))
  tmle_crp <- terra::app(x = mle_crp, fun = function(x){
    x <- as.numeric(na.omit(x))
    if(length(x) > 1){
      y <- trend::sens.slope(x)$estimates
    } else {
      y <- NA
    }
    return(y)
  })
  tmle_crp <- terra::resample(x = tmle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/education/trnd_male_edu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = tmle_crp, out) }
  
  fle_crp <- terra::crop(x = fle, terra::ext(shp))
  tfle_crp <- terra::app(x = fle_crp, fun = function(x){
    x <- as.numeric(na.omit(x))
    if(length(x) > 1){
      y <- trend::sens.slope(x)$estimates
    } else {
      y <- NA
    }
    return(y)
  })
  tfle_crp <- terra::resample(x = tfle_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/education/trnd_female_edu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = tfle_crp, out) }
  
  dff_crp <- terra::crop(x = dff, terra::ext(shp))
  tdff_crp <- terra::app(x = dff_crp, fun = function(x){
    x <- as.numeric(na.omit(x))
    if(length(x) > 1){
      y <- trend::sens.slope(x)$estimates
    } else {
      y <- NA
    }
    return(y)
  })
  tdff_crp <- terra::resample(x = tdff_crp, y = ref) %>% terra::mask(x = ., mask = shpr)
  out <- paste0(root,'/data/',iso,'/education/trnd_difference_edu.tif')
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = tdff_crp, out) }
  
  return(cat('Done\n'))
  
}
#'run
get_edu_vars(iso = iso)


#'Function to get agricultural area variable
#' @param iso (character): Country ISO-3 code
#' @return Raster file for agricultural area
#' .../agricultural_area/ag_area.tif
get_agricultural_area_var <- function(iso){
  shp <- terra::vect(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
  ref  <- terra::rast(extent = terra::ext(shp), crs = terra::crs(shp), resolution = c(0.008333334, 0.008333334))
  shpr <- terra::rasterize(x = shp, y = ref, field = 'NAME_0')
  
  smm <- smm %>%
    terra::crop(shpr) %>%
    terra::resample(shpr) %>%
    terra::mask(shpr)
  
  out <- paste0(root,'/data/',iso,'/agricultural_area/ag_area.tif')
  
  dir.create(dirname(out), showWarnings = F, recursive = T)
  if(!file.exists(out)){ terra::writeRaster(x = smm, out, overwrite=TRUE) }
}
#'run
get_agricultural_area_var(iso = iso)

#' Function to get Absolute Wealth index
#' @param iso (character): Conuntry ISO-3 code
#' @return Raster file for Absolute wealth index
#' .../wealth_index/{iso}_AWE.tif
get_AWI_var <- function(iso, shp_country){
  ISO3 <- iso
  gdp  <- readr::read_csv(paste0(root,'/data/_global/wealth_index/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488.csv'))
  gini <- readr::read_csv(paste0(root,'/data/_global/wealth_index/API_SI.POV.GINI_DS2_en_csv_v2_2445276/API_SI.POV.GINI_DS2_en_csv_v2_2445276.csv'))
  mask <- raster::raster(paste0(root,'/data/_global/masks/mask_world_1km.tif'))
  
  wealth_dir <- paste0(root,'/data/_global/wealth_index')
  rwi_out_dir <-  paste0(root,'/data/', iso,'/wealth_index')
  
  if(!dir.exists(rwi_out_dir)){dir.create(rwi_out_dir, F, T)}
  
  country_gini <- gini %>% 
    dplyr::filter(`Country Code` == ISO3) %>% 
    dplyr::select(where(~!all(is.na(.)))) %>% 
    dplyr::pull(ncol(.))
  country_gini <- country_gini/100
  
  country_gdp <- gdp %>% 
    filter(`Country Code` == ISO3) %>% 
    dplyr::pull(`2009`)
  
  shp_country <- as(shp, "Spatial")#raster::shapefile(paste0(root,'/data/',ISO3,'/_shps/',ISO3,'.shp'))
  
  country_rwi <- read_csv(paste0(wealth_dir,'/',ISO3,'_relative_wealth_index.csv')) %>% 
    dplyr::select(rwi, longitude, latitude)
  
  country_rwi$rank <- rank(country_rwi$rwi, ties.method = "random")/(max(rank(country_rwi$rwi, ties.method = "random"))+1)
  country_rwi$icdf <- ICDF(x     = country_rwi$rank, 
                           theta = 1, 
                           gini  = country_gini, 
                           gdppc = country_gdp) # 41986 - OCED mean
  
  country_rwi$AWE <- country_rwi$icdf*(country_gdp/(mean(country_rwi$icdf)) )
  cap <- quantile(country_rwi$AWE, probs = 0.98)
  country_rwi <- country_rwi %>% 
    dplyr::mutate(AWE = ifelse(AWE > cap, cap, AWE))
  
  coordinates(country_rwi) <-  ~longitude+latitude
  crs(country_rwi) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  r <- raster(resolution = 0.04166667, ext = extent(shp_country))
  crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
  r <- raster::mask(r, shp_country)
  r_f <- rasterize(country_rwi, r, fun = mean, field = "AWE")
  r_f <- raster::resample(r_f, mask %>% raster::crop(., extent(shp_country)) )
  
  raster::writeRaster(r_f, paste0(rwi_out_dir,"/",ISO3,"_AWE.tif"), overwrite = T)
  
  
}
#' run
get_AWI_var(iso = iso, shp_country = shp)


#' Function to get Relative Wealth index
#' @param iso (character): Conuntry ISO-3 code
#' @return Raster file for Relative wealth index
#' .../wealth_index/{iso}_rwi.tif
get_RWI_var <- function(iso, shp_country){
  ISO3 <- iso
  wealth_dir  <- paste0(root,'/data/_global/wealth_index')
  rwi_out_dir <- paste0(root,'/data/',ISO3,'/wealth_index' )
  
  if(!dir.exists(rwi_out_dir)){dir.create(rwi_out_dir,F,T)}
  
  wealth_df   <- readr::read_csv(list.files(wealth_dir, pattern = paste0("^",ISO3), full.names = T))
  mask        <- raster::raster(paste0(root,'/data/_global/masks/mask_world_1km.tif'))
  shp_country <- as(shp_country, "Spatial")#raster::shapefile(paste0(root,'/data/',ISO3,'/_shps/',ISO3,'.shp'))
  
  coordinates(wealth_df) <- ~longitude+latitude
  crs(wealth_df) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  r <- raster(resolution = 0.04166667, ext = extent(shp_country))
  crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
  r <- raster::mask(r, shp_country)
  r_f <- rasterize(wealth_df, r, fun = mean, field = "rwi")
  r_f <- raster::resample(r_f, mask %>% raster::crop(., extent(shp_country)) )
  
  raster::writeRaster(r_f, paste0(rwi_out_dir,"/",ISO3,"_rwi.tif"), overwrite = T)
  
  return("ok")
}
#'Run
get_RWI_var(iso = iso, shp_country = shp)


