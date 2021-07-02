# CPC temperature: regression downscaling using SRTM 5 km
# By: H. Achicanoy
# CIAT, 2021

options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,raster,sf,RSAGA,parallel,doSNOW))

cpc_downscaling <- function(iso = 'SDN', country = 'Sudan'){
  # Load country shapefile
  shp <- raster::getData(name = 'GADM', country = iso, level = 0, download = T)
  
  grep2 <- Vectorize(grep, vectorize.args = 'pattern')
  climate <- list.files(path = paste0(root,'/cpc_data/50km/stack'), full.names = T) %>%
    grep2(pattern = 1981:2020, x = ., value = T) %>%
    as.vector()
  clnames <- basename(climate)
  
  trmv  <- list.files(path = paste0(root,'/cpc_data/50km/individuals'), all.files = T, full.names = T)
  trmv %>% purrr::map(file.remove)
  
  # Parallelization
  clusterExport <- local({
    gets <- function(n, v) { assign(n, v, envir = .GlobalEnv); NULL }
    function(cl, list, envir = .GlobalEnv) {
      ## do this with only one clusterCall--loop on slaves?
      for (name in list) {
        clusterCall(cl, gets, name, get(name, envir = envir))
      }
    }
  })
  createCluster <- function(noCores, logfile = "/dev/null", export = NULL, lib = NULL) {
    require(doSNOW)
    cl <- makeCluster(noCores, type = "SOCK", outfile = logfile)
    if(!is.null(export)) clusterExport(cl, export)
    if(!is.null(lib)) {
      plyr::l_ply(lib, function(dum) { 
        clusterExport(cl, "dum", envir = environment())
        clusterEvalQ(cl, library(dum, character.only = TRUE))
      })
    }
    registerDoSNOW(cl)
    return(cl)
  }
  env <- rsaga.env(path = 'C:/sagaGIS')
  
  cl  <- createCluster(10, export = list("climate","clnames","shp","root","env","country"), lib = list("tidyverse","raster","sf","RSAGA"))
  1:length(climate) %>% parallel::parLapply(cl, ., function(l){
    
    # Load raster stack
    rstck <- raster::stack(climate[l])
    # Rotate raster stack
    if(!raster::rotated(rstck)){
      rstck <- raster::rotate(rstck)
    }
    # Cropping to Vietnam and Indonesia region
    rstck <- raster::crop(rstck, raster::extent(shp))
    # Individualize rasters
    rstck <- raster::unstack(rstck)
    days  <- rstck %>% purrr::map(names) %>% unlist %>% substr(., start = 1, stop = 11) %>% gsub('X','',.)
    
    1:length(rstck) %>%
      purrr::map(., .f = function(i){raster::writeRaster(rstck[[i]], paste0(root,'/cpc_data/50km/individuals/',gsub('.nc','',clnames[l]) %>% substr(.,1,4),'.',days[i],'.tif'), overwrite = T)})
    
    1:length(rstck) %>%
      purrr::map(., .f = function(i){
        
        if(!dir.exists(dirname(paste0(root,'/cpc_data/5km/',tolower(country),'/',gsub('.nc','',clnames[l]) %>% substr(.,1,4),'.',days[i],'.tif')))){
          dir.create(path = dirname(paste0(root,'/cpc_data/5km/',tolower(country),'/',gsub('.nc','',clnames[l]) %>% substr(.,1,4),'.',days[i],'.tif')),
                     FALSE,
                     recursive = TRUE)
        }
        
        rsaga.geoprocessor(lib    = 'statistics_regression',
                           module = 'GWR for Grid Downscaling',
                           param  = list(PREDICTORS = paste0(root,'/cpc_data/srtm/',country,'/srtm_5km.tif'),
                                         REGRESSION = paste0(root,'/cpc_data/5km/',tolower(country),'/',gsub('.nc','',clnames[l]) %>% substr(.,1,4),'.',days[i],'.tif'),
                                         DEPENDENT  = paste0(root,'/cpc_data/50km/individuals/',gsub('.nc','',clnames[l]) %>% substr(.,1,4),'.',days[i],'.tif')),
                           env    = env)
        
      })
    
  })
  parallel::stopCluster(cl)
}
