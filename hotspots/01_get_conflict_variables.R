options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

# -------------------------------------- #
# Obtain conflict variables
# -------------------------------------- #

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA')
cnty <- c('Sudan','Zimbabwe','Senegal','Mali','Nigeria','Kenya','Uganda')

get_conflict_vars <- function(iso = 'KEN', country = 'Kenya'){
  
  # Create conflict file for the specific country
  out <- paste0(root,'/data/',iso,'/conflict/',iso,'_conflict.csv')
  dir.create(path = dirname(out), F, T)
  if(!file.exists(out)){
    # Filter African conflict to the specific country
    cnf <- readxl::read_excel(paste0(root,'/data/_global/conflict/Africa_1997-2021_Apr02.xlsx'), sheet = 1)
    cnf <- cnf %>% dplyr::filter(COUNTRY == country)
    readr::write_csv(cnf, out)
    conflict <- cnf; rm(cnf, out)
  } else {
    # Load country conflict
    conflict <- readr::read_csv(out); rm(out)
  }; rm(out)
  
  # Load the country lowest administrative level shapefile
  if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
    dir.create(path = dirname(paste0(root,'/data/',iso,'/_shps/',iso,'.shp')), recursive = TRUE)
    shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
  } else {
    shp <- raster::shapefile(x = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
  }
  sft <- as(shp, 'SpatVector') # Shapefile in terra format
  
  # Extract administrative names using reported conflict coordinates by ACLED
  conflict  <- cbind(conflict,raster::extract(x = shp, y = conflict[,c('LONGITUDE','LATITUDE')]))
  
  # Get conflict summarization variables per coordinates
  cnf_summ <- conflict %>%
    dplyr::group_by(LONGITUDE, LATITUDE) %>%
    dplyr::summarise(EVENTS           = dplyr::n(),
                     TYPE_RICHNESS    = EVENT_TYPE %>% unique() %>% length(),
                     SUBTYPE_RICHNESS = SUB_EVENT_TYPE %>% unique() %>% length(),
                     ACTOR1_RICHNESS  = ACTOR1 %>% unique() %>% length(),
                     ACTOR2_RICHNESS  = ACTOR2 %>% unique() %>% length(),
                     FATALITIES       = sum(FATALITIES)) %>%
    dplyr::ungroup()
  
  # Create reference rasters in terra and raster formats
  ref <- terra::rast(extent = terra::ext(sft), crs = terra::crs(sft), resolution = c(0.008333334, 0.008333334)) # Terra
  t.empty <- terra::rast(ref)
  sfr <- terra::rasterize(x = sft, y = ref, field = 'key')
  
  # Get cellID for ACLED reported coordinates
  cnf_summ$cellID <- terra::cellFromXY(object = ref, xy = base::as.matrix(cnf_summ[,c('LONGITUDE','LATITUDE')]))
  
  # Get conflict summarization variables per coordinates
  # summing up the duplicated coordinates where close conflict events occur
  cnf_summ <- cnf_summ %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cellID) %>%
    dplyr::summarise(EVENTS = sum(EVENTS),
                     TYPE_RICHNESS = max(TYPE_RICHNESS),
                     SUBTYPE_RICHNESS = max(SUBTYPE_RICHNESS),
                     ACTOR1_RICHNESS = max(ACTOR1_RICHNESS),
                     ACTOR2_RICHNESS = max(ACTOR2_RICHNESS),
                     FATALITIES = sum(FATALITIES)) %>%
    dplyr::ungroup()
  # Add central coordinates from reference raster
  cnf_summ <- cbind(base::as.data.frame(terra::xyFromCell(object = ref, cell = cnf_summ$cellID)), cnf_summ)
  cnf_summ <- cnf_summ %>% tidyr::drop_na()
  
  # Add pseudo non-conflict reported coordinates
  # 1 . Create density raster using all conflict points
  roi     <- spatstat.geom::as.owin(sf::st_bbox(sf::st_as_sf(shp))) # region of interest
  cnf_crd <- sp::SpatialPointsDataFrame(coords = cnf_summ[,c('x','y')], proj4string = raster::crs(shp), data = cnf_summ)
  pts     <- sp::coordinates(cnf_crd)
  p       <- spatstat.geom::ppp(pts[,1], pts[,2], window = roi); rm(roi)
  ds      <- spatstat.core::density.ppp(p); rm(pts,p)
  d       <- raster::raster(ds); rm(ds)
  d       <- as(d, 'SpatRaster')
  terra::crs(d) <- terra::crs(sfr)
  d       <- terra::resample(x = d, y = ref) %>% terra::mask(mask = sfr)
  terra::values(d) <- terra::values(d)/max(terra::values(d), na.rm = T)
  d <- 1-d # Give more importance to pixels without conflict presence
  
  # Create auxiliary data.frame
  aux.df <- d %>% terra::as.data.frame(xy = T, na.rm = T)
  aux.df$cellID <- terra::cellFromXY(d, xy = as.matrix(aux.df[,c('x','y')]))
  
  # Create a sample of 1000 pixels
  set.seed(1235); smp <- sample(x = aux.df$cellID, size = 1000, replace = F, prob = aux.df[,3])
  # Remove possible overlaps between conflict and pseudo non-conflict pixels
  smp <- base::setdiff(smp, cnf_summ$cellID)
  # Create pseudo-coords with no conflict
  psd <- aux.df[aux.df$cellID %in% smp,c('x','y','cellID')] %>%
    dplyr::mutate(EVENTS           = 0,
                  TYPE_RICHNESS    = 0,
                  SUBTYPE_RICHNESS = 0,
                  ACTOR1_RICHNESS  = 0,
                  ACTOR2_RICHNESS  = 0,
                  FATALITIES       = 0)
  cnf_summ <- dplyr::bind_rows(cnf_summ, psd); rm(psd, aux.df, d, smp)
  
  # Interpolation with Random Forest Spatial Interpolation (RFSI)
  # https://github.com/AleksandarSekulic/RFSI
  rfsi.interpolation <- function(df = cnf_summ, vr = 'TYPE_RICHNESS', ref = ref){
    
    # Extract predictors information in coordinates
    prd <- c(list.files(path = paste0(root,'/data/',iso,'/accessibility'), pattern = '*.tif$', full.names = T),
             list.files(path = paste0(root,'/data/',iso,'/nightlights'), pattern = '*.tif$', full.names = T)) %>%
      purrr::map(.f = function(f){
        r <- terra::rast(f)
        vct <- terra::extract(x = r, y = df[,c('x','y')])[,2]
        df  <- data.frame(vct); names(df) <- names(r)
        return(df)
      }) %>%
      purrr::reduce(.x = ., .f = cbind)
    names(prd) <- c('accessibility','friction','cvar_lights','medn_lights','trnd_lights')
    
    df <- cbind(df, prd); rm(prd)
    
    # Training sample
    # 1. Create a small data frame with needed variables
    dt <- df %>%
      dplyr::mutate(time = 0) %>%
      dplyr::select('cellID','x','y','time','accessibility','friction','cvar_lights','medn_lights','trnd_lights',vr) %>%
      base::as.data.frame()
    rownames(dt) <- 1:nrow(dt)
    colnames(dt)[ncol(dt)] <- 'Layer'
    # 2. Transform to Spatial Pixels Data Frame format
    sdt <- sp::SpatialPixelsDataFrame(points = as.matrix(dt[,c('x','y')]), data = dt, proj4string = raster::crs(ref %>% raster::raster()))
    
    # 3. Hyper-parameter grid to optimize
    n.obs           <- seq(2,10,2)
    min.node.size   <- seq(2,10,2)
    sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
    splitrule       <- 'variance'
    ntree           <- 1000
    mtry            <- 3:(2+2*max(n.obs))
    tgrid <- expand.grid(min.node.size   = min.node.size,
                         num.trees       = ntree,
                         mtry            = mtry,
                         n.obs           = n.obs,
                         sample.fraction = sample.fraction)
    
    # Define formula
    fml <- as.formula('Layer ~ accessibility + friction + cvar_lights + medn_lights + trnd_lights')
    # 4. RFSI Hyper-parameter tuning
    rfsi_fit <- meteo::tune.rfsi(formula    = fml,
                                 data       = sdt,
                                 zero.tol   = 0,
                                 use.idw    = FALSE,
                                 s.crs      = terra::crs(ref),
                                 t.crs      = terra::crs(ref),
                                 tgrid      = tgrid,
                                 tgrid.n    = length(tgrid),
                                 tune.type  = 'LLO',
                                 k          = 5,
                                 seed       = 1235,
                                 cpus       = detectCores()-1,
                                 progress   = TRUE,
                                 acc.metric = 'RMSE',
                                 importance = 'impurity')
    
    # Predict pixels
    rax <- ref
    terra::values(rax) <- 1:(terra::ncell(rax))
    rax <- rax %>% terra::mask(sfr)
    crd_rast <- rax %>% terra::as.data.frame(xy = T) %>% .[,1:2]
    crd_rast$cellID <- terra::cellFromXY(object = rax, xy = as.matrix(crd_rast)); rm(rax)
    crd_rast <- crd_rast[crd_rast$cellID %in% base::setdiff(crd_rast$cellID, df$cellID[df$EVENTS > 0]),]
    
    # Extract predictors information in coordinates
    prd <- c(list.files(path = paste0(root,'/data/',iso,'/accessibility'), pattern = '*.tif$', full.names = T),
             list.files(path = paste0(root,'/data/',iso,'/nightlights'), pattern = '*.tif$', full.names = T)) %>%
      purrr::map(.f = function(f){
        r <- terra::rast(f)
        vct <- terra::extract(x = r, y = crd_rast[,c('x','y')])[,2]
        df  <- data.frame(vct); names(df) <- names(r)
        return(df)
      }) %>%
      purrr::reduce(.x = ., .f = cbind)
    names(prd) <- c('accessibility','friction','cvar_lights','medn_lights','trnd_lights')
    
    crd_rast <- cbind(crd_rast, prd); rm(prd)
    crd_rast <- crd_rast %>% tidyr::drop_na()
    
    ndt <- crd_rast %>% dplyr::mutate(time = 0) %>% dplyr::select(cellID,x,y,time,accessibility,friction,cvar_lights,medn_lights,trnd_lights) %>% base::as.data.frame()
    rownames(ndt) <- 1:nrow(ndt)
    sndt <- sp::SpatialPixelsDataFrame(points = as.matrix(ndt[,c('x','y')]), data = ndt, proj4string = raster::crs(ref %>% raster::raster()))
    
    rfsi_pred <- meteo::pred.rfsi(model         = rfsi_fit$final.model,
                                  data          = sdt,
                                  zcol          = 'Layer',
                                  data.staid.x.y.time = 1:4,
                                  newdata       = sndt,
                                  newdata.staid.x.y.time = 1:4,
                                  output.format = "SpatialPixelsDataFrame",
                                  zero.tol      = 0,
                                  s.crs         = terra::crs(ref),
                                  newdata.s.crs = terra::crs(ref),
                                  t.crs         = terra::crs(ref),
                                  cpus          = detectCores()-1,
                                  progress      = TRUE)
    
    tbl1 <- df %>% dplyr::filter(EVENTS > 0) %>% dplyr::select('cellID','x','y',vr)
    names(tbl1)[ncol(tbl1)] <- 'Layer'
    tbl2 <- base::as.data.frame(rfsi_pred)
    tbl2 <- tbl2[,c('x','y','pred')]
    names(tbl2)[ncol(tbl2)] <- 'Layer'
    tbl <- dplyr::bind_rows(tbl1, tbl2); rm(tbl1, tbl2)
    
    rfsi.rst <- raster::rasterFromXYZ(xyz = tbl[,c('x','y','Layer')], res = raster::res(raster::raster(ref)), crs = raster::crs(raster::raster(ref)))
    # rfsi.rst < rfsi.rst %>% raster::mask(x = ., mask = rs)
    raster::writeRaster(x = rfsi.rst, filename = paste0(root,'/data/',iso,'/conflict/',vr,'_rfsi.tif'), overwrite = TRUE)
    return(rfsi.rst)
    
  }
  rfsi.interpolation(df = cnf_summ, vr = 'EVENTS',           ref = ref)
  rfsi.interpolation(df = cnf_summ, vr = 'TYPE_RICHNESS',    ref = ref)
  rfsi.interpolation(df = cnf_summ, vr = 'SUBTYPE_RICHNESS', ref = ref)
  rfsi.interpolation(df = cnf_summ, vr = 'ACTOR1_RICHNESS',  ref = ref)
  rfsi.interpolation(df = cnf_summ, vr = 'ACTOR2_RICHNESS',  ref = ref)
  rfsi.interpolation(df = cnf_summ, vr = 'FATALITIES',       ref = ref)
  
  return(cat('Done\n'))
  
}
purrr::map2(.x = isos, .y = cnty, .f = function(iso, country){
  get_conflict_vars(iso = iso, country = country)
})
