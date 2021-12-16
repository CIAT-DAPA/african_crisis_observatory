options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend, meteo, parallel, spatstat, maptools))
#install.packages("meteo", repos="http://R-Forge.R-project.org")


root <- 'D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory' 
#'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

# -------------------------------------- #
# Obtain conflict variables
# -------------------------------------- #

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA')
cnty <- c('Sudan','Zimbabwe','Senegal','Mali','Nigeria','Kenya','Uganda')

iso <- "ZWE"


shp <- raster::shapefile(paste0(root,"/data/", iso, "/_shps/",iso,".shp" )) 

world_mask <- raster::raster(paste0(root,"/data/_global/masks/mask_world_1km.tif")) %>% 
  raster::crop(., extent(shp))

################################################################################
################################################################################
################################################################################

get_conflic_data <- function(root, iso, country = 'Senegal'){
  
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
  sft <- shp # Shapefile in terra format
  
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
    dplyr::ungroup() %>% 
    dplyr::select(x= LONGITUDE, y = LATITUDE, everything(.))
  
  
  return(list(cnf_summ, sft))
}

rfsi.interpolation <- function(root, iso , df = cnf_summ, vr = 'TYPE_RICHNESS', ref = ref, crop_layer = NULL, k = NA){
  
  # Extract predictors information in coordinates
  prd <- c(list.files(path = paste0(root,'/data/',iso,'/accessibility'), pattern = '*.tif$', full.names = T),
           list.files(path = paste0(root,'/data/',iso,'/nightlights'), pattern = '*.tif$', full.names = T)) %>%
    purrr::map(.f = function(f){
      r <- terra::rast(f)
      vct <- terra::extract(x = r, y = df[,c('x','y')])[,1]
      df  <- data.frame(vct); names(df) <- names(r)
      return(df)
    }) %>%
    purrr::reduce(.x = ., .f = cbind)
  names(prd) <- c('accessibility','friction','cvar_lights','medn_lights','trnd_lights')
  
  df <- cbind(df, prd)
  
  # Training sample
  # 1. Create a small data frame with needed variables
  dt <- df %>%
    dplyr::mutate(time = 0) %>%
    dplyr::select('cellID','x','y','time','accessibility','friction','cvar_lights','medn_lights','trnd_lights',vr) %>%
    base::as.data.frame()
  rownames(dt) <- 1:nrow(dt)
  colnames(dt)[ncol(dt)] <- 'Layer'
  # 2. Transform to Spatial Pixels Data Frame format
  sdt <- sp::SpatialPixelsDataFrame(points = as.matrix(dt[,c('x','y')]), 
                                    data = dt, 
                                    proj4string = raster::crs(ref %>% raster::raster()),
                                    tolerance =  0.916421)
  
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
                               data       = dt%>% dplyr::select(cellID,x, y,  everything(.)),
                               zero.tol   = 0,
                               use.idw    = FALSE,
                               s.crs      = CRS("+proj=longlat +datum=WGS84 +no_def"),
                               t.crs      = CRS("+proj=longlat +datum=WGS84 +no_def"),
                               tgrid      = tgrid,
                               tgrid.n    = length(tgrid),
                               tune.type  = 'LLO',
                               k          = 5,
                               seed       = 1235,
                               cpus       = parallel::detectCores()-1,
                               progress   = TRUE,
                               acc.metric = 'RMSE',
                               importance = 'impurity')
  
  metrics <- base::data.frame(rfsi_fit$tuned.parameters)
  rownames(metrics) <- 1:nrow(metrics)
  metrics$R.squared <- rfsi_fit$final.model$r.squared
  metrics$Variable  <- vr
  metrics$iso       <- iso
  metrics <- metrics %>% dplyr::select(iso,Variable,dplyr::everything(.))
  metrics$k <- k
  
  ###########
  # rfsi_cvf <- meteo::cv.rfsi(formula    = fml,
  #                            data       = sdt,
  #                            zero.tol   = 0,
  #                            use.idw    = FALSE,
  #                            s.crs      = terra::crs(ref),
  #                            t.crs      = terra::crs(ref),
  #                            tgrid      = tgrid,
  #                            tgrid.n    = length(tgrid),
  #                            tune.type  = "LLO",
  #                            seed       = 1235,
  #                            k          = 5,
  #                            acc.metric = 'RMSE',
  #                            cpus       = detectCores()-1,
  #                            progress   = TRUE)
  ###########
  
  # Predict pixels
  rax <- ref
  terra::values(rax) <- 1:(terra::ncell(rax))
  rax <- rax %>% terra::mask(roiCtr)
  crd_rast <- rax %>% terra::as.data.frame(xy = T) %>% .[,1:2]
  crd_rast$cellID <- terra::cellFromXY(object = rax, xy = as.matrix(crd_rast)); rm(rax)
  crd_rast <- crd_rast[crd_rast$cellID %in% base::setdiff(crd_rast$cellID, df$cellID[df$EVENTS > 0]),]
  
  # Extract predictors information in coordinates
  prd <- c(list.files(path = paste0(root,'/data/',iso,'/accessibility'), pattern = '*.tif$', full.names = T),
           list.files(path = paste0(root,'/data/',iso,'/nightlights'), pattern = '*.tif$', full.names = T)) %>%
    purrr::map(.f = function(f){
      r <- terra::rast(f)
      vct <- terra::extract(x = r, y = crd_rast[,c('x','y')])
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
                                data          = dt%>% dplyr::select(cellID,x, y,  everything(.)),
                                zcol          = 'Layer',
                                data.staid.x.y.time = 1:4,
                                newdata       = ndt,
                                newdata.staid.x.y.time = 1:4,
                                output.format = "SpatialPixelsDataFrame",
                                zero.tol      = 0,
                                s.crs         = CRS("+proj=longlat +datum=WGS84 +no_def"),
                                newdata.s.crs = CRS("+proj=longlat +datum=WGS84 +no_def"),
                                t.crs         = CRS("+proj=longlat +datum=WGS84 +no_def"),
                                cpus          = detectCores()-1,
                                progress      = TRUE)
  
  tbl1 <- df %>% dplyr::filter(EVENTS > 0) %>% dplyr::select('cellID','x','y',vr)
  names(tbl1)[ncol(tbl1)] <- 'Layer'
  tbl2 <- base::as.data.frame(rfsi_pred)
  tbl2 <- tbl2[,c('x','y','pred')]
  names(tbl2)[ncol(tbl2)] <- 'Layer'
  tbl <- dplyr::bind_rows(tbl1, tbl2) %>% 
    dplyr::mutate(Layer = round(Layer)); rm(tbl1, tbl2)
  
  
  #rfsi.rst <- raster::rasterFromXYZ(xyz = tbl[,c('x','y','Layer')], res = raster::res(raster::raster(ref)), crs = raster::crs(raster::raster(ref)))
  rfsi.rst <- raster::rasterize(x = tbl[,c('x','y')], y = raster::raster(ref), field  = tbl$Layer)
  
  #writeRaster(rfsi.rst, "D:/OneDrive - CGIAR/Attachments/Desktop/mapas/rsfi1.tif")
  # raster::writeRaster(x = rfsi.rst, filename = paste0(root,'/data/',iso,'/conflict/',vr,'_rfsi.tif'), overwrite = TRUE)
  return(list(metrics = metrics, raster = rfsi.rst))
  
}

conflict <- get_conflic_data(root = root,
                             iso = iso)


cnf_summ_raw <- conflict[[1]]
  
ref <- raster::rasterize(conflict[[2]], world_mask)
ref[!is.na(ref[])]<- 1
plot(ref)
## Define conflict influence area

# 1. Calculate distance matrix among conflict reported points
crd_terra <- terra::vect(as.matrix(cnf_summ_raw[,c('x','y')]), crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
crd <- sf::st_as_sf(x = cnf_summ_raw, coords = c('x','y'), crs = sf::st_crs(4326))
dst <- terra::distance(x = crd_terra)

# 2. Perform hierarchical clustering with Ward metric to test a range of k-values
# Hierarchical clustering
hcl <- hclust(d = dst, method = 'ward.D2')

  k <- 16
  grp <- stats::cutree(tree = hcl, k = k)
  # Group together the different convex hulls
  cvxL <- 1:k %>%
    purrr::map(.f = function(i){
      cvx <- sf::st_convex_hull(sf::st_union(crd[which(grp == i),]))
      cvx <- as(cvx, 'Spatial')
      return(cvx)
    })
  cvxL <- cvxL %>% raster::bind()
   
  
  if(class(cvxL) == 'list'){cvxL <- cvxL[[1]]}
  # Buffer the convex hulls with an extra range of 20 km
  cvxL <- raster::buffer(cvxL, width = 0.5)
  
  raster::shapefile(cvxL, paste0(root, "/data/", iso, "/conflict/conflict_area.tif"))
  
  x11();plot(cvxL); points(cnf_summ_raw[,c("x", "y")])
  
  
  w <- spatstat.geom::owin(xrange=c(raster::extent(cvxL)@xmin,
                               raster::extent(cvxL)@xmax),
                      yrange =c(raster::extent(cvxL)@ymin,
                                raster::extent(cvxL)@ymax))
  
  occurrences_ppp <-ppp(cnf_summ_raw$x, cnf_summ_raw$y, window=w) 
  bw_dig <- bw.diggle(occurrences_ppp)
  kernel <- density.ppp(x=occurrences_ppp,sigma=bw_dig,at="pixels",verbose=F,diggle=T)
  kernel <- raster::raster(kernel) %>% 
    raster::mask(., shp)
  kernel[] <- kernel[]/max(kernel[] , na.rm = T)
  writeRaster(kernel, paste0(root, "/data/", iso, "/conflict/conflict_kernel_densit.tif" )) 
  
  # Intersect the buffered convex hulls with the country shapefile
  # roiC <- raster::intersect(cvxL, rgeos::gUnaryUnion(shp))
  # roiC <- rgeos::gUnaryUnion(roiC)
  # roiC@bbox <- shp@bbox
  # 
  roiCtr <- ref #terra::rasterize(x = terra::vect(roiC), y = ref)
  
   
  cnf_summ <- cnf_summ_raw %>% 
    dplyr::mutate(cellID = raster::cellFromXY(roiCtr, xy = as.matrix(.[, c("x", "y")]))) %>% 
    dplyr::select(cellID, x, y , everything(.))
  
  
  
  bf_sp <- raster::buffer(SpatialPoints(coords = cnf_summ[, c("x","y")], 
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")), 
                       width = 30000)
 
  
  
  pseudo <- raster::mask(x = raster::raster(roiCtr), mask = bf_sp, inverse = T ) %>% 
    raster::as.data.frame(., xy = T) %>% 
    dplyr::slice_sample(., prop = 0.05) %>% 
    dplyr::mutate(cellID = raster::cellFromXY(roiCtr, xy = as.matrix(.[, c("x", "y")]))) %>% 
    dplyr::select(cellID, x, y ) 
   
  cnf_summ <-  bind_rows(cnf_summ, pseudo) %>%
    dplyr::mutate(across(where(is.numeric), function(i){ifelse(is.na(i), 0, i)})) 
  
 var_names <- names(cnf_summ %>% dplyr::select(-cellID, -x, -y))
 mtrs <- data.frame()  
 for(i in var_names){
   
   res1 <- rfsi.interpolation(root = root, 
                              iso = iso , 
                              df = cnf_summ, 
                              vr = i, 
                              ref = ref, 
                              crop_layer = NULL)
   
   writeRaster(res1$raster, paste0(root, "/data/", iso, "/conflict/", i, "_rfsi.tif"), overwrite =T)
   #riteRaster(res1$raster, "D:/OneDrive - CGIAR/Attachments/Desktop/mapas/rsfi5p_psuedo_knl.tif", overwrite =T)
   
   mtrs <- bind_rows(mtrs, res1$metrics)
 }
  
 write.csv(mtrs, paste0(root, "/data/", iso, "/conflict/rfsi_metrics.csv"))
  
  
  # lapply(list.files("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/ZWE/conflict/old/",
  #                   full.names = T, pattern = ".tif$"), function(i){
  #                     print(is(roiC))
  #                     r <- raster(i)  %>% 
  #                       raster::mask(., roiC)
  #                     
  #                     writeRaster(r, paste0("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/ZWE/conflict/", names(r), ".tif"))
  #                   })
  # 
  # 
  ## Add pseudo non-conflict reported coordinates
 