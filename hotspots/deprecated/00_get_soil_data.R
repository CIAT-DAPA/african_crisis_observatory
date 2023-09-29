rm(list=ls(all=TRUE))
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/AWCPTF.R')

# Input parameters:
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




library(future)
library(furrr)
library(tidyverse)
future::plan(multisession, workers = 11)

isos <- c('SDN','ZWE','SEN','MLI','NGA','KEN','UGA', "ETH", "GTM", "PHL", "ZMB")

isos %>% 
  furrr::future_map(.x = ., .f = function(.x){
    iso <- .x
    cat("Processing: ", iso, "/n")
    shp_fl <- paste0('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/', iso,'/_shps/', iso, '.shp')
    
    out_dir_soil <- paste0('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/', iso,'/climatic_indexes/temp/')
    
    if(!dir.exists(out_dir_soil)){dir.create(out_dir_soil)}
    
    soil_cp  <- paste0('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/', iso,'/climatic_indexes/temp/soilcp.tif')
    soil_sat <- paste0('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/', iso,'/climatic_indexes/temp/soilsat.tif')
    
    if(file.exists(soil_cp) & file.exists(soil_sat)){
      print("file already exists")
    }else{
      get_soil(shp_fl = shp_fl, root_depth = 60, outfiles = c(soil_cp,
                                                              soil_sat))
      
    }
    
   return(NULL)  
  })

future::plan(sequential)
