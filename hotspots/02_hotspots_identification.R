# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain all socio-economic hotspots filtered by mask area and climate-conflict clusters intersection
# Steps:
# 1. Read impact pathways table to identify the most vulnerable socio-economic
#    conditions per country/region
# 2 Sum them all applying the corresponding threshold (10% or 90% of the distribution)
# 2. Execute this script to obtain:
#    Raster file with hotspots
# Author: Harold Achicanoy, Andres Mendez
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'C:/Users/acmendez/OneDrive - CGIAR/African_Crisis_Observatory'

iso <- 'KEN'
country <- 'Kenya'

# Load and identify impact pathways
summ <- readxl::read_excel(path = paste0(root,'/Country_pathways.xlsx'), sheet = 2) %>%
  dplyr::filter(Country == country & Dimension != 'Climate')
ip_id <- unique(summ$IP_id)
vars <- as.character(summ$Variable)
summ <- summ[which(!is.na(vars)),]; rm(vars)

# Global mask 1 km resolution
msk  <- raster::raster(paste0(root,'/data/_global/masks/mask_world_1km.tif'))

pth  <- paste0(root, '/data/',iso)

# Country shapefile
shp  <- raster::shapefile(paste0(pth,'/_shps/',iso,'.shp'))

# Raster template
tmp  <- msk %>% 
  raster::crop(x = ., y = raster::extent(shp)) %>% 
  raster::mask(mask = shp)

ip_id %>%
  purrr::map(.f = function(ip){
    tb <- summ %>% dplyr::filter(IP_id == ip)
    
    
    regions <- stringr::str_split(unique(tb$Region_value), ";") %>% 
      unlist %>% 
      stringr::str_trim()
    
    var_name <- unique(tb$Region_key)
    
    shp_region <- shp[shp@data %>% pull(!!var_name) %in% regions, ] 
    
    tmp  <- msk %>% 
      raster::crop(x = ., y = raster::extent(shp_region)) 
    
    
    # Load raster variables
    htp <- tb$Variable %>%
      purrr::map(.f = function(var){
        
        # Raster template
        r <- raster::raster(list.files(path = pth, pattern = paste0(var,'.tif$'), full.names = T, recursive = T))
        r <- r %>% 
          raster::crop(x = ., y = raster::extent(tmp)) %>% 
          raster::resample(x = ., y = tmp)%>% 
          raster::mask(mask = shp_region)
        
        thr <- tb$Threshold[tb$Variable == var]
        prc <- tb$Percentile[tb$Variable == var]
        
        if(!is.na(thr)){
          eval(parse(text = paste0('r[!(r[] ',thr,')] <- NA')))
        }
        
        qtl <- raster::quantile(r, probs = prc)
        
        if(prc > 0.5){
          r[r[] < qtl] <- NA
          r[r[] >= qtl] <- 1
        } else {
          r[r[] > qtl] <- NA
          r[r[] <= qtl] <- 1
        }
        
        return(r)
        
      })
    rst <- htp %>% raster::stack() %>% sum(na.rm = T)
    rst[rst[] == 0] <- NA
    out <- paste0(root,'/data/',iso,'/_results/hotspots/',iso,'_all_hotspots_',ip,'.tif')
    dir.create(path = dirname(out), showWarnings = F, recursive = T)
    raster::writeRaster(rst, filename = out, overwrite = T)
  })

# SDN: check soil_carbon
