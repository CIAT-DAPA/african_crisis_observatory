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

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'

iso <- 'KEN'
country <- 'Kenya'

# Load and identify impact pathways
summ <- read.csv(file = paste0(root,'/data/',iso,'/_results/hotspots/soc_eco_selected_variables.csv'))
summ$Code <- gsub(pattern = '{iso}', replacement = iso, x = summ$Code, fixed = T)
ip_id <- unique(summ$IP_id)

# Global mask 1 km resolution
msk  <- terra::rast(paste0(root,'/data/_global/masks/mask_world_1km.tif'))

pth  <- paste0(root, '/data/',iso)

# Country shapefile
shp  <- terra::vect(paste0(pth,'/_shps/',iso,'.shp'))

# Raster template
tmp  <- msk %>% 
  terra::crop(x = ., y = terra::ext(shp)) %>% 
  terra::mask(mask = shp)

ip_id %>%
  purrr::map(.f = function(ip){
    tb <- summ %>% dplyr::filter(IP_id == ip)
    
    regions <- stringr::str_split(unique(tb$Region_value), ";") %>% 
      unlist %>%
      stringr::str_trim()
    
    var_name <- unique(tb$Region_key)
    
    stmp <- raster::shapefile(paste0(pth,'/_shps/',iso,'.shp'))
    shp_region <- stmp[stmp@data %>% pull(!!var_name) %in% regions, ] 
    shp_region <- as(shp_region, 'SpatVector'); rm(stmp)
    
    tmp  <- msk %>% terra::crop(x = ., y = terra::ext(shp_region)) 
    
    # Load raster variables
    htp <- tb$Code %>%
      purrr::map(.f = function(var){
        
        # Raster template
        r <- terra::rast(list.files(path = pth, pattern = paste0(var,'.tif$'), full.names = T, recursive = T))
        r <- r %>%
          terra::crop(x = ., y = terra::ext(tmp)) %>%
          terra::resample(x = ., y = tmp) %>%
          terra::mask(mask = shp_region)
        
        thr <- tb$Threshold[tb$Code == var]
        prc <- tb$Percentile[tb$Code == var]
        
        if(!is.na(thr)){
          eval(parse(text = paste0('r[!(r ',thr,')] <- NA')))
        }
        
        qtl <- global(x = r, fun = quantile, probs = prc, na.rm = T) %>% as.numeric()
        
        if(prc > 0.5){
          r[r < qtl] <- NA
          r[r >= qtl] <- 1
        } else {
          r[r > qtl] <- NA
          r[r <= qtl] <- 1
        }
        
        return(r)
        
      })
    rst <- htp %>% terra::rast() %>% sum(na.rm = T)
    rst[rst == 0] <- NA
    out <- paste0('D:/',iso,'_all_hotspots_',ip,'.tif')
    # out <- paste0(root,'/data/',iso,'/_results/hotspots/',iso,'_all_hotspots_',ip,'.tif')
    dir.create(path = dirname(out), showWarnings = F, recursive = T)
    terra::writeRaster(rst, filename = out, overwrite = T)
  })
