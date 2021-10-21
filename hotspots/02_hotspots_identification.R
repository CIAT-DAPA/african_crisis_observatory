options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

iso <- 'SDN'
country <- 'Sudan'

# Load and identify impact pathways
summ <- readxl::read_excel(path = paste0(root,'/Country_pathways.xlsx'), sheet = 2) %>%
  dplyr::filter(Country == country & Dimension != 'Climate')
ip_id <- unique(summ$IP_id)
vars <- as.character(summ$Variable)
summ <- summ[which(!is.na(vars)),]; rm(vars)

# Global mask 1 km resolution
msk  <- raster::raster('D:/OneDrive - CGIAR/African_Crisis_Observatory/data/_global/masks/mask_world_1km.tif')

pth  <- paste0('D:/OneDrive - CGIAR/African_Crisis_Observatory/data/',iso)

# Country shapefile
shp  <- raster::shapefile(paste0(pth,'/_shps/',iso,'.shp'))

# Raster template
tmp  <- msk %>% raster::crop(x = ., y = raster::extent(shp)) %>% raster::mask(mask = shp)

ip_id %>%
  purrr::map(.f = function(ip){
    tb <- summ %>% dplyr::filter(IP_id == ip)
    
    # Load raster variables
    htp <- tb$Variable %>%
      purrr::map(.f = function(var){
        
        r <- raster::raster(list.files(path = pth, pattern = paste0(var,'.tif$'), full.names = T, recursive = T))
        r <- r %>% raster::crop(x = ., y = raster::extent(tmp)) %>% raster::resample(x = ., y = tmp)
        
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
