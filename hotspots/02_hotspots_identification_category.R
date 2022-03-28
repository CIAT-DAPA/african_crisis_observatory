# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain socio-economic hotspots per category filtered by mask area and
# climate-conflict clusters intersection
# Steps:
# 1. Read impact pathways table to identify the most vulnerable socio-economic
#    conditions per country/region
# 2 Sum them all applying the corresponding threshold (10% or 90% of the distribution)
# 2. Execute this script to obtain:
#    Raster file with socio-economic hotspots summed by category
# Author: Harold Achicanoy, Andres Mendez
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO"

iso <- 'KEN'
country <- 'Kenya'

# Load and identify impact pathways
summ <- read.csv(file = paste0(root,'/data/',iso,'/_results/hotspots/soc_eco_all_variables.csv')) # soc_eco_all_variables.csv # soc_eco_selected_variables.csv
summ$Code <- gsub(pattern = '{iso}', replacement = iso, x = summ$Code, fixed = T)
ip_id <- unique(summ$IP_id)

# Global mask 1 km resolution
msk <- terra::rast(paste0(root, "/data/", iso, "/mask/", iso, "_mask.tif")) 
pth <- paste0(root, "/data/",iso)

# Country shapefile
shp <- terra::vect(paste0(pth,'/_shps/',iso,'.shp'))

# Raster template
tmp  <- msk # %>% raster::crop(x = ., y = raster::extent(shp)) %>% raster::mask(mask = shp)
tmp[!is.na(tmp)] <- 1

# Iterate by impact pathway
smm_df <- ip_id %>%
  purrr::map(.f = function(ip){
    # Filter full table by impact pathway
    tb <- summ %>% dplyr::filter(IP_id == ip)
    
    # Identify regions if it's the case
    regions <- stringr::str_split(unique(tb$Region_value), ";") %>% 
      unlist %>%
      stringr::str_trim()
    
    var_name <- unique(tb$Region_key)
    
    stmp <- raster::shapefile(paste0(pth,'/_shps/',iso,'.shp'))
    shp_region <- stmp[stmp@data %>% pull(!!var_name) %in% regions, ] 
    shp_region <- as(shp_region, 'SpatVector'); rm(stmp)
    
    tmp <- tmp %>% terra::crop(x = ., y = terra::ext(shp_region)) %>% terra::mask(mask = shp_region)
    
    # Identify hotspots categories
    ct <- sort(unique(tb$Classification))
    # Iterate by category
    cts <- 1:length(ct) %>%
      purrr::map(.f = function(i){
        # Load raster variables within the category
        htp <- tb$Code[tb$Classification == ct[i]] %>%
          purrr::map(.f = function(var){
            
            r <- terra::rast(list.files(path = pth, pattern = paste0(var,'.tif$'), full.names = T, recursive = T))
            r <- r %>% terra::crop(x = ., y = terra::ext(tmp)) %>% terra::resample(x = ., y = tmp) %>% 
              terra::mask(., mask = shp_region)
            
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
        # Identify vulnerable areas per category
        rst <- htp %>% terra::rast() %>% sum(na.rm = T)
        rst[rst == 0] <- NA
        rst[!is.na(rst)] <- 10^(i-1) # Values assignation with scientific notation
        out <- paste0(root,'/data/',iso,'/_results/hotspots/',iso,'_',ip,'_',tolower(abbreviate(ct[i])),'_hotspots.tif')
        dir.create(path = dirname(out), showWarnings = F, recursive = T)
        terra::writeRaster(rst, filename = out, overwrite = T)
        df <- data.frame(iso = iso, ip = ip, category = ct[i], raster_value = 10^(i-1))
        return(df)
      }) %>%
      dplyr::bind_rows()
    
    # Load category's hotspots
    htp_ct <- list.files(path = paste0(root,'/data/',iso,'/_results/hotspots'), pattern = paste0(iso,'_',ip,'_*[a-z]*_hotspots.tif$'), full.names = T)
    hotspots <- terra::rast(htp_ct)
    # Remove natural resources category from the identification of the hotspots
    # cond <- grep(pattern = 'rsrs', x = terra::sources(hotspots))
    # if(length(cond) > 0){hotspots <- hotspots[[base::setdiff(1:terra::nlyr(hotspots),cond)]]}
    # Compute general hotspots
    hotspots <- hotspots %>% sum(na.rm = T)
    hotspots[hotspots == 0] <- NA
    out <- paste0(root,'/data/',iso,'/_results/hotspots/',iso,'_all_cat_hotspots_',ip,'.tif')
    terra::writeRaster(hotspots, filename = out, overwrite = T)
    
    return(cts)
  }) %>%
  dplyr::bind_rows()

write.csv(x = smm_df, file = paste0(root,'/data/',iso,'/_results/hotspots/',iso,'_hotspots_values.csv'), row.names = F)
