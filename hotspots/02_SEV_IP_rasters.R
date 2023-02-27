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

<<<<<<< HEAD:hotspots/02_hotspots_identification.R
iso <- 'SDN'
=======
iso <- 'GTM'
>>>>>>> cd51984c600e751d181699138d24c4782033365d:hotspots/02_SEV_IP_rasters.R
country <-  switch (iso,
                    "KEN" = "Kenya",
                    "SEN" = "Senegal",
                    "NGA" = "Nigeria",
                    "UGA" = "Uganda",
                    "ZWE" = "Zimbabwe",
                    "MLI" = "Mali",
                    "SDN"  = "Sudan",
                    "GTM"  = "Guatemala",
                    "PHL" = "Philippines"
)
#' This scrips makes SEV maps
# Load and identify impact pathways
ip_var_list <- read_csv(paste0(root, "/data/", iso, "/_results/hotspots/soc_eco_all_variables.csv")) %>% 
  dplyr::mutate(Code = ifelse(grepl("\\{iso\\}", Code), gsub("\\{iso\\}", iso, Code), Code),
                Code = tolower(Code))




#summ <- read.csv(file = paste0(root,'/data/',iso,'/_results/hotspots/soc_eco_all_variables.csv')) # soc_eco_all_variables.csv # soc_eco_selected_variables.csv
#summ$Code <- gsub(pattern = '{iso}', replacement = iso, x = summ$Code, fixed = T)

ip_id <- unique(ip_var_list$IP_id)

# Global mask 1 km resolution
msk <- terra::rast(paste0(root,'/data/_global/masks/mask_world_1km.tif'))
pth <- paste0(root, '/data/',iso)
out_pth <- paste0(pth, "/_results/hotspots/ip_maps/")
if(!dir.exists(out_pth)){dir.create(out_pth)}

# Country shapefile
shp <- terra::vect(paste0(pth,'/_shps/',iso,'.shp'))
stmp <- raster::shapefile(paste0(pth,'/_shps/',iso,'.shp'))

# Raster template
tmp <- msk %>% 
  terra::crop(x = ., y = terra::ext(shp)) %>% 
  terra::mask(mask = shp)

n_vars <- 10

ip_id %>%
  purrr::map(.f = function(ip){
    
    fl <- list.files(paste0(root, "/data/",iso,  "/_results/hotspots/"), pattern = paste0("_sorted_",unique(shp$NAME_0), '_', ip, ".xlsx"), full.names = T)
    
    
    tb <- readxl::read_excel(fl) %>% 
      dplyr::mutate(Code = tolower(Code)) %>% 
      dplyr::left_join(., ip_var_list %>% 
                         dplyr::filter(IP_id == ip ) %>% 
                         dplyr::select(-Variable, - Classification), by = c("Code" = "Code")) %>% 
      dplyr::mutate(Code = ifelse(grepl("_awe", Code), paste0(iso,"_AWE"), Code),
                    Code = ifelse(grepl("_rwi", Code), paste0(iso, "_rwi"), Code)) %>% 
      dplyr::slice(1:n_vars)
    
    
    regions <- stringr::str_split(unique(tb$Region_value), ";") %>% 
      unlist %>%
      stringr::str_trim()
    
    var_name <- unique(tb$Region_key)
    
    
    shp_region <- stmp[stmp@data %>% pull(!!var_name) %in% regions, ] 
    shp_region <- as(shp_region, 'SpatVector')
    
    tmp  <- msk %>% 
      terra::crop(x = ., y = terra::ext(shp_region)) %>% 
      terra::mask(mask = shp_region)
    
    # Load raster variables
    htp <- 1:length(tb$Code) %>%
      purrr::map(.f = function(i){
        
        var <- tb$Code[i]
        
        cat(">>> Processing: ", var, "\n")
        # Raster template
        r <- terra::rast(list.files(path = pth, pattern = paste0(var,'.tif$'), full.names = T, recursive = T))
        
        thr <- tb$Threshold[i] 
        prc <- tb$Percentile[i] 
        
        if(!is.na(thr)){
          eval(parse(text = paste0('r[!(r ',thr,')] <- NA')))
        }
         
        mn <- median(r[], na.rm = T)
        qtl_country <- global(x = r, fun = quantile, probs = prc, na.rm = T) %>% as.numeric()
        
        
        r <- r %>%
          terra::crop(x = ., y = terra::ext(tmp)) %>%
          terra::resample(x = ., y = tmp) %>%
          terra::mask(mask = shp_region)
        
        to_check <- (min(r[], na.rm = T) == quantile(r[], probs = 0.25, na.rm = T) & prc == 0.1) |
          (max(r[], na.rm = T) == quantile(r[], probs = 0.75, na.rm = T) & prc == 0.9)
          
        
        stopifnot("Raster values distribution very Skewed" = !to_check)
    
          if(prc > 0.5){
          
              r[r < qtl_country] <- NA
              r[r >= qtl_country] <- 1
            
          } else {
              r[r > qtl_country] <- NA
              r[r <= qtl_country] <- 1
           
          }
        
     
         return(r)
        
      })
    rst <- htp %>% terra::rast() %>% sum(na.rm = T)
    rst[rst == 0] <- NA
    out <- paste0(root,'/data/',iso,'/_results/hotspots/ip_maps/',iso,'_all_hotspots_',ip,'.tif') # out <- paste0('D:/',iso,'_all_hotspots_',ip,'.tif')
    dir.create(path = dirname(out), showWarnings = F, recursive = T)
    terra::writeRaster(rst, filename = out, overwrite = T)
  })


