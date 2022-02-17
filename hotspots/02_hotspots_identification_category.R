options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- "C:/Users/acmendez/OneDrive - CGIAR/African_Crisis_Observatory"

iso <- 'KEN'
country <- 'Kenya'

# Load and identify impact pathways
summ <- readxl::read_excel(path = paste0(root,'/Country_pathways.xlsx'), sheet = 2) %>%
  dplyr::filter(Country == country & Dimension != 'Climate')
ip_id <- unique(summ$IP_id)
vars <- as.character(summ$Variable)
summ <- summ[which(!is.na(vars)),]; rm(vars)

# Global mask 1 km resolution
msk  <- raster(paste0(root, "/data/", iso, "/mask/", iso, "_mask.tif")) 
  #raster::raster('D:/OneDrive - CGIAR/African_Crisis_Observatory/data/_global/masks/mask_world_1km.tif')

pth  <- paste0(root, "/data/",iso)

# Country shapefile
shp  <- raster::shapefile(paste0(pth,'/_shps/',iso,'.shp'))

# Raster template
tmp  <- msk # %>% raster::crop(x = ., y = raster::extent(shp)) %>% raster::mask(mask = shp)
tmp[!is.na(tmp)] <- 1
# Iterate by impact pathway
smm_df <- ip_id %>%
  purrr::map(.f = function(ip){
    # Filter full table by impact pathway
    tb <- summ %>% dplyr::filter(IP_id == ip)
    # Identify hotspot categories
    ct <- sort(unique(tb$Categories))
    # Iterate by category
    cts <- 1:length(ct) %>%
      purrr::map(.f = function(i){
        # Load raster variables within the category
        htp <- tb$Variable[tb$Categories == ct[i]] %>%
          purrr::map(.f = function(var){
            
            r <- raster::raster(list.files(path = pth, pattern = paste0(var,'.tif$'), full.names = T, recursive = T))
            r <- r %>% raster::crop(x = ., y = raster::extent(tmp)) %>% raster::resample(x = ., y = tmp) %>% 
              raster::mask(., mask = tmp)
            
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
        # Identify vulnerable areas per category
        rst <- htp %>% raster::stack() %>% sum(na.rm = T)
        rst[rst[] == 0] <- NA
        rst[!is.na(rst[])] <- 10^(i-1) # Values assignation with scientific notation
        out <- paste0(root,'/data/',iso,'/_results/hotspots/',iso,'_',ip,'_',tolower(abbreviate(ct[i])),'_hotspots.tif')
        dir.create(path = dirname(out), showWarnings = F, recursive = T)
        raster::writeRaster(rst, filename = out, overwrite = T)
        df <- data.frame(iso = iso, ip = ip, category = ct[i], raster_value = 10^(i-1))
        return(df)
      }) %>%
      dplyr::bind_rows()
    
    # Load category's hotspots
    htp_ct <- list.files(path = paste0(root,'/data/',iso,'/_results/hotspots'), pattern = paste0(iso,'_',ip,'_*[a-z]*_hotspots.tif$'), full.names = T)
    hotspots <- raster::stack(htp_ct)
    # Remove natural resources category from the identification of the hotspots
    cond <- grep(pattern = 'rsrs', x = names(hotspots))
    if(length(cond) > 0){hotspots <- hotspots[[base::setdiff(1:nlayers(hotspots),cond)]]}
    # Compute general hotspots
    hotspots <- hotspots %>% sum(na.rm = T)
    hotspots[hotspots[] == 0] <- NA
    out <- paste0(root,'/data/',iso,'/_results/hotspots/',iso,'_all_cat_hotspots_',ip,'.tif')
    raster::writeRaster(hotspots, filename = out, overwrite = T)
    
    return(cts)
  }) %>%
  dplyr::bind_rows()

write.csv(x = smm_df, file = paste0(root,'/data/',iso,'/_results/hotspots/',iso,'_hotspots_values.csv'), row.names = F)
