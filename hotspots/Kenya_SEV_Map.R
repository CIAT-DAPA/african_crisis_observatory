# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain SEV maps for Kenya
# Author: Andres Mendez
# Modified by : Victor Korir
# Alliance Bioversity International - CIAT, 2024
# ----------------------------------------------------------------------------------- #

# R options
.rs.restartR()                      # Restart R session
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,geojsonsf, readxl, geojsonlint, RColorBrewer, writexl, raster,terra, sf, stringr, stringi, lattice, rasterVis, maptools,
                                latticeExtra, RColorBrewer,cowplot, grid,tmap, tmaptools, geojson, geojsonio, MetBrewer, paletteer, exactextractr))

#' Variable definition
#'

iso <- "KEN"
baseDir <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/"
root <- paste0(baseDir, iso, "/")
scale_bar_pos <- switch( iso, "ZWE" = "left", "KEN" = "left", "UGA" = "right", "MLI" = "left", "SEN" = "left", "NGA" = "right", "SDN" = "right", 'PHL'="right", 'GTM'="right", "NER" = "right", "BFA" = "right", "SOM" = "right", "RWA" = "right", "MOZ"= "right")
scale_bar_top <- switch( iso, "ZWE" = "bottom", "KEN" = "bottom", "UGA" = "bottom", "MLI" = "bottom", "SEN" = "top", "NGA" = "bottom", "SDN" = "bottom", 'PHL'="top", 'GTM'="bottom", "NER" = "bottom", "BFA" = "bottom","SOM" = "bottom","RWA" = "bottom", "MOZ" = "bottom")


create_labels <- function(text, type = c("short", "long")){
  
  if(type == "short"){
    
    to_ret <- sapply(text, function(i){
      
      if(grepl(" ", i)){
        
        ret1 <- stringr::str_split(str_to_upper(i), " ") %>% 
          unlist %>%  
          str_extract(., "[A-Z]") %>% 
          paste0( ., collapse = '') 
        
        
      }else{
        ret1 <- i %>%   
          str_extract(., "[A-Z]")
        
      }
      return(ret1)
    }, simplify = "array" )
    
  }else{
    
    to_ret <- sapply(text, function(i){
      if(grepl(" ", i)){
        ret1 <- stringr::str_split(str_to_upper(i), " ") %>% 
          unlist %>%  
          str_extract(., "[A-Z]") %>% 
          paste0( ., collapse = '') 
        
        ret1 <- ifelse(!is.na(i), paste0(i, "(", ret1, ")"), NA)
        
      }else{
        
        ret1 <- i %>%   
          str_extract(., "[A-Z]")
        
        ret1 <- ifelse(!is.na(i), paste0(i, "(",  i %>%   
                                           str_extract(., "[A-Z]"), ")" ), NA )  
      }
      return(ret1)
    })
    
  }
  
  return(to_ret)
}

fix_label <- function( rast_labs, labs = av_labs){
  
  to_compare <- rast_labs %>% 
    dplyr::filter(!grepl("\\+", final_short_lab)) %>% 
    dplyr::pull(final_short_lab) %>% 
    str_split(., "\\+", simplify = T) %>% 
    c %>% 
    .[nchar(.) != 0] %>% 
    unique
  
  vec <- sapply(rast_labs$final_short_lab, function(txt){
    
    spt <- str_split(txt, "\\+", simplify = T) %>% c
    
    if(!all(spt %in% to_compare)){
      
      lab_diff <- setdiff(spt, to_compare)
      to_replace <- labs %>% 
        dplyr::filter(final_short_lab == lab_diff) %>% 
        dplyr::pull(long_labs)
      ret <- str_replace(txt, patter = lab_diff, replacement = to_replace)
      
    }else if(!grepl("\\+", txt)){
      ret <- labs %>% 
        dplyr::filter(final_short_lab == txt) %>% 
        dplyr::pull(long_labs)
    }else{
      ret <- txt
    }
    
    return(ret)
  }, simplify = T)
  
  return(vec)
  
}

cats_extract <- function(values, coverage_fractions){
  
  
  df_in <- data.frame(values, coverage_fractions) %>% 
    drop_na()
  
  ret <- tibble(ID = unique(df_in$values)) %>% 
    dplyr::left_join(., rast_labs %>% dplyr::select(ID, final_label ), by = c("ID" = "ID"))
  
  
  return(list(ret)) 
}

ip_text_description <- function(shp_object ,ip, df, n_vars = 10){
  
  
  ip_var_list <- read_csv(paste0(root,"_results/hotspots/soc_eco_all_variables.csv")) %>% 
    dplyr::mutate(Code = ifelse(grepl("\\{iso\\}", Code), gsub("\\{iso\\}", iso, Code), Code),
                  Code = tolower(Code))
  
  
  fl <- list.files(paste0(root,"_results/hotspots/"), pattern = paste0("_sorted_",unique(shp$NAME_0), '_', ip, ".xlsx"), full.names = T)
  tb <- readxl::read_excel(fl) %>% 
    dplyr::mutate(Code = tolower(Code)) %>% 
    dplyr::left_join(., ip_var_list %>% 
                       dplyr::filter(IP_id == ip ) %>% 
                       dplyr::select(-Variable, - Classification), by = c("Code" = "Code")) %>% 
    dplyr::mutate(Code = ifelse(grepl("_awe", Code), paste0(iso,"_AWE"), Code),
                  Code = ifelse(grepl("_rwi", Code), paste0(iso, "_rwi"), Code)) %>% 
    dplyr::slice(1:n_vars) %>% 
    dplyr::mutate(Variable = stringr::str_replace(Variable, pattern = "[0-9]+:", replacement = ""),
                  Variable = stringr::str_trim(Variable)) %>% 
    dplyr::select(Code, Variable, Percentile)
  
  df_proc <- df %>% 
    dplyr::filter(IP_id == ip) %>%
    dplyr::left_join(., tb , by = c("Code" = "Code")) %>% 
    drop_na(Variable) %>%
    dplyr::mutate(prefix = ifelse(Percentile == 0.9, paste0("[Values greather than ", Percentile, " distribution percentile]"),
                                  paste0("[Values lower than ", Percentile, " distribution percentile]") ),
                  txt_descp = paste(Variable, prefix, paste0("[",source_file ,"]") )) 
  
  text_raw <- lapply(1:nrow(df_proc), function(k){
    r <- raster(df_proc$file_path[k])
    lst <- exactextractr::exact_extract(r, sf::st_as_sf(shp_object))
    
    txt_out <- sapply(lst, function(i){
      if(length(na.omit(i$value))> 0){
        
        ret <- df_proc$txt_descp[k]
      }else{
        ret <- ""
      }
      
      return(ret)
    })
    
    return(txt_out)
    
    
    
  }) %>% 
    do.call(cbind, .) %>% 
    apply(., 1, function(l){
      
      if(all(l == "")){
        ret <- ""
      }else{
        ret <- paste(l[nchar(l) > 0], collapse = ";")
      }
      
    })
  
  return(text_raw)
}#end function

w_mask <- raster::raster(paste0(baseDir, "_global/masks/mask_world_1km.tif"))



to_share_dir <- paste0(root, "_results/", iso, "_to_share")
if(!dir.exists(to_share_dir)){dir.create(to_share_dir)}

##################################################################
######### Generate climate and conflict intersection #############
##################################################################

shp_c  <- st_read(paste0(root, "_shps/", iso, ".shp"))
shp <- shp_c

c_mask <- w_mask %>% 
  terra::crop(., extent(shp_c)) %>% 
  terra::mask(., shp_c)

conf_clust <- st_read(paste0(root, "_results/cluster_results/conflict/conflict_regular_clust.shp"))
clim_clust <- st_read(paste0(root, "_results/cluster_results/climate/climate_regular_clust.shp"))

st_crs(conf_clust) <- st_crs(c_mask)
st_crs(clim_clust) <- st_crs(c_mask)



clim_clust_labs <- readr::read_csv(paste0(baseDir, iso, "/_results/cluster_results/climate/climate_reg_cluster_text_description.csv")) %>% 
  dplyr::mutate(across(everything(.), as.character))

conf_clust_labs <- read_csv(paste0(root, "_results/cluster_results/conflict/conflict_cluster_text_description.csv")) %>%
  dplyr::mutate(across(everything(.), as.character))


conf_clust <-  conf_clust %>% 
  left_join(., conf_clust_labs , by = c("clust")) %>% 
  dplyr::rename(clust_km = short_label)  %>% 
  rename_with(., function(i){return("conflict_cluster_text_description")}, starts_with("clutert")) %>% 
  #dplyr::rename("label" = "clust") %>% 
  dplyr::mutate(short_label = stringr::str_extract(string = label, pattern = "[A-Za-z]+"),
                label = factor(label, levels = c("High conflict", "Moderate conflict", "Limited conflict"))) %>% 
  as_tibble()



clim_clust <- clim_clust %>%
  dplyr::mutate(clust = as.character(clust)) %>%
  dplyr::left_join(., clim_clust_labs , by = c("clust" = "clust")) %>% 
  as_tibble()

conf_clust <- st_as_sf(conf_clust)
clim_clust <- st_as_sf(clim_clust)
#conf_clust@data$clim_cluster <- as.character(sp::over(conf_clust, clim_clust, returnList = F)$text_output)
conf_clim <- st_intersection(conf_clust, clim_clust)
conf_clim <- conf_clim %>% dplyr::rename(clim_cluster = text_output)
conf_clim <- conf_clim %>% dplyr::rename(clim_cluster_short_label = label)

#conf_clust$clim_cluster <- as.character(unlist(st_intersects(st_as_sf(conf_clust), st_as_sf(clim_clust))))

#sf_ov <- as.character(unlist(sf_ov))
#conf_clust@data$clim_cluster <- sf_ov


#conf_clust@data$clim_cluster_short_label <- as.character(sp::over(conf_clust, clim_clust, returnList = F)$label)
#conf_clust@data$clim_cluster_order <- as.numeric(sp::over(conf_clust, clim_clust, returnList = F)$order)

conf_clim$intersect_conf_clim <- paste0(conf_clim$label, "-[", conf_clim$clim_cluster_short_label,"]")


conf_data <- read_csv(paste0(root, "conflict/", iso, "_conflict.csv"))

conf_occ <- conf_data %>% 
  dplyr::select(LONGITUDE, LATITUDE, FATALITIES) %>% 
  dplyr::filter(FATALITIES > 0) %>% 
  dplyr::group_by(LONGITUDE, LATITUDE) %>% 
  dplyr::summarise(FATALITIES = sum(FATALITIES)) %>% 
  dplyr::ungroup()

coordinates(conf_occ) <- ~LONGITUDE+LATITUDE
crs(conf_occ) <- crs(w_mask)
conf_occ <- st_as_sf(conf_occ)
inters <-st_intersects(st_as_sf(conf_occ), shp_c)
for(i in 1:nrow(conf_occ)){
  conf_occ[, 'over'] <- lapply(inters[i], function(x){shp_c$NAME_0[x]}) %>% unlist() %>% unique() %>% paste(., collapse = ",")
}
#conf_occ$over <- st_intersection(st_as_sf(conf_occ) ,shp_c )$NAME_0
conf_occ <- conf_occ[!is.na(conf_occ$over),]


get_ip_names <- "ip_all"


ip_x <- get_ip_names
  
cat(">>> creating maps for: ", ip_x, "\n")
  
ip_codes <- read_csv(paste0(root, "_results/hotspots/ip_maps/", iso, "_hotspots_values.csv")) %>% 
    dplyr::filter(ip == ip_x) 
  
  
labs_tbl <- lapply(1:length(unique(ip_codes$category)), function(i){
    
    x <- utils::combn(ip_codes$category, m = i ) %>% 
      t(.) %>% 
      tibble::as_tibble(.)
    
    return(x)
  }) %>% 
    dplyr::bind_rows()
  
  
vals_tbl <- apply(labs_tbl, 2, function(i){
    data.frame(V1 = i) %>% 
      left_join(., ip_codes %>% dplyr::select(category, raster_value), by = c("V1" = "category")) %>% 
      dplyr::pull(raster_value)
    
  }) %>% 
    as_tibble() %>% 
    dplyr::mutate(rast_values = rowSums(., na.rm = T))
  
short_labels_tbl <- apply(labs_tbl, 2, create_labels, type = "short") %>% 
    as_tibble() %>% 
    dplyr::mutate(final_short_lab = apply(., 1, function(i){
      
      ret <- paste0(i[!is.na(i)], collapse = "+")
      
      return(ret)
    }))
  
long_labels_tbl <- apply(labs_tbl , 2, create_labels, type = "long") %>% 
    as_tibble() %>% 
    dplyr::select(long_labs = V1) %>% 
    dplyr::bind_cols(., vals_tbl %>% dplyr::select(rast_values)) %>% 
    dplyr::slice(1:length(unique(ip_codes$category)))
  
final_label_tbl <- short_labels_tbl %>% 
    dplyr::select(final_short_lab) %>% 
    bind_cols(., vals_tbl %>% dplyr::select(rast_values)) %>% 
    left_join(., long_labels_tbl) %>%
    left_join(., ip_codes %>% dplyr::select(category, raster_value), by = c("rast_values" = "raster_value"))
  
tb <- read.csv(file = paste0(baseDir,iso,'/_results/hotspots/soc_eco_all_variables.csv')) # soc_eco_all_variables.csv # soc_eco_selected_variables.csv
tb$Code <- gsub(pattern = '{iso}', replacement = iso, x = tb$Code, fixed = T)
tb <- tb %>% dplyr::filter(IP_id == ip_x)
rg <- stringr::str_split(unique(tb$Region_value), ";") %>% unlist() %>% stringr::str_trim()
var_name <- unique(tb$Region_key)
shp_c <- shp[shp %>% dplyr::pull(!!var_name) %in% rg,]
  
  
ht_rast <- terra::rast(paste0(root,"_results/hotspots/ip_maps/",iso, '_all_hotspots_', ip_x, ".tif" )) %>% 
    terra::crop(., extent(shp)) %>% 
    terra::mask(., shp)
  
terra::writeRaster(ht_rast, paste0(to_share_dir, "/", ip_x, "_hotspots_map.tif"), overwrite = T)
  
rast_labs <- tibble( ID = unique(ht_rast[]) ) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(., final_label_tbl, by = c("ID" = "rast_values")) %>%
    dplyr::mutate(final_label = fix_label(rast_labs = ., labs =  final_label_tbl %>%
                                            dplyr::filter(!is.na(category)) ),chars = nchar(final_label)) %>% 
    dplyr::arrange(desc(chars))%>%
    dplyr::mutate(seq = 1:nrow(.))
  
  
rast_labs$ID <-as.numeric(rast_labs$ID)
rast_labs <-rast_labs %>%
    dplyr::select(value_ID = ID, label = final_label) %>% 
    write_csv(., paste0(to_share_dir, "/", ip_x, "_hotspots_labels.csv"))
  
  

  
names(rast_labs) <- c('ID', 'final_label')
ext<- exactextractr::exact_extract(ht_rast, sf::st_as_sf(conf_clim), fun = cats_extract)
  
conf_clim <- conf_clim %>% 
    dplyr::mutate(!!paste0(ip_x, "_rast_values") :=  lapply(ext, function(vec){paste(vec$ID, collapse = ";") }) %>% do.call( rbind, .),
                  !!paste0(ip_x, "_category") :=  lapply(ext, function(vec){paste(vec$final_label, collapse = ";") }) %>% do.call( rbind, .) ) 
  
    
cat_ids <- rast_labs %>% 
      dplyr::mutate(seq = nchar(final_label)) %>% 
      dplyr::arrange(desc(seq)) %>% 
      dplyr::mutate(seq = 1:nrow(.))
    
n_colors <- length(cat_ids$final_label)
    
    # Convert cat_ids to a matrix, containing only the relevant columns
cat_ids_matrix <- as.matrix(cat_ids[, c("ID", "seq")])
manual_labels <- c('Resource Scarcity(RS)', 'Undernutrition(U)','Inequality(I)',
                       'I+RS+U','I+RS', 'RS+U','I+U')
    # Use classify instead of subs
ht_rast_f <- terra::classify(ht_rast, cat_ids_matrix)
    
unique_values <- unique(values(ht_rast_f, na.rm = TRUE))
    
raster_levels <- data.frame(ID = unique_values, 
                                final_label = manual_labels)
    
    # Assign these levels to the raster
levels(ht_rast_f) <- list(raster_levels)
    

all_ip_map <- tmap::tm_shape(shp)+
    tm_borders(col = "gray50")+
    tm_shape(ht_rast_f)+
    tm_raster( style = "cat", palette = get_brewer_pal('Pastel1', n = n_colors, contrast = NA, stretch = TRUE, plot = TRUE), title =  paste("Hotspots", toupper(ip_x)))+
    tm_compass(type = "8star", position = c("right", "top")) +
    tm_scale_bar(breaks = c(0, 100, 200), text.size = 1, position = c(scale_bar_pos, scale_bar_top))+
    tm_layout(legend.outside=T, 
                legend.text.size = 1.3,
                legend.title.size= 1.3,
                legend.frame=F, 
                #legend.position=c(0.985, 0.985),
                #legend.just = c("left", "top"), 
                legend.outside.position = "right",
                legend.outside.size = 0.45,
                #legend.width= 1,
                outer.margins = c(0.01,0,0.01,0),
                inner.margins = c(0.01,0.01,0.01,0.01),
                legend.height= -0.2 )
    
tmap_save(all_ip_map,
              filename= paste0(root, "_results/hotspots/ip_maps/full_ip_categories_map_SEV.png"),
              dpi=300, 
              #insets_tm= insetmap,
              #insets_vp=viewport(x= 0.44, y= 0.78, width= 0.2, height= 0.2),
              height=9,
              width=16,
              units="in")
    

  
 