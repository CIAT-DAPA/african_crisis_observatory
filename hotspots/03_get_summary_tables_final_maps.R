# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain maps, graphs and final tables
# Author: Andres Mendez
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
.rs.restartR()                      # Restart R session
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,geojsonsf, readxl, geojsonlint, RColorBrewer, writexl, raster,terra, sp, sf, stringr, stringi, lattice, rasterVis, maptools,
                                latticeExtra, RColorBrewer,cowplot, grid,tmap, tmaptools, geojson, geojsonio, MetBrewer, paletteer, exactextractr))

#' Variable definition
#'

iso <- "BFA"
baseDir <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/"
root <- paste0(baseDir, iso, "/")
scale_bar_pos <- switch( iso, "ZWE" = "left", "KEN" = "left", "UGA" = "right", "MLI" = "left", "SEN" = "left", "NGA" = "right", "SDN" = "right", 'PHL'="right", 'GTM'="right", "NER" = "right", "BFA" = "right")
scale_bar_top <- switch( iso, "ZWE" = "bottom", "KEN" = "bottom", "UGA" = "bottom", "MLI" = "bottom", "SEN" = "top", "NGA" = "bottom", "SDN" = "bottom", 'PHL'="top", 'GTM'="bottom", "NER" = "bottom", "BFA" = "bottom")


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

shp_c  <- raster::shapefile(paste0(root, "_shps/", iso, ".shp"))
shp <- shp_c

c_mask <- w_mask %>% 
  raster::crop(., extent(shp_c)) %>% 
  raster::mask(., shp_c)

conf_clust <- raster::shapefile(paste0(root, "_results/cluster_results/conflict/conflict_regular_clust.shp"))
clim_clust <- raster::shapefile(paste0(root, "_results/cluster_results/climate/climate_regular_clust.shp"))

crs(conf_clust) <- crs(c_mask)
crs(clim_clust) <- crs(c_mask)



clim_clust_labs <- readr::read_csv(paste0(baseDir, iso, "/_results/cluster_results/climate/climate_reg_cluster_text_description.csv")) %>% 
  dplyr::mutate(across(everything(.), as.character))

conf_clust_labs <- read_csv(paste0(root, "_results/cluster_results/conflict/conflict_cluster_text_description.csv")) %>%
  dplyr::mutate(across(everything(.), as.character))


conf_clust@data <-  conf_clust@data %>% 
  left_join(., conf_clust_labs , by = c("clust")) %>% 
  dplyr::rename(clust_km = short_label)  %>% 
  rename_with(., function(i){return("conflict_cluster_text_description")}, starts_with("cltrt")) %>% 
  #dplyr::rename("label" = "clust") %>% 
  dplyr::mutate(short_label = stringr::str_extract(string = label, pattern = "[A-Za-z]+"),
                label = factor(label, levels = c("High conflict", "Moderate conflict", "Limited conflict"))) %>% 
  as_tibble()



clim_clust@data <- clim_clust@data %>%
  dplyr::mutate(clust = as.character(clust)) %>%
  dplyr::left_join(., clim_clust_labs , by = c("clust" = "clust")) %>% 
  as_tibble()

conf_clust@data$clim_cluster <- as.character(sp::over(conf_clust, clim_clust, returnList = F)$text_output)
conf_clust@data$clim_cluster_short_label <- as.character(sp::over(conf_clust, clim_clust, returnList = F)$label)
conf_clust@data$clim_cluster_order <- as.numeric(sp::over(conf_clust, clim_clust, returnList = F)$order)

conf_clust@data$intersect_conf_clim <- paste0(conf_clust@data$label, "-[", conf_clust@data$clim_cluster_short_label,"]")


###################################################################
################## Generate conflict graphs #######################
###################################################################

conf_data <- read_csv(paste0(root, "conflict/", iso, "_conflict.csv"))

conf_occ <- conf_data %>% 
  dplyr::select(LONGITUDE, LATITUDE, FATALITIES) %>% 
  dplyr::filter(FATALITIES > 0) %>% 
  dplyr::group_by(LONGITUDE, LATITUDE) %>% 
  dplyr::summarise(FATALITIES = sum(FATALITIES)) %>% 
  dplyr::ungroup()

coordinates(conf_occ) <- ~LONGITUDE+LATITUDE
crs(conf_occ) <- crs(w_mask)
conf_occ@data$over <- conf_occ %over%  shp_c %>% dplyr::pull(NAME_0)
conf_occ <- conf_occ[!is.na(conf_occ$over),]
#x11()
mainmap <- tmap::tm_shape(shp_c)+
  tm_borders(col = "black")+
  tm_shape(conf_clust)+
  tm_fill(col = "label", palette = c("#d7191c", "#e5a03e", "#ffffbf"), alpha = 0.7, title = expression("Conflict clusters"))+
  #tm_borders(col ="black") + 
  tm_compass(type = "8star", position = c("left", "top")) +
  tm_scale_bar(breaks = c(0, 100, 200), text.size = 1, position = c(scale_bar_pos, scale_bar_top))+
  tm_shape(conf_occ)+
  tm_symbols(col = "#31b225", border.col = "black", size = "FATALITIES", scale = 3, title.size = "Number of fatalities")+
  tm_layout(legend.outside=T, 
            legend.text.size = 1.3,
            legend.title.size=1.3,
            legend.frame=F, 
            #legend.position=c(0.985, 0.985),
            legend.just = c("left", "top"), 
            #legend.width=-0.25,
            #outer.margins = c(0,0,0,0),
            #inner.margins = c(0,0,0,0)
            legend.height= -0.3 )

x11();mainmap
tmap_save(mainmap,
          filename= paste0(root, "_results/cluster_results/conflict/geographic_distr_conflict.png"),
          dpi=300, 
          #insets_tm=insetmap, 
          #insets_vp=vp,
          height=8,
          width=15,
          units="in")

##############################################
#### Barplot of number of EVENTS TYPE ########
##############################################

events_bp <- conf_data %>% 
  dplyr::select(LONGITUDE, LATITUDE, EVENT_TYPE, FATALITIES) %>% 
  dplyr::mutate(sp::over( sp::SpatialPoints(.[, c("LONGITUDE", "LATITUDE")], proj4string = crs(w_mask)), conf_clust, returnList = F)) %>% 
  dplyr::filter(!is.na(short_label)) %>% 
  dplyr::group_by(EVENT_TYPE, short_label) %>% 
  dplyr::tally() %>% 
  # dplyr::summarise(FATALITIES = sum(FATALITIES)) %>% 
  dplyr::ungroup() %>% 
  na.omit()

write_csv(events_bp %>%
            dplyr::group_by(EVENT_TYPE) %>% 
            dplyr::summarise(short_label, n, freq = prop.table(n)), paste0(root, "_results/cluster_results/conflict/EVENTS_TYPE_barplot_df.csv"))


g1 <- events_bp %>% 
  dplyr::mutate(short_label = factor(x = short_label, levels = c('High','Moderate','Limited'))) %>%
  ggplot( aes(x = short_label, y = n, fill = EVENT_TYPE))+
  geom_bar( stat = "identity")+
  scale_fill_brewer(palette = "Set3")+
  xlab("")+
  ylab("Number of events [counts]")+
  labs(fill = "Events type")+
  #theme(text = element_text(size = 20))+
  theme_bw(base_size = 20)

ggsave(g1,
       filename= paste0(root, "_results/cluster_results/conflict/EVENTS_TYPE_barplot.png"),
       dpi=300, 
       #insets_tm=insetmap, 
       #insets_vp=vp,
       height=6,
       width=8,
       units="in")

##############################################
#### Barplot of number of FATALITITES ########
##############################################

fata_bp <- conf_data %>% 
  dplyr::select(LONGITUDE, LATITUDE, EVENT_TYPE, FATALITIES) %>% 
  dplyr::mutate(sp::over( sp::SpatialPoints(.[, c("LONGITUDE", "LATITUDE")], proj4string = crs(w_mask)), conf_clust, returnList = F)) %>% 
  dplyr::filter(!is.na(short_label)) %>% 
  dplyr::group_by(EVENT_TYPE, short_label) %>% 
  dplyr::summarise(counts = sum(FATALITIES)) %>% 
  # dplyr::summarise(FATALITIES = sum(FATALITIES)) %>% 
  dplyr::ungroup() %>% 
  na.omit()

write_csv(fata_bp %>% 
            dplyr::group_by(EVENT_TYPE) %>% 
            dplyr::summarise(short_label, counts, freq = prop.table(counts))
          , paste0(root, "_results/cluster_results/conflict/FATALITIES_barplot_df.csv"))


g2 <- fata_bp %>% 
  dplyr::mutate(short_label  = factor(short_label , levels = c("High", "Moderate", "Limited"))) %>% 
  ggplot( aes(x = short_label, y = counts, fill = EVENT_TYPE))+
  geom_bar( stat = "identity")+
  scale_fill_brewer(palette = "Set3")+
  xlab("")+
  ylab("Number of fatalitites [counts]")+
  labs(fill = "Events type")+
  #theme(text = element_text(size = 20))+
  theme_bw(base_size = 20)

ggsave(g2,
       filename= paste0(root, "_results/cluster_results/conflict/FATALITIES_barplot.png"),
       dpi=300, 
       #insets_tm=insetmap, 
       #insets_vp=vp,
       height=6,
       width=8,
       units="in")



#######################################################################################
############## create clim cluster maps with labels ##################################
#####################################################################################


lvls_clim <- clim_clust@data %>% 
  dplyr::select(label, order) %>% 
  dplyr::filter(!duplicated(label)) %>% 
  dplyr::arrange(desc(order)) %>% 
  pull(label) 

clim_clust@data$label <- factor(clim_clust@data$label, levels = lvls_clim)

clim_map <- tmap::tm_shape(shp)+
  tm_borders(col = "black")+
  tm_shape(clim_clust)+
  tm_fill(col = "label", palette = "-YlOrRd" , alpha = 0.7, title = "Climatic clusters")+
  #tm_borders(col ="black") + 
  tm_compass(type = "8star", position = c("left", "top")) +
  tm_scale_bar(breaks = c(0, 100, 200), text.size = 1, position = c(scale_bar_pos, scale_bar_top))+
  tm_layout(legend.outside=T, 
            legend.text.size = 1,
            legend.title.size= 1.3,
            legend.frame=F,
            #legend.position=c(0.985, 0.985),
            legend.just = c("left", "top"), 
            legend.width= 1,
            #outer.margins = c(0,0,0,0),
            #inner.margins = c(0,0,0,0)
            legend.height= -0.2 )
#x11();clim_map

tmap_save(clim_map,
          filename= paste0(root, "_results/cluster_results/climate/climate_clusters_map.png" ),
          dpi=300, 
          #insets_tm=insetmap, 
          #insets_vp=vp,
          height=10,
          width=24,
          units="in")

######################################################################################
################## Create intersection of climate-conflict map #######################
######################################################################################
                                                      

#raster::shapefile(conf_clust, paste0(root, "_results/cluster_results/conflict_climate_intersection.shp"), overwrite = T)


#x11()

lapply(unique(conf_clust@data$label), function(i){
  
  cat("Making maps for: ", as.character(i) , "\n")
  
  tmp_df <- conf_clust
  tmp_df@data$label <-  as.character(tmp_df@data$label)
  tmp_df@data$intersect_conf_clim[tmp_df@data$label != i] <- "Others combinations"
  
  lvls <- tmp_df@data %>% 
    dplyr::select(clim_cluster_short_label, clim_cluster_order) %>% 
    dplyr::filter(!duplicated(clim_cluster_short_label)) %>% 
    dplyr::arrange(desc(clim_cluster_order)) %>% 
    pull(clim_cluster_short_label) %>% 
    paste0(i,"-[", .,"]") %>% 
    c(., "Others combinations" )
  
  tmp_df@data$intersect_conf_clim <- factor(tmp_df@data$intersect_conf_clim, levels = lvls)
  
  main_col <- switch(i, "High conflict" = "YlOrRd",
                     "Moderate conflict" =  "Oranges",
                     "Limited conflict" = "YlGn") ##DFDFDF
  
  other_color <- "#DADADA"

  if(i == "Limited conflict"){
    pal_g <- brewer.pal(length(unique(tmp_df@data$intersect_conf_clim))-1 , main_col)
  }else{
    pal_g <- rev(brewer.pal(length(unique(tmp_df@data$intersect_conf_clim))-1 , main_col))
  }
  pal_g[length(pal_g)+1] <- other_color
  
  
  inter_map <- tmap::tm_shape(shp)+
    tm_borders(col = "black")+
    tm_shape(tmp_df)+
    tm_fill(col = "intersect_conf_clim", palette = pal_g , alpha = 0.7, title = paste0(i, " cluster"))+
    #tm_borders(col ="black") + 
    tm_compass(type = "8star", position = c("right", "top")) +
    tm_scale_bar(breaks = c(0, 100, 200), text.size = 1, position = c(scale_bar_pos, scale_bar_top))+
    tm_layout(legend.outside=T, 
              legend.text.size = 1,
              legend.title.size= 1.3,
              legend.frame=F,
              #legend.position=c(0.985, 0.985),
              legend.just = c("left", "top"), 
              legend.width= 1,
              #outer.margins = c(0,0,0,0),
              #inner.margins = c(0,0,0,0)
              legend.height= -0.2 )
  #x11();inter_map
  
  tmap_save(inter_map,
            filename= paste0(root, "_results/cluster_results/",tolower(gsub(" ", "_", i)),"_climate_intersection.png"),
            dpi=300, 
            #insets_tm=insetmap, 
            #insets_vp=vp,
            height=10,
            width=24,
            units="in")
  
})


#################################################
###### Hotspots map for IP's ####################
#################################################

get_ip_names <- "ip_all"
  
  # list.files(paste0(root, "_results/hotspots/ip_maps")) %>% 
  # str_extract(.,"ip[0-9]") %>% 
  # unique() %>% 
  # na.omit()


for(i in get_ip_names){
  
  ip_x <- i
  
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
  shp_c <- shp[shp@data %>% dplyr::pull(!!var_name) %in% rg,]
  
  
  ht_rast <- raster::raster(paste0(root,"_results/hotspots/ip_maps/",iso, '_all_cat_hotspots_', ip_x, ".tif" )) %>% 
    raster::crop(., extent(shp_c)) %>% 
    raster::mask(., shp_c)
  
  raster::writeRaster(ht_rast, paste0(to_share_dir, "/", ip_x, "_hotspots_map.tif"), overwrite = T)
  
  rast_labs <- tibble( ID = unique(ht_rast[]) ) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(., final_label_tbl, by = c("ID" = "rast_values")) %>%
    dplyr::mutate(final_label = fix_label(rast_labs = ., labs =  final_label_tbl %>%
                                            dplyr::filter(!is.na(category)) ),
                  chars = nchar(final_label)) %>% 
    dplyr::arrange(desc(chars))%>%
    dplyr::mutate(seq = 1:nrow(.))
  
  
  rast_labs %>%
    dplyr::select(value_ID = ID, label = final_label) %>% 
    write_csv(., paste0(to_share_dir, "/", ip_x, "_hotspots_labels.csv"))
  

  ext<- exactextractr::exact_extract(ht_rast, sf::st_as_sf(conf_clust), fun = cats_extract)
  
  conf_clust@data <- conf_clust@data %>% 
    dplyr::mutate(!!paste0(ip_x, "_rast_values") :=  lapply(ext, function(vec){paste(vec$ID, collapse = ";") }) %>% do.call( rbind, .),
                  !!paste0(ip_x, "_category") :=  lapply(ext, function(vec){paste(vec$final_label, collapse = ";") }) %>% do.call( rbind, .) ) 
  
  #plotear mapa con todas las categorias
 if(TRUE){
   
   cat_ids <- rast_labs %>% 
     dplyr::mutate(seq = nchar(final_label)) %>% 
     dplyr::arrange(desc(seq)) %>% 
     dplyr::mutate(seq = 1:nrow(.))
   
   n_colors <- length(cat_ids$final_label)
   
   ht_rast_f <- raster::subs(ht_rast, cat_ids[, c("ID", "seq")])
   ht_rast_f <- as.factor(ht_rast_f)
  
   levels(ht_rast_f) <- data.frame(id = levels(ht_rast_f)[[1]], x = cat_ids$final_label, row.names = 1:length(cat_ids$final_label)   ) 
     
    
    
   #paletteer_d("wesanderson::Rushmore", n_colors)
   all_ip_map <- tmap::tm_shape(shp_c)+
     tm_borders(col = "gray50")+
     tm_shape(ht_rast_f)+
     tm_raster( style = "cat", palette = paletteer_dynamic("ggthemes_ptol::qualitative", n_colors), title =  paste("Hotspots", toupper(ip_x)))+
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
               outer.margins = c(0,0,0,0),
               inner.margins = c(0,0,0,0),
               legend.height= -0.2 )
   
   tmap_save(all_ip_map,
             filename= paste0(root, "_results/hotspots/ip_maps/full_ip_categories_map.png"),
             dpi=300, 
             #insets_tm= insetmap,
             #insets_vp=viewport(x= 0.44, y= 0.78, width= 0.2, height= 0.2),
             height=9,
             width=16,
             units="in")
   
 }
 
  
  for(k in  na.omit(unique(rast_labs$category)) ){
  

  sht <- na.omit(rast_labs$final_short_lab[rast_labs$category == k])
  
  repl <- rast_labs %>% 
    dplyr::filter(!is.na(category)) %>% 
    dplyr::select(ID, final_short_lab, long_labs, seq)
  
  
  
  cat_ids <- rast_labs %>% 
    dplyr::filter(category == k | grepl(sht, rast_labs$final_short_lab )) %>% 
    dplyr::select(ID, final_label, seq) 
  
  lst <- lapply(cat_ids$final_label, function(p){unlist(str_split(p, "\\+"))}) 

  
  for(j in 2:length(lst)){
    
    for(i in 1:length(lst[[j]])){
      
      to_test <- rast_labs %>% 
        dplyr::filter(final_short_lab == lst[[j]][i]) %>% 
        dplyr::pull(long_labs)
      
      if(any(to_test %in% unlist(lst[1:(j-1)])   )  ){
        lst[[j]][i] <-   lst[[j]][i]
      }else{
        lst[[j]][i] <- to_test
      }
      
    }
    
} 
 
  tmp_rst <- ht_rast
  tmp_rst[!ht_rast[] %in% cat_ids$ID]<- NaN
  
  
  cat_ids$final_label <- sapply(lst, function(l){paste(l, collapse = "+")})
  cat_ids$seq <- 1:nrow(cat_ids)
  
  cat_ids %>% 
    dplyr::select(ID, final_label) %>% 
    write_csv(., paste0(root, "_results/hotspots/ip_maps/", ip_x, "_hotspots_", k, "_label.csv"))
  
  ht_rast_f <- raster::subs(tmp_rst, cat_ids[, c("ID", "seq")])
  ht_rast_f <- as.factor(ht_rast_f)
  levels(ht_rast_f) <- data.frame(id = levels(ht_rast_f)[[1]], x = cat_ids$final_label, row.names = 1:length(cat_ids$final_label)   )
 
  n_colors <- length(cat_ids$final_label)
  #paletteer_d("wesanderson::Rushmore", n_colors)
  hots_map <- tmap::tm_shape(shp_c)+
    tm_borders(col = "gray50")+
    tm_shape(ht_rast_f)+
    tm_raster( style = "cat", palette = paletteer_dynamic("ggthemes_ptol::qualitative", n_colors), title =  paste("Hotspots", toupper(ip_x)))+
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
              outer.margins = c(0,0,0,0),
              inner.margins = c(0,0,0,0),
              legend.height= -0.2 )
  
  #x11();hots_map
  
 sg <-  tmaptools::bb_poly(shp_c)
  
  insetmap = tm_shape(shp) +
    tm_fill(col="lightgrey") +
    tm_shape(shp_c) + 
    tm_borders(lwd = 1, col="darkgrey") +
    tm_shape(sg) +
    tm_borders(lw=2, col="red") +
    tm_layout(inner.margins = c(0.04,0.04,0.04,0.04), outer.margins=c(0,0,0,0))
 
  if(var_name == "NAME_0" | length(unique(shp_c@data[, var_name  ])) > 2){
    tmap_save(hots_map,
              filename= paste0(root, "_results/hotspots/ip_maps/", ip_x, "_", k,"_map.png"),
              dpi=300, 
              #insets_tm= insetmap,
              #insets_vp=viewport(x= 0.44, y= 0.78, width= 0.2, height= 0.2),
              height=9,
              width=16,
              units="in")
    
  }else{
    tmap_save(hots_map,
              filename= paste0(root, "_results/hotspots/ip_maps/", ip_x, "_", k,"_map.png"),
              dpi=300, 
              insets_tm= insetmap,
              insets_vp=viewport(x= 0.44, y= 0.78, width= 0.2, height= 0.2),
              height=9,
              width=16,
              units="in")
  }
  
  
  ###############################################################
  ###### Intersection between hotspots and climate-conflict #####
  ###############################################################
  
  lapply(unique(conf_clust@data$label), function(var){
    
    
    tmp_df <- conf_clust
    tmp_df@data$label <-  as.character(tmp_df@data$label)
    tmp_df@data$intersect_conf_clim[tmp_df@data$label != var] <- "Others combinations"
    
    lvls <- tmp_df@data %>% 
      dplyr::select(clim_cluster_short_label, clim_cluster_order) %>% 
      dplyr::filter(!duplicated(clim_cluster_short_label)) %>% 
      dplyr::arrange(desc(clim_cluster_order)) %>% 
      pull(clim_cluster_short_label) %>% 
      paste0(var,"-[", .,"]") %>% 
      c(., "Others combinations" )
    
    tmp_df@data$intersect_conf_clim <- factor(tmp_df@data$intersect_conf_clim, levels = lvls)
    
    tmp_df <- tmp_df[!is.na(sp::over(tmp_df, shp_c, returnList = F)$GID_0), ] 
    tmp_df@data$intersect_conf_clim <- droplevels(tmp_df@data$intersect_conf_clim)
     
    
    main_col <- switch(var, "High conflict" = "YlOrRd",
                       "Moderate conflict" =  "Oranges",
                       "Limited conflict" = "YlGn") ##DFDFDF
    
    other_color <- "#DADADA"
    
    if(var == "Limited conflict"){
      pal_g <- brewer.pal(length(unique(tmp_df@data$intersect_conf_clim))-1 , main_col)
    }else{
      pal_g <- rev(brewer.pal(length(unique(tmp_df@data$intersect_conf_clim)) , main_col))
    }
    pal_g[length(pal_g)] <- other_color
    
    
    int_hotst_conf <- tmap::tm_shape(shp_c)+
      tm_borders(col = "black")+
      tm_shape(ht_rast_f)+
      tm_raster( style = "cat", palette = paletteer_dynamic("ggthemes_ptol::qualitative", n_colors), title =  paste("Hotspots", toupper(ip_x)))+
      tm_shape(tmp_df)+
      tm_fill(col = "intersect_conf_clim", palette = pal_g , alpha = 0.5, title = paste0(var, " cluster"))+
      #tm_borders(col ="black") + 
      tm_compass(type = "8star", position = c("right", "top")) +
      tm_scale_bar(breaks = c(0, 100, 200), text.size = 1, position = c(scale_bar_pos, scale_bar_top))+
      tm_layout(legend.outside=T, 
                legend.text.size = 1,
                legend.title.size= 1.3,
                legend.frame=F,
                #legend.position=c(0.985, 0.985),
                #legend.just = c("left", "top"), 
                #legend.width= 1,
                legend.outside.position = "right",
                legend.outside.size = 0.45,
                #outer.margins = c(0,0,0,0),
                #inner.margins = c(0,0,0,0)
                legend.height= -0.2 )
    #x11();int_hotst_conf
    
    if(var_name == "NAME_0"| length(unique(shp_c@data[, var_name  ])) > 2){
      tmap_save(int_hotst_conf,
                filename= paste0(root, "_results/cluster_results/",tolower(gsub(" ", "_", var)), "_",ip_x, "_", k , "_climate_intersection.png"),
                dpi=300, 
                #insets_tm= insetmap,
                #insets_vp= viewport(x= 0.45, y= 0.78, width= 0.2, height= 0.2),
                height=10,
                width=16,
                units="in")
    }else{
      
      tmap_save(int_hotst_conf,
                filename= paste0(root, "_results/cluster_results/",tolower(gsub(" ", "_", var)), "_",ip_x, "_", k , "_climate_intersection.png"),
                dpi=300, 
                insets_tm= insetmap,
                insets_vp= viewport(x= 0.45, y= 0.78, width= 0.2, height= 0.2),
                height=10,
                width=16,
                units="in")
    }
   
    
  })
  
  
  
 
  
}#end for ips categories
  
 
  
  }#end for


##################################################################################
############## conflict-climte-ip's overlays ####################################
################################################################################

#' Donwload from https://fews.net/fews-data/335
livelihood_pth <-  switch (iso,
                           "KEN" =  "livelihood/KE_LHZ_2011.shp",
                           "SEN" =  "livelihood/SN_LHZ_2021.shp",
                           "NGA" =  "livelihood/NG_LHZ_2018.shp",
                           "MLI" =  "livelihood/ML_LHZ_2014.shp",
                           "GTM" =  "livelihood/GT_LHZ_2016.shp",
                           "PHL" =  "livelihood/ML_LHZ_2014.shp",
                           "NER" =  "livelihood/NE_LHZ_2011.shp",
                           "BFA" =  "livelihood/BF_LHZ_2014.shp"
)

mf_diff <- raster::raster(paste0(root, "education/medn_difference_edu.tif"))

m_edu <- raster::raster(paste0(root, "education/medn_male_edu.tif"))

f_edu <- raster::raster(paste0(root, "education/medn_female_edu.tif"))

eth <- raster::shapefile("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/_global/ethnicity/GREG.shp") 

m_pop <- raster::raster(paste0(root, "gender_population/", iso,"_male_population.tif"))
m_pop[m_pop[] <= 1] <- NA

f_pop <- raster::raster(paste0(root, "gender_population/", iso,"_female_population.tif"))
f_pop[f_pop[] <= 1] <- NA


  
livelihoods <- raster::shapefile(paste0(root, livelihood_pth ))

fips_country <- switch (iso,
                        "KEN" = "KE",
                        "UGA" = "UG",
                        "GTM" = "GT",
                        "MLI" = "ML",
                        "NGA" = "NI",
                        "PHL" = "RP",
                        "SDN" = "SU",
                        "SEN" = "SG",
                        "ETH" = "ET",
                        "ZMB" = "ZA",
                        "ZWE" = "ZI",
                        "PHL" = "PH",
                        "NER" = "NE"
)

gwis_country <- switch (iso,
                        "KEN" = "501",
                        "UGA" = "500",
                        "GTM" = "90",
                        "MLI" = "432",
                        "NGA" = "475",
                        "PHL" = "840",
                        "SDN" = "625",
                        "SEN" = "433",
                        "ETH" = "530",
                        "ZMB" = "551",
                        "ZWE" = "552"
)



eth_c<- eth[eth@data$FIPS_CNTR == fips_country,]

eth_c@data$eth_short_name <- eth_c@data %>% 
  dplyr::select(contains("SHORTNAM")) %>%
  apply(., 1, function(rw){paste(na.omit(rw), collapse = ";")})

eth_c@data$eth_long_name <- eth_c@data %>% 
  dplyr::select(contains("LONGNAM")) %>%
  apply(., 1, function(rw){paste(na.omit(rw), collapse = ";")})

###### save conflict-climate - ip's intersection

#extract data from raster files

#clusts_to_share <- as(sf::st_read(paste0(root, "/_results/clim_conflict_ips_overlays.geojson")), "Spatial")

clusts_to_share <- conf_clust

rs <- sp::over(clusts_to_share, shp, returnList = T) 


mf_diff_ext <- exactextractr::exact_extract(mf_diff, sf::st_as_sf(clusts_to_share), fun = "median")
m_edu_ext <- exactextractr::exact_extract(m_edu, sf::st_as_sf(clusts_to_share), fun = "median")
f_edu_ext <- exactextractr::exact_extract(f_edu, sf::st_as_sf(clusts_to_share), fun = "median")
m_pop_ext <- exactextractr::exact_extract(m_pop, sf::st_as_sf(clusts_to_share), fun = "sum")
f_pop_ext <- exactextractr::exact_extract(f_pop, sf::st_as_sf(clusts_to_share), fun = "sum")
eth_ext <- sp::over(clusts_to_share, eth_c, returnList = T) 
liveext <- sp::over(clusts_to_share, livelihoods, returnList = T)


##extraer todas las variables de clima
all_rasts <- data.frame(var_name = list.files(paste0(root), pattern = ".tif$", recursive = T), 
           full_path = list.files(paste0(root), pattern = ".tif$", recursive = T, full.names = T) ) %>% 
  dplyr::filter(!grepl("old|_results", var_name))


clim_selc <- read_csv(paste0(root, "_results/cluster_results/climate/climate_most_imp_clim_vars.csv")) %>%
  unlist %>% 
  paste0(., collapse = "|")

clim_rasts <- all_rasts %>% 
  dplyr::filter(grepl( clim_selc, var_name)) %>% 
  purrr::pmap(., .f = function(var_name, full_path ){
    
    to_ret <- data.frame(vr = exactextractr::exact_extract(raster(full_path), sf::st_as_sf(clusts_to_share), fun = "median") )
     names(to_ret) <- stringr::str_extract(var_name, pattern = "[a-zA-Z0-9_]+.tif" ) %>% 
       stringr::str_replace(., pattern = ".tif", replacement = "") %>% 
       paste0("climvar_", .)
    return(to_ret)
    }) %>% 
  bind_cols(.)

#extraer todas las variables de los ips

ip <- "ip_all"

n_vars <- 10

ip_var_list <- read_csv(paste0(root,"_results/hotspots/soc_eco_all_variables.csv")) %>%
  dplyr::mutate(Code = ifelse(grepl("\\{iso\\}", Code), gsub("\\{iso\\}", iso, Code), Code),
                Code = tolower(Code))

ip_vars <- list.files(paste0(root,"_results/hotspots/"), pattern = paste0("_sorted_",unique(shp$NAME_0), '_', ip, ".xlsx"), full.names = T) %>%  
  readxl::read_excel(.) %>% 
  dplyr::mutate(Code = tolower(Code)) %>% 
  dplyr::left_join(., ip_var_list %>% 
                     dplyr::filter(IP_id == ip ) %>% 
                     dplyr::select(-Variable, - Classification), by = c("Code" = "Code")) %>% 
  dplyr::mutate(Code = ifelse(grepl("_awe", Code), paste0(iso,"_AWE"), Code),
                Code = ifelse(grepl("_rwi", Code), paste0(iso, "_rwi"), Code)) %>% 
  dplyr::slice(1:n_vars)


ip_rasts <- all_rasts %>% 
  dplyr::filter(grepl( ip_vars$Code %>% paste0(., collapse = "|"), var_name)) %>% 
  dplyr::mutate(Code =  stringr::str_extract(var_name, pattern = "[a-zA-Z0-9_]+.tif" ) %>% 
                  stringr::str_replace(., pattern = ".tif", replacement = "")) %>% 
  dplyr::left_join(., ip_vars) %>% 
  dplyr::mutate(final_name = paste0(Classification, "_", Code, "_thr[", Threshold, "]_perc[", Percentile,"]")) %>% 
  purrr::pmap(., .f = function(full_path, final_name, Threshold, Percentile, ...){
    
    r <- terra::rast(full_path)
    
    thr <- Threshold
    prc <- Percentile %>% as.numeric(.)
    
    if(!is.na(thr)){
      eval(parse(text = paste0('r[!(r ',thr,')] <- NA')))
    }
    
    mn <- median(r[], na.rm = T)
    qtl_country <- global(x = r, fun = quantile, probs = prc, na.rm = T) %>% as.numeric()
    
    c_mask <- terra::rast(c_mask)
    
    r <- r %>%
      terra::crop(x = ., y = terra::ext(c_mask)) %>% 
      terra::resample(x = ., y = c_mask) %>% 
      terra::mask(., mask = terra::vect(shp_c))
    
    to_check <- (min(r[], na.rm = T) == quantile(r[], probs = 0.25, na.rm = T) & prc == 0.1) |
      (max(r[], na.rm = T) == quantile(r[], probs = 0.75, na.rm = T) & prc == 0.9)
    
    stopifnot("Raster values distribution very Skewed" = !to_check)
    
    if(prc > 0.5){
      
      r[r < qtl_country] <- NA
     
      
    } else {
      
      r[r > qtl_country] <- NA
      
      
    }
    
    to_ret <- data.frame(vr = exact_extract(r, sf::st_as_sf(clusts_to_share), fun = "median") )
    names(to_ret) <- final_name
    return(to_ret)
    
  }) %>% 
  dplyr::bind_cols()


clusts_to_share@data <- clusts_to_share@data %>% 
  dplyr::mutate(NAME_1 = lapply(rs, function(df){df %>% pull(NAME_1) %>% unique(.) %>% paste(., collapse = ";")}) %>% unlist,
                NAME_2 = lapply(rs, function(df){df %>% pull(NAME_2) %>% unique(.) %>% paste(., collapse = ";")}) %>% unlist,
                NAME_3 = lapply(rs, function(df){df %>% pull(NAME_3) %>% unique(.) %>% paste(., collapse = ";")}) %>% unlist,
                livelihoods = lapply(liveext, function(df){df %>% pull(LZNAMEEN) %>% unique(.) %>% paste(., collapse = ";") }) %>%  unlist,
                median_male_female_edu_diff  = mf_diff_ext,
                median_male_edu = m_edu_ext,
                median_female_edu = f_edu_ext,
                male_population = m_pop_ext,
                female_population = f_pop_ext,
                ethnicity_short_name = lapply(eth_ext, function(df){df %>% pull(eth_short_name) %>% unique(.) %>% paste(., collapse =";")}) %>%  unlist,
                ethnicity_long_name = lapply(eth_ext, function(df){df %>% pull(eth_long_name)%>% unique(.) %>% paste(., collapse =";") }) %>%  unlist ) 


# extract information from ip variables
rc_fls <- list.files(path = paste0(root, "_results/hotspots/ip_maps/"), pattern = "rc.tif$")

ip_raw_info <- tibble(file_names = rc_fls,
                      source_file = stringr::str_extract(rc_fls , pattern = "_([A-Z](\\w| )+)_") %>% 
                        stringr::str_replace_all("_", "") %>% 
                        stringr::str_replace("ip", ""),
                      IP_id = stringr::str_extract(rc_fls , pattern = "ip_all"),
                      Code = str_replace(rc_fls , pattern = "_(ip_all[0-9]rc.tif)", "") %>% 
                        str_replace(. , pattern = "_[A-Z][a-z]+ [a-z]+", "") %>% 
                        str_replace(. , pattern = "_[A-Z][a-z]+", "") ,
                      file_path = paste0(root, "_results/hotspots/ip_maps/",rc_fls))

  
  

for(ip_x in get_ip_names){
  
  v_name <- paste0(ip_x, "_text_description")
  clusts_to_share@data <- clusts_to_share@data %>% 
    dplyr::mutate(!!v_name := ip_text_description(shp_object = clusts_to_share,
                                                  ip = ip_x,
                                                  df = ip_raw_info,
                                                  n_vars = 10 ))
  
}



clusts_to_share@data <- clusts_to_share@data %>% 
  dplyr::mutate(across(everything(.), as.vector)) %>% 
  dplyr::select(  EVENTS,
                  TYPE_RICHNESS = TYPE_RI,
                  SUBTYPE_RICHNESS = SUBTYPE,
                  ACTOR1_RICHNESS = ACTOR1_,
                  ACTOR2_RICHNESS = ACTOR2_,
                  FATALITIES = FATALIT,
                  conflict_clust_label = label, 
                  conflict_clust_short_label = clust_km, 
                  conflict_cluster_text_description = clutert_text_description,
                  clim_cluster_text_description = clim_cluster,
                  clim_cluster_short_label, 
                  intersect_conf_clim,
                  starts_with("ip"),
                  NAME_1,
                  NAME_2,
                  NAME_3,
                  clim_cluster_order,
                  starts_with("median"),
                  starts_with("ethnicity"),
                  female_population,
                  male_population,
                  livelihoods) 


clusts_to_share@data <- clusts_to_share@data %>% 
  dplyr::bind_cols(.,clim_rasts ) %>% 
  dplyr::bind_cols(., ip_rasts)
#clusts_to_share@data$conflict_cluster_text_description[1]

row.names(clusts_to_share@data) <- sapply(slot(clusts_to_share, "polygons"), function(x) slot(x, "ID"))

geojsonio::geojson_json(clusts_to_share) %>% 
  geojsonio::geojson_write(., file = paste0(root, "/_results/clim_conflict_ips_overlays.geojson") )

clusts_to_share@data %>% writexl::write_xlsx(paste0(root, "_results/clim_conflict_ips_overlays.xlsx"))

#raster::shapefile(clusts_to_share, paste0(root, "_results/clim_conflict_ips_overlays.shp"), overwrite = T)
base::saveRDS(clusts_to_share, paste0(root, "_results/clim_conflict_ips_overlays.rds"))


# 
# hot_high_conflict_areas <- clusts_to_share %>%
#   sf::st_as_sf() %>%
#   dplyr::filter(conflict_clust_short_label == "High",
#                 clim_cluster_order %in% c(3,4)) %>% 
#   dplyr::filter(if_any(.cols = contains("category"), .fns = ~ nchar(.) != 0  ))  
# 
# 
# hot_moderate_conflict_areas <- clusts_to_share %>%
#   sf::st_as_sf() %>%
#   dplyr::filter(conflict_clust_short_label == "Moderate",
#                 clim_cluster_order %in% c(3,4)) %>% 
#   dplyr::filter(if_any(.cols = contains("category"), .fns = ~ nchar(.) != 0  ))  
# 
# 
# 
# geojsonsf::sf_geojson(hot_high_conflict_areas) %>% 
#   geojsonio::geojson_write(., file = paste0(root, "_results/hot_areas_high_conflict.geojson") )
# 
# 
# geojsonsf::sf_geojson(hot_moderate_conflict_areas) %>% 
#   geojsonio::geojson_write(., file = paste0(root, "_results/hot_areas_moderate_conflict.geojson") )
# 

#raster::shapefile(as(hot_high_conflict_areas, "Spatial"), paste0(root, "_results/hot_areas_high_conflict.shp"), overwrite  = T)
#raster::shapefile(as(hot_moderate_conflict_areas, "Spatial"), paste0(root, "_results/hot_areas_moderate_conflict.shp"), overwrite  = T)
