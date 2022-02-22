# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain maps, graphs and final tables
# Author: Andres Mendez
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, readxl, writexl, raster,terra, sp, sf, stringr, stringi, lattice, rasterVis, maptools,
                                latticeExtra, RColorBrewer, tmap, geojson, geojsonio, MetBrewer))

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

baseDir <- "C:/Users/acmendez/OneDrive - CGIAR/African_Crisis_Observatory/data/"
w_mask <- raster::raster(paste0(baseDir, "_global/masks/mask_world_1km.tif"))

iso <- "KEN"

root <- paste0(baseDir, iso, "/")

to_share_dir <- paste0(root, "_results/", iso, "_to_share")
if(!dir.exists(to_share_dir)){dir.create(to_share_dir)}

scale_bar_pos <- switch( iso, "ZWE" = "left", "KEN" = "left", "UGA" = "right", "MLI" = "left", "SEN" = "left", "NGA" = "right", "SDN" = "right")
scale_bar_top <- switch( iso, "ZWE" = "bottom", "KEN" = "bottom", "UGA" = "bottom", "MLI" = "bottom", "SEN" = "top", "NGA" = "bottom", "SDN" = "bottom")

##################################################################
######### Generate climate and conflict intersection #############
##################################################################

shp_c  <- raster::shapefile(paste0(root, "_shps/", iso, ".shp"))

c_mask <- w_mask %>% 
  raster::crop(., extent(shp_c)) %>% 
  raster::mask(., shp_c)

conf_clust <- raster::shapefile(paste0(root, "_results/cluster_results/conflict/conflict_regular_clust.shp"))
clim_clust <- raster::shapefile(paste0(root, "_results/cluster_results/climate/climate_regular_clust.shp"))

crs(conf_clust) <- crs(c_mask)
crs(clim_clust) <- crs(c_mask)

conf_clust@data <- conf_clust@data %>% 
  dplyr::rename("label" = clst_km) %>% 
  dplyr::mutate(short_label = stringr::str_extract(string = label, pattern = "[A-Za-z]+"),
                label = factor(label, levels = c("High conflict", "Moderate conflict", "Limited conflict"))) 

# conf_clust_mts<- readxl::read_excel(paste0(root, "_results/cluster_results/conflict/conflict_cluster_summary_metrics.xlsx"), 
#                                     sheet = "reg_rel_change")

# conf_cluts_labs <- conf_clust_mts %>% 
#   dplyr::slice(grep("kernel", Variables)) %>% 
#   dplyr::select(paste0("clust_", 1:3)) %>% 
#   tidyr::pivot_longer(cols = everything(.), names_to = "clust", values_to = "vals") %>% 
#   arrange(desc(vals)) %>% 
#   dplyr::mutate(label = c("High conflict", "Moderate conflict", "Limited conflict"),
#                 short_label = c("High", "Moderate", "Limited"),
#                 clust_num = stringr::str_extract(clust, "[0-9]")) %>% 
#   dplyr::mutate(across(everything(.), as.character),
#                 label = factor(label, levels = c("High conflict","Moderate conflict",  "Limited conflict")))
# 
# 
# 
# 
# conf_clust@data <- conf_clust@data %>%
#   dplyr::mutate(clust = as.character(clust)) %>% 
#   dplyr::left_join(., conf_cluts_labs %>% dplyr::select(label, short_label, clust_num), by = c("clust" = "clust_num"))

clim_clust_labs <- readxl::read_excel(paste0(baseDir, "temp_climate_clusters_labels.xlsx")) %>% 
  dplyr::filter(Country == iso) %>% 
  dplyr::mutate(across(everything(.), as.character))

clim_clust@data <- clim_clust@data %>% 
  dplyr::mutate(clust = as.character(clust)) %>% 
  dplyr::left_join(., clim_clust_labs %>% dplyr::select(clust = Cluster, label = Label), by = c("clust" = "clust"))

conf_clust@data$clim_cluster <- as.character(sp::over(conf_clust, clim_clust, returnList = F)$label)

conf_clust@data$intersect_conf_clim <- paste0(conf_clust@data$label, "-", conf_clust@data$clim_cluster)

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
conf_occ@data$over <- conf_occ %over%  shp_c %>% dplyr::pull(GID_0)
conf_occ <- conf_occ[!is.na(conf_occ$over),]
#x11()
mainmap <- tmap::tm_shape(shp_c)+
  tm_borders(col = "black")+
  tm_shape(conf_clust)+
  tm_fill(col = "label", palette = c("#d7191c", "#e5a03e", "#ffffbf"), alpha = 0.7, title = expression("Conflict clusters"))+
  tm_borders(col ="black") + 
  tm_compass(type = "8star", position = c("right", "top")) +
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

write_csv(events_bp %>% 
            dplyr::group_by(EVENT_TYPE) %>% 
            dplyr::summarise(short_label, n, freq = prop.table(n)), paste0(to_share_dir, "/EVENTS_TYPE_conflict_barplot_df.csv"))

g1 <- events_bp %>% 
  dplyr::mutate(short_label  = factor(short_label , levels = c("High conflict",  "Moderate conflict", "Limited conflict") )) %>%
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

write_csv(fata_bp %>% 
            dplyr::group_by(EVENT_TYPE) %>% 
            dplyr::summarise(short_label, counts, freq = prop.table(counts)), paste0(to_share_dir, "/FATALITIES_conflict_barplot_df.csv"))

g2 <- fata_bp %>% 
  dplyr::mutate(short_label  = factor(short_label , levels = c("High",  "Moderate", "Limited") )) %>% 
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

######################################################################################
################## Create intersection of climate-conflict map #######################
######################################################################################

conf_clust@data <- conf_clust@data %>% 
  dplyr::mutate(inter_short_label = case_when(
    grepl("High", intersect_conf_clim) & grepl("Harsh", intersect_conf_clim) ~ "High conflict-Harsh climate",
    grepl("Moderate", intersect_conf_clim) & grepl("Harsh", intersect_conf_clim)~ "Moderated conflict - Harsh climate",
    grepl("Limited", intersect_conf_clim) & grepl("Good", intersect_conf_clim) ~ "Limited conflict-Good climate",
    TRUE ~ "Other combinations"
    
  ),
  inter_short_label = factor(inter_short_label, levels = c("High conflict-Harsh climate",
                                                           "Moderated conflict - Harsh climate",
                                                           "Limited conflict-Good climate",
                                                           "Other combinations"
                                                           )))

raster::shapefile(conf_clust, paste0(root, "_results/cluster_results/conflict_climate_intersection.shp"), overwrite = T)

clusts_to_share <- conf_clust
clusts_to_share@data <- clusts_to_share@data %>% 
  dplyr::select(id,  conflict_clust = label, clim_cluster, intersect_conf_clim)
row.names(clusts_to_share@data) <- sapply(slot(clusts_to_share, "polygons"), function(x) slot(x, "ID"))

geojsonio::geojson_json(clusts_to_share) %>% 
  geojsonio::geojson_write(county_json_clipped, file = paste0(to_share_dir, "/", iso, "_conflict_climate_clusters.geojson") )

#x11()
inter_map <- tmap::tm_shape(shp_c)+
  tm_borders(col = "black")+
  tm_shape(conf_clust)+
  tm_fill(col = "inter_short_label", palette = c("#d7191c", "#fec980", "#2F740F", "#DFDFDF"), alpha = 0.7, title = expression("Conflict clusters"))+
  tm_borders(col ="black") + 
  tm_compass(type = "8star", position = c("right", "top")) +
  tm_scale_bar(breaks = c(0, 100, 200), text.size = 1, position = c(scale_bar_pos, scale_bar_top))+
  tm_layout(legend.outside=T, 
            legend.text.size = 1.3,
            legend.title.size= 1.3,
            legend.frame=F, 
            #legend.position=c(0.985, 0.985),
            legend.just = c("left", "top"), 
            #legend.width=-0.25,
            #outer.margins = c(0,0,0,0),
            #inner.margins = c(0,0,0,0)
            legend.height= -0.2 )
x11();inter_map

tmap_save(inter_map,
          filename= paste0(root, "_results/cluster_results/conflict_climate_intersection.png"),
          dpi=300, 
          #insets_tm=insetmap, 
          #insets_vp=vp,
          height=9,
          width=16,
          units="in")

#################################################
###### Hotspots map for IP1 #####################
#################################################

get_ip_names <- list.files(paste0(root, "_results/hotspots/")) %>% 
  str_extract(.,"ip[0-9]") %>% 
  unique() %>% 
  na.omit()

for(i in get_ip_names){
  
  ip_x <- i
  
  ip_codes <- read_csv(paste0(root, "_results/hotspots/", iso, "_hotspots_values.csv")) %>% 
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
    }) )
  
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
  
  ht_rast <- raster::raster(paste0(root,"_results/hotspots/",iso, '_all_cat_hotspots_', ip_x, ".tif" )) %>% 
    raster::crop(., extent(conf_clust)) %>% 
    raster::mask(., conf_clust)
  
  writeRaster(ht_rast, paste0(to_share_dir, "/", ip_x, "_hotspots_map.tif"), overwrite = T)
  
  rast_labs <- tibble( ID = unique(ht_rast[]) ) %>% 
    tidyr::drop_na() %>% 
    dplyr::left_join(., final_label_tbl, by = c("ID" = "rast_values")) %>% 
    dplyr::mutate(final_label = fix_label(rast_labs = ., labs =  final_label_tbl %>% 
                                            dplyr::filter(!is.na(category)) ),
                  chars = nchar(final_short_lab),
                  seq = 1:nrow(.)) %>% 
    arrange(chars) %>% 
    dplyr::select(ID, final_label, seq)
  
  rast_labs %>% 
    dplyr::select(value_ID = ID, label = final_label) %>% 
    write_csv(., paste0(to_share_dir, "/", ip_x, "_hotspots_labels.csv"))
  
  ht_rast_f <- raster::subs(ht_rast, rast_labs[, c("ID", "seq")])
  ht_rast_f <- as.factor(ht_rast_f)
  levels(ht_rast_f) <- data.frame(id = levels(ht_rast_f)[[1]], x =   rast_labs$final_label   )
  
  hots_map <- tmap::tm_shape(shp_c)+
    tm_borders(col = "gray50")+
    tm_shape(ht_rast_f)+
    tm_raster( style = "cat", palette = "-magma", title =  paste("Hotspots", toupper(ip_x)))+
    tm_compass(type = "8star", position = c("right", "top")) +
    tm_scale_bar(breaks = c(0, 100, 200), text.size = 1, position = c(scale_bar_pos, scale_bar_top))+
    tm_layout(legend.outside=T, 
              legend.text.size = 1.3,
              legend.title.size= 1.3,
              legend.frame=F, 
              #legend.position=c(0.985, 0.985),
              legend.just = c("left", "top"), 
              #legend.width=-0.25,
              #outer.margins = c(0,0,0,0),
              #inner.margins = c(0,0,0,0)
              legend.height= -0.2 )
  
  tmap_save(hots_map,
            filename= paste0(root, "_results/hotspots/", ip_x, "_map.png"),
            dpi=300, 
            #insets_tm=insetmap, 
            #insets_vp=vp,
            height=9,
            width=16,
            units="in")
  
  ###############################################################
  ###### Intersection between hotspots and climate-conflict #####
  ###############################################################
  
  #conf_clust_f <- conf_clust[conf_clust@data$inter_short_label == "High/Moderate conflict-Harsh climate",]
  
  int_hotst_conf <- tmap::tm_shape(shp_c)+
    tm_borders(col = "gray50")+
    tm_shape(ht_rast_f)+
    tm_raster( style = "cat", palette = "-magma", title =  paste("Hotspots", toupper(ip_x)))+
    tm_shape(conf_clust)+
    tm_fill(col = "inter_short_label", palette = c("#d7191c", "#fec980", "#2F740F", "#DFDFDF"), alpha = 0.5, title = expression("Conflict clusters"))+
    tm_borders(col ="black") + 
    tm_compass(type = "8star", position = c("right", "top")) +
    tm_scale_bar(breaks = c(0, 100, 200), text.size = 1, position = c(scale_bar_pos, scale_bar_top))+
    tm_layout(legend.outside=T, 
              legend.text.size = 1.3,
              legend.title.size=1.3,
              legend.frame=F, 
              #legend.position=c(0.985, 0.985),
              legend.just = c("left", "top"), 
              #legend.width=-0.25,
              #outer.margins = c(0,0,0,0),
              #inner.margins = c(0,0,0,0)
              legend.height= -0.2 )
  
  tmap_save(int_hotst_conf,
            filename= paste0(root, "_results/intersection_hotspots_conf_clim_", ip_x, "_map.png"),
            dpi=300, 
            #insets_tm=insetmap, 
            #insets_vp=vp,
            height=9,
            width=16,
            units="in")
  
  }#end for
