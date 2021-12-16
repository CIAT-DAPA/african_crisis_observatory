require(pacman)
pacman::p_load(tidyverse, raster, sf, sp, stringr, stringi, tmap, FactoMineR)

get_conflic_data <- function(root, iso, country = 'Senegal'){
  
  out <- paste0(root,'/data/',iso,'/conflict/',iso,'_conflict.csv')
  dir.create(path = dirname(out), F, T)
  if(!file.exists(out)){
    # Filter African conflict to the specific country
    cnf <- readxl::read_excel(paste0(root,'/data/_global/conflict/Africa_1997-2021_Apr02.xlsx'), sheet = 1)
    cnf <- cnf %>% dplyr::filter(COUNTRY == country)
    readr::write_csv(cnf, out)
    conflict <- cnf; rm(cnf, out)
  } else {
    # Load country conflict
    conflict <- readr::read_csv(out); rm(out)
  }; rm(out)
  
  # Load the country lowest administrative level shapefile
  if(!file.exists(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))){
    dir.create(path = dirname(paste0(root,'/data/',iso,'/_shps/',iso,'.shp')), recursive = TRUE)
    shp <- lowest_gadm(iso = iso, out = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
  } else {
    shp <- raster::shapefile(x = paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
    adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
    shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
  }
  sft <- shp # Shapefile in terra format
  
  # Extract administrative names using reported conflict coordinates by ACLED
  conflict  <- cbind(conflict,raster::extract(x = shp, y = conflict[,c('LONGITUDE','LATITUDE')]))
  
  # Get conflict summarization variables per coordinates
  cnf_summ <- conflict %>%
    dplyr::group_by(LONGITUDE, LATITUDE) %>%
    dplyr::summarise(EVENTS           = dplyr::n(),
                     TYPE_RICHNESS    = EVENT_TYPE %>% unique() %>% length(),
                     SUBTYPE_RICHNESS = SUB_EVENT_TYPE %>% unique() %>% length(),
                     ACTOR1_RICHNESS  = ACTOR1 %>% unique() %>% length(),
                     ACTOR2_RICHNESS  = ACTOR2 %>% unique() %>% length(),
                     FATALITIES       = sum(FATALITIES)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(x= LONGITUDE, y = LATITUDE, everything(.))
  
  
  return(list(cnf_summ, sft))
}

get_cluster_labels <- function(df){
  
  conf_cluts_vals <- df %>% 
    dplyr::select(EVENTS:FATALITIES, starts_with("clust")) %>% 
    group_by(clust) %>% 
    dplyr::summarise(m_events = mean(EVENTS)) %>% 
    arrange(desc(m_events)) %>% 
    dplyr::mutate(label = c("High conflict", "Moderate conflict", "Limited conflict"),
                  short_label = c("High", "Moderate", "Limited")) %>% 
    dplyr::mutate(across(everything(.), as.character),
                  label = factor(label, levels = c("High conflict","Moderate conflict",  "Limited conflict")))
  
  return(conf_cluts_vals)
}

make_cluster_plots <- function(df){

 g <- df %>% 
   dplyr::select(EVENTS:FATALITIES, starts_with("clust")) %>% 
   tidyr::pivot_longer(-clust, names_to = "var", values_to = "vals") %>% 
    ggplot(aes(y= vals, x = clust))+
    geom_violin(trim=FALSE, fill="gray")+
    geom_boxplot(width=0.1, fill="white")+
    facet_wrap(~var, scales = "free")+
    xlab("")+
    theme_bw(base_size = 10)+ 
   scale_y_log10()
 
return(g)
}

root <- 'C:/Users/acmendez/OneDrive - CGIAR/African_Crisis_Observatory'


iso <- "KEN"

pop_dens <- raster::raster(paste0(root,"/data/",iso, "/population_density/medn_popd.tif"))

shp <- raster::shapefile(paste0(root,"/data/", iso, "/_shps/",iso,".shp" )) %>% 
  sf::st_as_sf() %>% 
  dplyr::mutate(id = 1:nrow(.))


world_mask <- raster::raster(paste0(root,"/data/_global/masks/mask_world_1km.tif")) %>% 
  raster::crop(., extent(shp))

knl <- raster::raster(paste0(root, "/data/", iso, "/conflict/conflict_kernel_density.tif"))
crs(knl) <- crs(world_mask)


conflict_raw <-  get_conflic_data(root = root,
                                             iso = iso) %>% 
  purrr::pluck(1)




coordinates(conflict_raw) <- ~x+y
crs(conflict_raw) <- crs(world_mask)
conflict_sf <- sf::st_as_sf(conflict_raw)

grd <- st_make_grid(st_bbox(extent(shp)), cellsize = 0.2, square =  T) %>% 
  st_as_sf(.) %>%
  dplyr::mutate(id = 1:nrow(.))

st_crs(grd) <- st_crs(conflict_sf)

grd <- grd %>% 
  dplyr::mutate(ov = sapply(st_intersects(grd, shp) , function(i){if(length(i)>0){return(1)}else{return(NA)}})) %>% 
  dplyr::filter(!is.na(ov)) %>% 
  dplyr::select(-ov) %>% 
  dplyr::mutate(id = 1:nrow(.))



to_cluster <- conflict_sf %>% 
  dplyr::mutate(ov = unlist(st_intersects(conflict_sf, grd)),
                pop_dens =raster::extract(pop_dens, st_as_sf(.) %>% st_cast("POINT")%>% st_coordinates()),
                EVENTS = ifelse(pop_dens == 0 , 1, EVENTS/pop_dens),
                knl = raster::extract(knl, st_as_sf(.) %>% st_cast("POINT")%>% st_coordinates())) %>%
  dplyr::group_by(ov) %>% 
  dplyr::summarise(EVENTS = sum(EVENTS, na.rm = T),
                   TYPE_RICHNESS = max(TYPE_RICHNESS, na.rm = T),
                   SUBTYPE_RICHNESS = max(SUBTYPE_RICHNESS, na.rm = T),
                   ACTOR1_RICHNESS = max(ACTOR1_RICHNESS, na.rm = T),
                   ACTOR2_RICHNESS = max(ACTOR2_RICHNESS, na.rm = T),
                   FATALITIES = sum(FATALITIES, na.rm = T),
                   knl = median(knl, na.rm = T)) %>% 
  dplyr::ungroup()


#                knl = raster::extract(knl, st_as_sf(to_cluster) %>% st_cast("POINT")%>% st_coordinates())
x <- to_cluster %>% 
  dplyr::mutate(knl = ifelse(is.na(knl), mean(knl, na.rm = T), knl)) %>% 
  dplyr::select(-ov) %>% 
  st_drop_geometry() %>%  
  FactoMineR::PCA(., scale.unit = T, ncp =2 )

x <- rbind(x$ind$coord)%>% 
  dist(., method = "euclidean")


dst<- dist(x, method = "euclidean")
hc <- hclust(dst, method = "complete") %>% 
  stats::cutree(tree = ., k = 3)

km <- stats::kmeans(x, 3)



to_cluster$clust_km <-  left_join( tibble(clust = as.character(km$cluster)), to_cluster %>% 
                                      st_drop_geometry() %>% 
                                      dplyr::select(EVENTS:FATALITIES) %>% 
                                      bind_cols(clust = km$cluster) %>% 
                                      get_cluster_labels() %>% 
                                      dplyr::select(clust, label) , by = c("clust"))%>% 
  pull(label)



to_plot <- grd %>% 
  left_join(., to_cluster %>% st_drop_geometry(), by = c("id" = "ov")) %>% 
  dplyr::filter(!is.na(clust_km))


sf::st_write(to_plot, paste0(root, "/data/", iso, "/_results/cluster_results/conflict/conflict_regular_clust.shp"), delete_dsn = T)
#st_intersects(grd, shp) %>% unlist %>% lenght


mainmap3<- tmap::tm_shape(shp)+
  tm_borders(col = "black")+
  tm_shape(to_plot)+
  tm_fill(col = "clust_km", palette = c("#d7191c", "#e5a03e", "#ffffbf"), alpha = 0.7, title = expression("Conflict clusters"))+
  tm_borders(col ="black") 


x11();mainmap3

x11();make_cluster_plots(df = to_cluster %>% dplyr::select(EVENTS:FATALITIES, clust = clust_km)%>% st_drop_geometry())


  tm_fill(col = "clust", palette = c("#d7191c", "#e5a03e", "#ffffbf"), alpha = 0.7, title = expression("Conflict clusters"))+
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



tmap_save(mainmap,
          filename= paste0(root, "_results/cluster_results/conflict/geographic_distr_conflict.png"),
          dpi=300, 
          #insets_tm=insetmap, 
          #insets_vp=vp,
          height=8,
          width=15,
          units="in")




