# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain climate and conflict clusters at megapixels level of 20 km resolution
# Original data source:
# 1. Conflict: https://acleddata.com/curated-data-files/
# Steps:
# 1. Download manually the Africa's conflict table from ACLED
# 2. Execute this script to obtain:
#    Shapefile of conflict clusters at 20 km resolution
#    Shapefile of climate clusters at 20 km resolution
# Author: Andres Mendez, Harold Achicanoy
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
.rs.restartR()                      # Restart R session
g <- gc(reset = T); 
rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,raster,sf,stars,motif,tmap,spdep))
suppressMessages(pacman::p_load(meteo,sp,spacetime,gstat,plyr,xts,snowfall,doParallel,CAST,ranger))
suppressMessages(pacman::p_load(spatstat,maptools, Rcpp, maptree, exactextractr))
set.seed(1000)
source("./base__lowest_gadm.R")

#' Reclassifies raster layers using quantiles and also resamples to the same spatial resultion. 
#' @param root base directory path
#' @param iso country ISO code
#' @param country country name 
#' @param world_mask world raster layer in 5 km spatial resolution 
#' @param fconf filename of conflict excel file.


fconf <- 'Africa_1997-2022_Jul08.xlsx' #Name of conflict file

yearRange <- 1997:2022#range of years to select in ACCLED data
recompute <- TRUE #Recompute Kernel densities?
country_iso2 <- iso <- "SSD"
country <- 'South Sudan'
root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/'#dir path to folder data storage
baseDir <- paste0(root, "data/",country_iso2)


reclass_raster <- function(rast_path , shp_ext, world_mask, shp_country, dimension, conflict_area){
  
  r <- raster(rast_path)
  #r <- raster(file_paths$path[13] )
  #if(as.character(res(r)[1]) != as.character(res(world_mask)[1])){
  cat(paste("processign: ", rast_path, "\n"))
  w_msk<- world_mask %>% 
    raster::crop(., shp_ext) 
  
  #extent(world_mask) <- extent(world_mask)+5
  
  r <- r %>% 
    raster::resample(., w_msk) 
  
  if(dimension == "conflict"){
    c_area <- raster::shapefile(conflict_area)
    r <- r %>% 
      raster::crop(., extent(c_area)) %>% 
      raster::mask(., c_area)
    
    # if(grepl("FATALITIES|EVENTS", names(r))){
    #   pop_dens <- raster(paste0(baseDir, "/population_density/medn_popd.tif")) %>% 
    #     raster::resample(., w_msk) %>% 
    #     raster::crop(., extent(r))
    #   
    #   r <- r/pop_dens
    #   r[is.infinite(r[])] <- 0
    # }
  }
  
  #}
  
  if(length(unique(r[])) > 10){
    cat(rast_path, " reclassifying \n")
    qtl <- raster::quantile(r[r[] !=0], seq(.2,1,.2))
    rclmat <- data.frame(x1 = c(min(r[], na.rm = T),qtl[-length(qtl)]), x2 = qtl, y = 1:5) %>% as.matrix()
    ret <- r %>% 
      raster::reclassify(., rclmat) 
    
    
    ret <- stars::st_as_stars(ret)
    
  }else{
    cat(rast_path, " Not need for reclassify \n")
    ret <- r %>% 
      stars::st_as_stars()
    
  }
  
  return(ret)
}


get_conflic_data <- function(root, iso, country = country, world_mask, fconf, recompute){
  
  out <- paste0(root,'/data/',iso,'/conflict/',iso,'_conflict.csv')
  dir.create(path = dirname(out), F, T)
  
  if(recompute){
    # Filter African conflict to the specific country
    cnf <- readxl::read_excel(paste0(root,'/data/_global/conflict/', fconf), sheet = 1)
    cnf <- cnf %>% dplyr::filter(COUNTRY == country)
    cnf <- cnf %>% dplyr::filter(YEAR %in% yearRange)
    readr::write_csv(cnf, out)
    conflict <- cnf; rm(cnf, out)
 } else {
    # Load country conflict
    conflict <- readr::read_csv(out); rm(out)
  }
  
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
  conflict  <- cbind(conflict, raster::extract(x = shp, y = conflict[,c('LONGITUDE','LATITUDE')]))
  
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
  
  if(recompute){
    cat("Calculating conflict kernel density \n")
    
    ext <- extent(shp)
    msk <- world_mask %>% 
      raster::crop(., extent(shp)) %>% 
      raster::mask(., shp)
    
    p_var <- spatstat.geom::ppp(cnf_summ$x, cnf_summ$y, 
                                window = spatstat.geom::owin(c(ext[1], ext[2]), c(ext[3], ext[4]),
                                                             mask = matrix(TRUE,dim(msk)[1],dim(msk)[2]) ))
    ds <- spatstat.core::density.ppp(p_var, at = "pixels", kernel = "epanechnikov")
    knl <- raster::raster(ds) %>% 
      raster::crop(., extent(shp)) %>% 
      raster::mask(., shp) %>% 
      raster::resample(., raster(resolution = c(0.008983153, 0.008983153)))
    
    raster::writeRaster(knl, paste0(root,'/data/',iso,'/conflict/conflict_kernel_density.tif'), overwrite = T)
  }
  
  return(list(cnf_summ, sft))
}


#' Creates text lables for each cluster based on all variables in accled  
#' @param df dataframe of conflict data of the corresponding country



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
  
  
  lbls <- readxl::read_excel(paste0(root, "Hostpots_data_dictionary.xlsx")) %>% 
    dplyr::filter(Component == "Conflict") %>% 
    dplyr::select(definition = Variable, Code, Units) %>% 
    dplyr::mutate(Code = toupper(Code))
  
  
  to_label <- df %>% 
    dplyr::select(EVENTS:FATALITIES, starts_with("clust")) %>% 
    dplyr::group_by(clust) %>% 
    dplyr::summarise(across(EVENTS:FATALITIES, median)) %>% 
    dplyr::select(clust, EVENTS, FATALITIES, everything(.)) %>% 
    tidyr::pivot_longer(., cols = -clust, names_to = "variable", values_to = "median") %>% 
    dplyr::left_join(., 
                     lbls, by = c("variable" = "Code")) %>% 
    dplyr::mutate(id = 1:nrow(.) )
  
  
  
  cats <- c("High", "Moderate", "Limited")
  
  for(k in unique(to_label$variable)){
    
    tmp <- to_label %>% 
      dplyr::filter(variable == k) %>% 
      arrange(desc(median)) %>%
      dplyr::mutate(prefix = paste0("[", cats, "]")) %>% 
      dplyr::select(id, prefix)
    
    to_label[ tmp$id, "prefix"] <- tmp$prefix
    
  }
  
  
  txt_desc <- lapply(unique(to_label$clust), function(i){
    # start_text <- paste("Conflict cluster ", i , " is characterized by: ")
    ret <- to_label %>% 
      dplyr::filter(clust == i) %>% 
      dplyr::mutate(text = purrr::pmap(.l = list(def = definition, un = Units, pr = prefix, me = median), .f= function(def, un, pr, me){
        
        
        paste( pr, "values of", def, paste0("(", round(me, 2), " median ", un,")" ))
      }) %>%  unlist) %>% 
      pull(text) %>% 
      paste(., collapse = ", ") 
    
    
    return(data.frame(clust = i, clutert_text_description = ret))
  }) %>% 
    bind_rows() 
  
  ret <- conf_cluts_vals %>% 
    left_join(., txt_desc, by = c("clust" = "clust")) %>% 
    dplyr::mutate(short_label = factor(short_label, levels = c("High", "Moderate", "Limited")))
  
  
  
  return(ret)
}

#' Plot box plots of different clusters (high, median and  low) each conflict variable.
#' @param df Dataframe for
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

#' Compute statistics for 
#' @param df 
#' 
get_cluster_statistics <- function(df){
  
  df <- df %>% 
    dplyr::select(-contains("ov")) %>% 
    dplyr::mutate(knl = ifelse(is.na(knl), 0, knl))
  
  clust_median <- df %>%
    dplyr::group_by(clust) %>% 
    dplyr::summarise(across(where(is.numeric), median)) %>% 
    dplyr::ungroup()
  
  
  g_median <- df %>% 
    dplyr::summarise(across(where(is.numeric), median)) %>% 
    dplyr::mutate(clust = "Global") %>% 
    dplyr::select(clust, everything())
  
  rel_change <- lapply(
    c("EVENTS", "TYPE_RICHNESS" ,"SUBTYPE_RICHNESS", "ACTOR1_RICHNESS", "ACTOR2_RICHNESS" ,"FATALITIES" ,   "knl"),
    function(i){
      
      v <- g_median %>% dplyr::pull(i)
      
      ret <- ((clust_median[i] -  v)/v )*100
      return(ret)
    }
  ) %>% 
    dplyr::bind_cols() %>% 
    dplyr::rename_with(., ~paste0(.x, "_rel_change"))
  
  ret <- clust_median %>% 
    dplyr::rename_with(., ~paste0(.x, "_median"))%>% 
    dplyr::bind_cols(rel_change) %>% 
    dplyr::add_row(g_median %>% 
                     dplyr::rename_with(., ~paste0(.x, "_median" ))) 
  
  return(ret)    
}


world_mask <- raster::raster(paste0(root,"/data/_global/masks/mask_world_1km.tif")) %>% 
  raster::crop(., extent(shp))

#' Gets country conflict data
#' @param root
#' @param iso
#' @param country#
#'@param world_mask
#' 
conflict_raw <-  get_conflic_data(root = root,
                                  iso = iso,
                                  country = country,
                                  world_mask = world_mask, fconf, recompute) %>% 
  purrr::pluck(1)


#' Read country shapefile
shp <- raster::shapefile(paste0(baseDir,"/_shps/",country_iso2,".shp" )) %>% 
  sf::st_as_sf(.) %>% 
  dplyr::mutate(id = 1:nrow(.))
#country <- unique(shp$NAME_0)
#' Create a cluster to represent the megapixels
grd <- st_make_grid(st_bbox(extent(shp)+2), cellsize = 0.2, square =  T) %>% 
  st_as_sf(.) %>%
  dplyr::mutate(id = 1:nrow(.))


#######################################################
###### NEW CONFLICT CLUSTERING #######################
#####################################################

#' Output directory
out_conflict_dir<- paste0(root, "/data/",iso, "/_results/cluster_results/conflict")
if(!dir.exists(out_conflict_dir)){dir.create(out_conflict_dir, recursive = T)}
#' Population density
pop_dens <- raster::raster(paste0(root,"data/",iso, "/population_density/medn_popd.tif"))

shp <- raster::shapefile(paste0(root,"data/", iso, "/_shps/",iso,".shp" )) %>% 
  sf::st_as_sf() %>% 
  dplyr::mutate(id = 1:nrow(.))


knl <- raster::raster(paste0(root, "/data/", iso, "/conflict/conflict_kernel_density.tif"))
crs(knl) <- crs(world_mask)

coordinates(conflict_raw) <- ~x+y
crs(conflict_raw) <- crs(world_mask)
conflict_sf <- sf::st_as_sf(conflict_raw)

st_crs(grd) <- st_crs(conflict_sf)

grd <- grd %>% 
  dplyr::mutate(ov = sapply(st_intersects(grd, shp) , function(i){if(length(i)>0){return(1)}else{return(NA)}})) %>% 
  dplyr::filter(!is.na(ov)) %>% 
  dplyr::select(-ov) %>% 
  dplyr::mutate(id = 1:nrow(.))


#' Dataframe of aggregated values per megapixels including the kernel density results
conflict_sf <- conflict_sf %>% 
  dplyr::mutate(NA_intersect = sapply(st_intersects(conflict_sf, shp), function(i){ifelse(length(i)==0,NA, i)})) %>%
  dplyr::filter(!is.na(NA_intersect)) %>%
  dplyr::select(-NA_intersect) 

to_cluster <- conflict_sf %>%
  dplyr::mutate(ov = sapply(st_intersects(conflict_sf, grd), function(i){ifelse(length(i)==0,NA, i)}), 
                pop_dens =raster::extract(pop_dens, st_as_sf(.) %>% st_cast("POINT")%>% st_coordinates()),
                pop_dens = ifelse(is.na(pop_dens), 0, pop_dens),
                EVENTS = ifelse(pop_dens == 0 , 1, EVENTS/pop_dens),
                FATALITIES = ifelse(pop_dens == 0 , 1, FATALITIES/pop_dens),
                knl = raster::extract(knl, st_as_sf(.) %>% st_cast("POINT")%>% st_coordinates())) %>%
  dplyr::filter(!is.na(knl))%>%
  dplyr::group_by(ov) %>% 
  dplyr::summarise(EVENTS = median(EVENTS, na.rm = T),
                   TYPE_RICHNESS = max(TYPE_RICHNESS, na.rm = T),
                   SUBTYPE_RICHNESS = max(SUBTYPE_RICHNESS, na.rm = T),
                   ACTOR1_RICHNESS = max(ACTOR1_RICHNESS, na.rm = T),
                   ACTOR2_RICHNESS = max(ACTOR2_RICHNESS, na.rm = T),
                   FATALITIES = median(FATALITIES, na.rm = T)
                   ,knl = median(knl, na.rm = T)
  ) %>% 
  dplyr::ungroup()


#                knl = raster::extract(knl, st_as_sf(to_cluster) %>% st_cast("POINT")%>% st_coordinates())

x <- to_cluster %>% 
  dplyr::mutate(knl = ifelse(is.na(knl), 0, knl)) %>% 
  dplyr::select(-ov) %>% 
  st_drop_geometry() %>%  
  FactoMineR::PCA(., scale.unit = T, ncp =2, graph = F )

#' Identification of outliers based on 1.5  times the IQR and assigning them a lower
#' weight
cords <- rbind(x$ind$coord) %>% 
  as_tibble() %>% 
  dplyr::mutate(dst = sqrt(( Dim.1 - 0)^2 + (Dim.2 - 0 )^2   ),
                rng = dst > quantile(dst)[3] + IQR(dst)*1.5) 

thr <- cords %>% dplyr::filter(rng) %>% pull(dst) %>% quantile %>% .[4]

#' Assigning weights based distance outside the IQR
#' wparam weighting parameter
wparam <- 10
wh <- cords %>% 
  dplyr::mutate(wh = ifelse(dst > thr, wparam, 1)) %>% 
  dplyr::pull(wh)

#' Recalculate PCA with new weights
pca_w <- to_cluster %>% 
  dplyr::mutate(knl = ifelse(is.na(knl), 0, knl)) %>% 
  dplyr::select(-ov) %>% 
  st_drop_geometry() %>%  
  FactoMineR::PCA(., scale.unit = T, ncp =2, graph = F , row.w = wh)


x_wh <- rbind(pca_w$ind$coord)%>% 
  dist(., method = "euclidean")
set.seed(1000) #https://towardsdatascience.com/three-versions-of-k-means-cf939b65f4ea
km <- kmeans(x_wh, centers = 3, iter.max = 1000, nstart = 20, algorithm="MacQueen")#stats::kmeans(x_wh, 3)
table(km$cluster)

original_df <- conflict_sf %>% 
  dplyr::mutate(ov = sapply(st_intersects(conflict_sf, grd), function(i){ifelse(length(i)==0,NA, i)}),
                knl = raster::extract(knl, st_as_sf(.) %>% st_cast("POINT")%>% st_coordinates())) %>%
  dplyr::group_by(ov) %>% 
  dplyr::summarise(EVENTS = sum(EVENTS, na.rm = T),
                   TYPE_RICHNESS = max(TYPE_RICHNESS, na.rm = T),
                   SUBTYPE_RICHNESS = max(SUBTYPE_RICHNESS, na.rm = T),
                   ACTOR1_RICHNESS = max(ACTOR1_RICHNESS, na.rm = T),
                   ACTOR2_RICHNESS = max(ACTOR2_RICHNESS, na.rm = T),
                   FATALITIES = sum(FATALITIES, na.rm = T)
                   ,knl = median(knl, na.rm = T)
  ) %>% 
  dplyr::filter(!is.na(knl)) %>% 
  dplyr::ungroup() %>% 
  dplyr::bind_cols(clust = as.character(km$cluster))


cluster_labels <- original_df %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(EVENTS:FATALITIES, contains("clust")) %>% 
  get_cluster_labels()

plot(original_df['clust'])

write_csv(cluster_labels, paste0(root,"/data/", iso, "/_results/cluster_results/conflict/conflict_cluster_text_description.csv"))


to_save <- original_df %>% 
  st_drop_geometry() %>% 
  dplyr::right_join(., grd, by = c("ov" = "id")) %>% 
  dplyr::filter(!is.na(clust))


sf::st_write(to_save, paste0(root, "/data/", iso, "/_results/cluster_results/conflict/conflict_regular_clust.shp"), delete_dsn = T)
#st_intersects(grd, shp) %>% unlist %>% lenght

original_df <-  original_df %>% 
  left_join(., cluster_labels , by = c("clust")) %>% 
  dplyr::rename(clust_km = short_label) 


to_plot <- grd %>% 
  left_join(., original_df %>% st_drop_geometry(), by = c("id" = "ov")) %>% 
  dplyr::filter(!is.na(clust_km))

plot(to_plot['clust'])

mainmap3<- tmap::tm_shape(shp)+
  tm_borders(col = "black")+
  tm_shape(to_plot)+
  tm_fill(col = "clust_km", palette = c("#d7191c", "#e5a03e", "#ffffbf"), alpha = 0.7, title = expression("Conflict clusters"))#+
#tm_borders(col ="black") 


x11();mainmap3

tmap_save(mainmap3,
          filename= paste0(root, "/data/", iso, "/_results/cluster_results/conflict/conflict_regular_clust_map.png"),
          dpi=300, 
          #insets_tm=insetmap, 
          #insets_vp=vp,
          height=8,
          width=15,
          units="in")



to_boxplot <- to_save %>% 
  dplyr::select(-ov)

g <- make_cluster_plots(df = to_boxplot)
x11();g

ggsave(g, filename= paste0(root,"/data/", iso, "/_results/cluster_results/conflict/conflict_regular_clust_boxplots.png"),
       dpi=300, 
       height=8,
       width=15,
       units="in")


get_cluster_statistics(df = to_boxplot) %>% 
  write.csv(., paste0(root, "/data/", iso, "/_results/cluster_results/conflict/conflict_regular_clust_rel_change.csv"), row.names = F)
