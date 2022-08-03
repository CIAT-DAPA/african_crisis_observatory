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
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,raster,sf,stars,motif,tmap))
suppressMessages(pacman::p_load(meteo,sp,spacetime,gstat,plyr,xts,snowfall,doParallel,CAST,ranger))
suppressMessages(pacman::p_load(spatstat,maptools, Rcpp, maptree, exactextractr))


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

get_conflic_data <- function(root, iso, country = 'Senegal', world_mask){
  
  out <- paste0(root,'/data/',iso,'/conflict/',iso,'_conflict.csv')
  dir.create(path = dirname(out), F, T)
  if(!file.exists(out)){
    # Filter African conflict to the specific country
    cnf <- readxl::read_excel(paste0(root,'/data/_global/conflict/Africa_1997-2022_Jul08.xlsx'), sheet = 1)
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
  
  if(!file.exists(paste0(root,'/data/',iso,'/conflict/conflict_kernel_density.tif'))){
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
    
    writeRaster(knl, paste0(root,'/data/',iso,'/conflict/conflict_kernel_density.tif'), overwrite = T)
  }
  
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


get_sum_cl_mtrs <- function(rast_paths, eco_grid_sf,shp_ext, world_mask){
  rs <- lapply(rast_paths, function(i){
    
    r <- raster(i) %>% 
      raster::resample(., world_mask %>% raster::crop(., shp_ext)) %>% 
      raster::crop(., shp_ext)
    
  }) %>% raster::stack(.)
  
  vals <- exactextractr::exact_extract(rs, eco_grid_sf, full_colnames = T)
  
  df_extracted <- eco_grid_sf %>% 
    mutate(rast_vals = purrr::map(.x = vals , .f = function(.x){return(.x)})) %>% 
    dplyr::select(clust, rast_vals) %>% 
    sf::st_drop_geometry()
  
  df_extracted_final <- lapply(unique(df_extracted$clust), function(i){
    ret <- df_extracted %>% 
      dplyr::filter(clust == i) %>% 
      dplyr::pull(rast_vals) %>% 
      bind_rows() %>% 
      dplyr::mutate(coverage_fraction = i) %>% 
      dplyr::rename("clust" = coverage_fraction) 
    
    return(ret)
  }) %>% 
    bind_rows() %>% 
    as_tibble
  return(df_extracted_final)
}

clust_descriptives <- function(clust_sum){
  
  results <- list()
  
  normal_vars <- clust_sum[, !grepl("(cv)|(cvar)|(trnd)", tolower(names(clust_sum))) ]
  cv_vars <- clust_sum[, grepl("(cv)|(cvar)", tolower(names(clust_sum)))]
  nms_no_cv <- names(clust_sum)[!grepl("(cv)|(cvar)", tolower(names(clust_sum)))]
  trnd_vars <- clust_sum[, nms_no_cv[grepl("trnd",nms_no_cv)]]
  rm(nms_no_cv)
  
  global_median <- normal_vars %>% 
    dplyr::select(-starts_with("clust")) %>% 
    apply(., 2, function(i){
      ret <- median(i, na.rm =T )
      return(ret)
    })
  
  
  normal_median <- normal_vars %>% 
    dplyr::group_by(clust) %>% 
    dplyr::summarise(across(everything(), function(i){
      ret <- median(i, na.rm = T)
      return(ret)
    }))%>%
    ungroup()  
  
  results$median_rel_change <- cbind(clust =normal_median$clust, t((t(normal_median[,-1]) - (global_median))/global_median)*100) %>% 
    as_tibble()
  
  results$cv_rel_change <- cv_vars %>% 
    bind_cols(., clust = clust_sum$clust) %>% 
    dplyr::group_by(clust) %>% 
    dplyr::summarise(across(everything(), function(i){
      ret <- ((median(i, na.rm = T) - 0.15)/0.15)*100
      return(ret)
    }))%>%
    ungroup() 
  
  results$trnd_rel_change <- trnd_vars %>% 
    bind_cols(., clust = clust_sum$clust) %>% 
    dplyr::group_by(clust) %>% 
    dplyr::summarise(across(everything(), function(i){
      ret <- ((median(i, na.rm = T) - 0.1)/0.1)*100
      return(ret)
    }))%>%
    ungroup() 
  
 ret <- purrr::reduce(results, left_join, by  = "clust")
 ret$clust <- paste0("cluster_", ret$clust)
 global_metrics <- c("global_metrics",global_median, rep(0.15, ncol(results$cv_rel_change)-1), rep(0.1, ncol(results$trnd_rel_change)-1))
 
 ret <- rbind(ret, global_metrics) %>% 
     tidyr::pivot_longer(., cols = -c("clust"), names_to = "variable", values_to ="value") %>% 
   tidyr::pivot_wider(names_from = "clust", values_from = "value")
 
  return(ret)
}


labeling_function <- function(db, n_vars){
  
  
  res_pca <-  db %>% 
    dplyr::select(-contains("flood")) %>% 
    #dplyr::filter(clust == 1) %>% 
    dplyr::select(-clust) %>% 
    FactoMineR::PCA(., scale.unit = T, ncp = 3, graph = F)
  
  #plot(res_pca, choix = "varcor" )
  
  selected_vars <- res_pca$var$contrib %>% 
    tibble::as_tibble( ., rownames = "var_name") %>%
    dplyr::mutate(Dim.1_w = Dim.1*res_pca$eig[1,2]/100, 
                  Dim.2_w = Dim.2*res_pca$eig[2,2]/100,
                  idx = Dim.1_w + Dim.2_w,
                  idx2 = idx/max(idx)) %>% 
    dplyr::arrange(., desc(idx2))  %>% 
    dplyr::slice(1:n_vars) %>% 
    pull(var_name)
  
  
  glb_df <- read_csv(paste0(dest_dir, dimension, "_reg_cluster_statistics.csv"))%>% 
    dplyr::filter(variable %in% selected_vars)
  
  
  lbls <- readxl::read_excel(paste0(root, "Hostpots_data_dictionary.xlsx")) %>% 
    dplyr::filter(Component == "Climate") %>% 
    dplyr::select(definition = Variable, Code, Units) %>% 
    dplyr::filter(Code %in% selected_vars)
  
  
  
  cats <- list( n_2 = c("High", "Low"),
                n_3 = c("High", "Moderate","Low"),
                n_4 = c("High", "High-Moderate", "Moderate-Low", "Low"),
                n_5 = c("Very High", "High", "Moderate", "Low", "Very Low"),
                n_6 = c("Very High", "High", "High-Moderate", "Moderate-Low", "Low", "Very Low"),
                n_7 = c("Very High", "High", "High-Moderate", "Moderate", "Moderate-Low", "Low", "Very Low"))
  
  fr <- dplyr::left_join(glb_df, lbls, by = c('variable' = 'Code')) %>% 
    dplyr::mutate(id = 1:nrow(.))
  
  n_clust <- length(unique(fr$clust))
  
  for(i in unique(fr$variable)){
    
    tmp <- fr %>% 
      dplyr::filter(variable == i) %>% 
      arrange(desc(median)) %>% 
      dplyr::mutate(prefix = unlist(cats[paste0("n_", n_clust)])) %>% 
      dplyr::select(id, prefix)
    
    
    
    fr[ tmp$id, "prefix"] <- tmp$prefix
  }
  
  
  
  
  gen_text <- function(df_clust){
    
    clust_num <- unique(df_clust$clust)
    
    
    txt_start <- paste("Climate cluster ", clust_num, "is characterized by: ")
    
    ret <- df_clust %>% 
      dplyr::mutate(text = purrr::pmap(.l = list(def = definition, un = Units, pr = prefix, me = median), .f= function(def, un, pr, me){
        
        paste( paste0("[", pr, "]"), "values of", def, paste0("(", round(me, 2), " ", un,")" ))
      }) %>%  unlist) %>% 
      pull(text) %>% 
      paste(., collapse = ", ")
    
    ret <- paste(txt_start, ret)
    return(ret)
    
  }
  
  text_output <- sapply(unique(fr$clust), function(i){
    
    fr %>% 
      dplyr::filter(clust == i) %>% 
      gen_text(.)
  }, simplify = T)
  
  
  ret <- tibble(clust = unique(fr$clust), text_output )
  
  return(ret)
}#end function





root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/'#dir path to folder data storage
country_iso2 <- iso <- "KEN"

baseDir <- paste0(root, "data/",country_iso2)

shp <- raster::shapefile(paste0(baseDir,"/_shps/",country_iso2,".shp" )) %>% 
  sf::st_as_sf() %>% 
  dplyr::mutate(id = 1:nrow(.))

country <- unique(shp$NAME_0)

source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/hotspots/01_link_IPinfo_climate_clusters.R') # Link IP text to identify climate variables



clm <- select_clim_vars(root = substr(root, start = 1, stop = nchar(root)-1 ), 
                        iso  = iso, 
                        cntr = country) %>% 
  dplyr::pull(Code) %>% 
  unique()


fls <- list.files(path = paste0(root,'/data/',country_iso2), pattern = 'tif$', full.names = T, recursive = T)
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')
fls <- fls[unlist(grep2(pattern = clm, x = fls))]; rm(clm)
fls <- unique(fls)
fls <- as.character(na.omit(fls))



file_paths <- tibble(path = fls,
                     type = "climate" ) %>% 
  add_row(path = c(list.files(paste0(baseDir, '/conflict'), pattern = ".tif$", full.names = T)), type = "conflict")


check_files <- lapply(file_paths$path, file.exists) %>% unlist()

world_mask <- raster(paste0(root, "/data/_global/masks/mask_world_1km.tif"))


pop_dens <- raster::raster(paste0(root,"/data/",iso, "/population_density/medn_popd.tif"))


knl <- raster::raster(paste0(root, "/data/", iso, "/conflict/conflict_kernel_density.tif"))
crs(knl) <- crs(world_mask)


grd <- st_make_grid(st_bbox(extent(shp)+2), cellsize = 0.2, square =  T) %>% 
  st_as_sf(.) %>%
  dplyr::mutate(id = 1:nrow(.))



stopifnot("File not found in paths. " = all(check_files))

### load and reclassify all rasters

dimension <- "climate"

cat(">>> starting process for: ", dimension, "\n")


dest_dir <- paste0(baseDir, "/_results/cluster_results/", dimension, "/")


if(!dir.exists(dest_dir)){dir.create(dest_dir, recursive = T)}

#file_paths %>% filter(type == dimension)  %>% pull(path)
r_files <- lapply(file_paths %>% filter(type == dimension)  %>% pull(path),
                  reclass_raster, 
                  shp_ext = extent(grd), 
                  world_mask= world_mask,
                  dimension  = dimension,
                  conflict_area = paste0(baseDir, "/conflict/conflict_area.shp"))

r_names <- lapply(file_paths %>% filter(type == dimension)  %>% pull(path), function(i){names(raster(i))})
names(r_files) <- unlist(r_names)

# =========================================================== #
# Spatial patterns clustering
# =========================================================== #
cat(">>> Calculating regular clusters \n")

eco_data <- r_files %>% purrr::reduce(., c, try_hard = T)


eco_signature <- motif::lsp_signature(eco_data, type = "incove", window = grd["id"])

eco_dist      <- motif::lsp_to_dist(eco_signature, dist_fun = "jensen-shannon")

eco_hclust <- hclust(eco_dist, method = "ward.D2")

###optimal number of cluster following 
###Kelley, L.A., Gardner, S.P., Sutcliffe, M.J. (1996) An automated approach for clustering an ensemble of NMR-derived protein structures into conformationally-related subfamilies, Protein Engineering, 9, 1063-1065.


if(dimension == "conflict"){
  c_optim_num <- 3
}else{
  optcl2 <- maptree::kgs(cluster = eco_hclust, diss = eco_dist, maxclust = 10)
  
  c_optim_num <- as.numeric(names(optcl2[which(optcl2 == min(optcl2))]))
  
}

cat(">>> Optim number of cluster is:",c_optim_num ,"/n")

clusters <- cutree(eco_hclust, k = c_optim_num)

ext <- c(xmin = extent(shp)[1], ymin = extent(shp)[3], xmax = extent(shp)[2], ymax = extent(shp)[4] )

eco_grid_sf <- motif::lsp_add_clusters(eco_signature, clusters, window = grd["id"])


#sf::st_write(eco_grid_sf, "D:/OneDrive - CGIAR/Attachments/Desktop/mapas/motfi_tr_conf_x.shp", delete_dsn = T)


metrics_clust <- motif::lsp_add_quality(eco_grid_sf, eco_dist, type = "cluster") %>%
  dplyr::group_by(clust) %>%
  dplyr::summarise(inhomogeneity = mean(inhomogeneity),
                   distinction = mean(distinction),
                   quality = mean(quality)) %>% 
  sf::st_drop_geometry() %>% 
  dplyr::mutate(cluster_type = "regular")

# # =========================================================== #
# # Irregular clusterin patters
# # =========================================================== #
# 
# cat(">>> Calculating irregular clusters /n")
# 
# irr_clust_lsp <- lsp_signature(eco_data, type = "incove",
#                                window = shp["id"], normalization = "pdf") 
# 
# irr_clust_dist <- lsp_to_dist(irr_clust_lsp, dist_fun = "jensen-shannon")
# irr_clust_dist_hc <-  hclust(irr_clust_dist, method = "ward.D2")
# 
# if(dimension == "conflict"){
#   c_optim_num_irr <- 3
# }else{
#   optcl_irr <- maptree::kgs(cluster = irr_clust_dist_hc, diss = irr_clust_dist, maxclust = 40)
#   
#   c_optim_num_irr <- as.numeric(names(optcl_irr[which(optcl_irr == min(optcl_irr))]))
# }
# 
# 
# cat(">>> Optim number of irregular cluster is:", c_optim_num_irr ,"/n")
# 
# clusters_irr <- cutree(irr_clust_dist_hc, k = c_optim_num_irr)
# 
# eco_grid_sf_irr <- motif::lsp_add_clusters(irr_clust_lsp, clusters_irr, window = shp["id"])
# 
# metrics_clust_irr <- motif::lsp_add_quality(eco_grid_sf_irr, irr_clust_dist, type = "cluster") %>%
#   dplyr::group_by(clust) %>%
#   dplyr::summarise(inhomogeneity = mean(inhomogeneity),
#                    distinction = mean(distinction),
#                    quality = mean(quality)) %>% 
#   sf::st_drop_geometry() %>% 
#   dplyr::mutate(cluster_type = "irregular")
# 
###================================================##
###============== SAVER RESULTS ===================##
###================================================##

cat(">>>Calculating summary metrics for each cluser \n")
clust_mtrs <- list()

clust_mtrs$reg_clust_values <- get_sum_cl_mtrs(rast_paths = file_paths %>% filter(type == dimension) %>% pull(path), eco_grid_sf = eco_grid_sf, world_mask = world_mask, shp_ext = extent(shp)) %>% 
  drop_na()


lapply(unique(clust_mtrs$reg_clust_values$clust), function(i){
  
  clust_mtrs$reg_clust_values %>%
    dplyr::filter(clust == i) %>% 
    dplyr::select(-clust) %>% 
    psych::describe(., na.rm = T) %>% 
    dplyr::select(-vars, -n) %>% 
    dplyr::mutate(clust = i) %>% 
    as_tibble(., rownames = "variable")
}) %>% 
  dplyr::bind_rows() %>% 
  write_csv(., paste0(dest_dir, dimension, "_reg_cluster_statistics.csv"))


clust_mtrs$reg_rel_change <- clust_descriptives(clust_sum = clust_mtrs$reg_clust_values)


if(length(unique(clust_mtrs$reg_clust_values$clust)) <= 7){
  
  cluster_text <- labeling_function(db = clust_mtrs$reg_clust_values, n_vars = 6)
  write_csv(cluster_text, paste0(dest_dir, dimension, "_reg_cluster_text_description.csv"))
  
  
}else{
  cat("Number of cluster greather than 7")
}


# clust_mtrs$irr_clust_values <- get_sum_cl_mtrs(rast_paths = file_paths %>% filter(type == dimension) %>% pull(path), eco_grid_sf = eco_grid_sf_irr, world_mask = world_mask, shp_ext = extent(shp))

# clust_mtrs$irr_rel_change <- clust_descriptives(clust_sum = clust_mtrs$irr_clust_values)



#writexl::write_xlsx(clust_mtrs[c("irr_rel_change", "reg_rel_change")], paste0(dest_dir, dimension, "_cluster_summary_metrics.xlsx"))
write_csv(clust_mtrs$reg_rel_change, paste0(dest_dir, dimension, "_cluster_rel_change.csv"))
write_csv(clust_mtrs$reg_clust_values, paste0(dest_dir, dimension, "_reg_cluster_values_extracted.csv"))

#write_csv(clust_mtrs$irr_clust_values, paste0(dest_dir, dimension, "_irr_cluster_values_extracted.csv"))



# cat(">>>saving clusters motif metrics to dest_dir \n")
# clust_mtrs_cl <- bind_rows(metrics_clust,  metrics_clust_irr)
# write_csv(clust_mtrs_cl, paste0(dest_dir, dimension, "_cluster_eval.csv"))
# 
cat(">>>writing clusters grid to dest_dir \n")
sf::st_write(eco_grid_sf, paste0(dest_dir, dimension,"_regular_clust.shp"), delete_dsn = T)
#sf::st_write(eco_grid_sf_irr, paste0(dest_dir, dimension,"_irregular_clust.shp"), delete_dsn = T)


cat(">>>Making boxplots \n")

g1 <- clust_mtrs$reg_clust_values %>% 
  pivot_longer(., cols =  -clust, names_to = "variable", values_to = "values") %>% 
  dplyr::mutate(clust = paste0("cluster_",clust)
                #,values = ifelse(values < 1, 1, values)
  ) %>% 
  #filter(values > 1) %>%
  ggplot(aes(y= values, x = clust))+
  geom_violin(trim=FALSE, fill="gray")+
  geom_boxplot(width=0.1, fill="white")+
  facet_wrap(~variable, scales = "free")+
  xlab("")+
  theme_bw(base_size = 10)

if(dimension != "climate"){
  g1 <- g1 +
    scale_y_log10()
 
}


ggsave(plot = g1, filename = paste0(dest_dir, dimension,'_regular_boxplots.png'), dpi = 400, width = 15, height = 8, units = "in")

ggsave(plot = g2, filename = paste0(dest_dir, dimension,'_irregular_boxplots.png'), dpi = 400, width = 15, height = 8, units = "in")



#######################################################
###### NEW CONFLICT CLUSTERING #######################
#####################################################

pop_dens <- raster::raster(paste0(root,"/data/",iso, "/population_density/medn_popd.tif"))

shp <- raster::shapefile(paste0(root,"/data/", iso, "/_shps/",iso,".shp" )) %>% 
  sf::st_as_sf() %>% 
  dplyr::mutate(id = 1:nrow(.))


world_mask <- raster::raster(paste0(root,"/data/_global/masks/mask_world_1km.tif")) %>% 
  raster::crop(., extent(shp))

knl <- raster::raster(paste0(root, "/data/", iso, "/conflict/conflict_kernel_density.tif"))
crs(knl) <- crs(world_mask)


conflict_raw <-  get_conflic_data(root = root,
                                  iso = iso,
                                  country = "Kenya") %>% 
  purrr::pluck(1)




coordinates(conflict_raw) <- ~x+y
crs(conflict_raw) <- crs(world_mask)
conflict_sf <- sf::st_as_sf(conflict_raw)

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
                FATALITIES = ifelse(pop_dens == 0 , 1, FATALITIES/pop_dens),
                knl = raster::extract(knl, st_as_sf(.) %>% st_cast("POINT")%>% st_coordinates())) %>%
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

cords <- rbind(x$ind$coord) %>% 
  as_tibble() %>% 
  dplyr::mutate(dst = sqrt(( Dim.1 - 0)^2 + (Dim.2 - 0 )^2   ),
                rng = dst > quantile(dst)[3] + IQR(dst)*1.5) 

thr <- cords %>% dplyr::filter(rng) %>% pull(dst) %>% quantile %>% .[4]

wh <- cords %>% 
  dplyr::mutate(wh = ifelse(dst > thr, 7, 1)) %>% 
  dplyr::pull(wh)


x_wh <- to_cluster %>% 
  dplyr::mutate(knl = ifelse(is.na(knl), 0, knl)) %>% 
  dplyr::select(-ov) %>% 
  st_drop_geometry() %>%  
  FactoMineR::PCA(., scale.unit = T, ncp =2, graph = F , row.w = wh)


x_wh <- rbind(x_wh$ind$coord)%>% 
  dist(., method = "euclidean")

km <- stats::kmeans(x_wh, 3)
table(km$cluster)


original_df <- conflict_sf %>% 
  dplyr::mutate(ov = unlist(st_intersects(conflict_sf, grd)),
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
  dplyr::ungroup() %>% 
  dplyr::bind_cols(clust = as.character(km$cluster))


cluster_labels <- original_df %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(EVENTS:FATALITIES, contains("clust")) %>% 
  get_cluster_labels()

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
  dplyr::select(-x, -clust_km)

g <- make_cluster_plots(df = to_boxplot)
x11();g

ggsave(g, filename= paste0(root,"/data/", iso, "/_results/cluster_results/conflict/conflict_regular_clust_boxplots.png"),
       dpi=300, 
       height=8,
       width=15,
       units="in")


get_cluster_statistics(df = to_boxplot) %>% 
  write.csv(., paste0(root, "/data/", iso, "/_results/cluster_results/conflict/conflict_regular_clust_rel_change.csv"), row.names = F)
