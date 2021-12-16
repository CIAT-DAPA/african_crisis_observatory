# Spatial patterns clustering
# African Crisis Observatory
# Alliance Bioversity-CIAT, 2021

options(warn = -1, scipen = 999)
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
  
  global_mtrs <- clust_sum %>% 
    dplyr::select(-clust) %>% 
    apply(., 2, function(i){
      x <- i[i>0]
      ret <- median(x, na.rm = T)
      return(ret)
    })
  
  rast_mts <- clust_sum %>% 
    dplyr::group_by(clust) %>% 
    dplyr::summarise(across(everything() , function(i){
      x <- i[i > 0 ]
    ret <- median(x, na.rm = T)
    return(ret)
    }  )) %>%
    ungroup() %>% 
    dplyr::select(- clust) %>% 
    t %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    as_tibble()
  
  names(rast_mts)<- c("Variables", paste0("clust_", 1:(ncol(rast_mts)-1) )) 
  
  ret <-  bind_cols(
    rast_mts,
    apply(rast_mts, 1, function(i){
      rast_nm <- i[1]
      fnl_nm  <- paste0( "_rel_change")
      glb_med <- global_mtrs[names(global_mtrs) == rast_nm]
      vals <- round((as.numeric(i[-1]) - glb_med)/glb_med*100, 2)
      names(vals) <- paste0(names(i)[-1],"_rel_change")
      return(vals)
    }) %>% t %>% 
      as_tibble(),
    global_metrics = global_mtrs
  ) 
  
  return(ret)
}

root <- 'D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/'
country_iso2 <- "SDN"
country <- 'Sudan'

baseDir <- paste0(root, "data/",country_iso2)


df <- readxl::read_excel(paste0(root,'/Country_pathways.xlsx'), sheet = 2)
df <- df %>% dplyr::filter(Country == country & Dimension == 'Climate')
clm <- df$Variable
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


stopifnot("File not found in paths. " = all(check_files))

### load and reclassify all rasters

for(dimension in unique(file_paths$type)){
  
  cat(">>> starting process for: ", dimension, "\n")
  
  
  dest_dir <- paste0(baseDir, "/_results/cluster_results/", dimension, "/")
  
  shp <- raster::shapefile(paste0(baseDir,"/_shps/",country_iso2,".shp" )) %>% 
    sf::st_as_sf() %>% 
    dplyr::mutate(id = 1:nrow(.))
  
  grd <- st_make_grid(st_bbox(extent(shp)+2), cellsize = 0.2, square =  T) %>% 
    st_as_sf(.) %>%
    dplyr::mutate(id = 1:nrow(.))
  
  
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
  
  # =========================================================== #
  # Irregular clusterin patters
  # =========================================================== #
  
  cat(">>> Calculating irregular clusters /n")
  
  irr_clust_lsp <- lsp_signature(eco_data, type = "incove",
                                 window = shp["id"], normalization = "pdf") 
  
  irr_clust_dist <- lsp_to_dist(irr_clust_lsp, dist_fun = "jensen-shannon")
  irr_clust_dist_hc <-  hclust(irr_clust_dist, method = "ward.D2")
  
  if(dimension == "conflict"){
    c_optim_num_irr <- 3
  }else{
    optcl_irr <- maptree::kgs(cluster = irr_clust_dist_hc, diss = irr_clust_dist, maxclust = 40)
    
    c_optim_num_irr <- as.numeric(names(optcl_irr[which(optcl_irr == min(optcl_irr))]))
  }
  
  
  cat(">>> Optim number of irregular cluster is:", c_optim_num_irr ,"/n")
  
  clusters_irr <- cutree(irr_clust_dist_hc, k = c_optim_num_irr)
  
  eco_grid_sf_irr <- motif::lsp_add_clusters(irr_clust_lsp, clusters_irr, window = shp["id"])
  
  metrics_clust_irr <- motif::lsp_add_quality(eco_grid_sf_irr, irr_clust_dist, type = "cluster") %>%
    dplyr::group_by(clust) %>%
    dplyr::summarise(inhomogeneity = mean(inhomogeneity),
                     distinction = mean(distinction),
                     quality = mean(quality)) %>% 
    sf::st_drop_geometry() %>% 
    dplyr::mutate(cluster_type = "irregular")
  
  ###================================================##
  ###============== SAVER RESULTS ===================##
  ###================================================##
  
  cat(">>>Calculating summary metrics for each cluser \n")
  clust_mtrs <- list()
  
  clust_mtrs$reg_clust_values <- get_sum_cl_mtrs(rast_paths = file_paths %>% filter(type == dimension) %>% pull(path), eco_grid_sf = eco_grid_sf, world_mask = world_mask, shp_ext = extent(shp)) %>% 
    drop_na()
  
  
  clust_mtrs$reg_rel_change <- clust_descriptives(clust_sum = clust_mtrs$reg_clust_values)
  
  
  clust_mtrs$irr_clust_values <- get_sum_cl_mtrs(rast_paths = file_paths %>% filter(type == dimension) %>% pull(path), eco_grid_sf = eco_grid_sf_irr, world_mask = world_mask, shp_ext = extent(shp))
  
  clust_mtrs$irr_rel_change <- clust_descriptives(clust_sum = clust_mtrs$irr_clust_values)
  
 
  
  writexl::write_xlsx(clust_mtrs[c("irr_rel_change", "reg_rel_change")], paste0(dest_dir, dimension, "_cluster_summary_metrics.xlsx"))
  write_csv(clust_mtrs$reg_clust_values, paste0(dest_dir, dimension, "_reg_cluster_values_extracted.csv"))
  write_csv(clust_mtrs$irr_clust_values, paste0(dest_dir, dimension, "_irr_cluster_values_extracted.csv"))
  
  
  
  cat(">>>saving clusters motif metrics to dest_dir \n")
  clust_mtrs_cl <- bind_rows(metrics_clust,  metrics_clust_irr)
  write_csv(clust_mtrs_cl, paste0(dest_dir, dimension, "_cluster_eval.csv"))
  
  cat(">>>writing clusters grid to dest_dir \n")
  sf::st_write(eco_grid_sf, paste0(dest_dir, dimension,"_regular_clust.shp"), delete_dsn = T)
  sf::st_write(eco_grid_sf_irr, paste0(dest_dir, dimension,"_irregular_clust.shp"), delete_dsn = T)
  
  
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
  
  g2 <- clust_mtrs$irr_clust_values %>% 
    pivot_longer(., cols =  -clust, names_to = "variable", values_to = "values") %>% 
    dplyr::mutate(clust = paste0("cluster_",clust)) %>% 
    ggplot(aes(y= values, x = clust))+
    geom_violin(trim=FALSE, fill="gray")+
    geom_boxplot(width=0.1, fill="white")+
    facet_wrap(~variable, scales = "free")+
    xlab("")+
    theme_bw(base_size = 10)
  
  if(dimension != "climate"){
    g1 <- g1 +
      scale_y_log10()
    g2 <- g2 +
      scale_y_log10()
  }
  
  
  ggsave(plot = g1, filename = paste0(dest_dir, dimension,'_regular_boxplots.png'), dpi = 400, width = 15, height = 8, units = "in")
  
  ggsave(plot = g2, filename = paste0(dest_dir, dimension,'_irregular_boxplots.png'), dpi = 400, width = 15, height = 8, units = "in")
  
  
}#END FOR



#####PUNTO 1
#se define la funcion para la camintata aleatoria con sus respectivos parametros
caminata_aleatoria <- function(n_pasos, # numero de pasos para la caminata aleatoria
                               probabilidad, # vector de dos valores c(arriba, abajo, izquierda, derecha) con las probabilidades para el eje x y eje y
                               origen, # puntos de origen caminta aleatoria z^2
                               graficar = TRUE  #si se debe generar o no la grafica TRUE o FALSE
                               ){
  
  #set.seed(1234) # se puede establecer una semilla para que los resultados den identicos en cada iteración
  
  #se crea una matrix con los posibles estados de la caminata aleatoria
  estados <-  c(sample(c(-1,1), n_pasos/2, replace = T, prob = c(probabilidad[1],probabilidad[2]) ), 
             sample(c(-1,1), n_pasos/2, replace = T, prob = c(probabilidad[3],  probabilidad[4])))
  
  #se crea una matrix vaica para almacenar los resultados de la caminata aleatoria
  caminata <- matrix(0, ncol = 2, nrow = n_pasos)
  
  
  #se adicionan las direcciones de la caminata dados los posibles estados 
  indx <- cbind(seq(n_pasos), sample(c(1, 2), n_pasos, TRUE))
  caminata[indx] <- estados
  # se calcula la camianta aleatoria como la suma acumulada de los valores de la matrix de estados 
  
  caminata[,1] <- cumsum(caminata[, 1])
  caminata[, 2] <- cumsum(caminata[, 2])
  #se adiciona el punto de origen a la caminata
  caminata <- rbind(origen, caminata)
  
 
  #se convierte la matrix a un dataframe
  caminata <- as.data.frame(caminata)
  
  #se define la tabla a retornar por la  function
  retornar <- data.frame(puntos_evolucion_eje_x  = caminata[,1],
             puntos_evolucion_eje_y = caminata[,2],
             sorteo_aleatorio_eje = c(0, estados) )
  #mediante un if se controla cuando plotear la grafica de la caminata
  if(graficar){
    caminata[ ,3]<- "blue"
    caminata[1,3] <- "red"
    caminata[nrow(caminata), 3] <- 'red'
    
    
    plot(caminata[,c(1,2)], xlab = "eje x", ylab = "eje y", main = "caminata aleatoria", col = caminata$V3, type = "b")
    lines(caminata[,c(1,2)])
  }
  
  return(retornar)
  
}

#se esejecuta la funcion con los parametros establecidos y se almacenan los reaultados en un objeto
tabla_resultado <- caminata_aleatoria(n_pasos = 50,
                   probabilidad = c(0.2,0.3, 0.2, 0.3),
                   origen =c(0,0))



###PUNTO 2


#funcion para aproximar la probabilidad de regresar la origen
simulacion<- function(n_sim, # numero de simulaciones
                      n_pasos,# numero de pasos para la simulacion,
                      probabilidad, # vector de dos valores c(p, 1-p) con las probabilidades
                      origen # puntos de origen caminta aleatoria z^2
                      ){
  #se define una vector vacio para alacenar las caminatas aleatorias atraves de las iteraciones
  resultados <- c()
  #se crea un for loop para iterar hasta el numero de similaciones n_sim
  for(i in 1:n_sim){
    #se ejecuta la funcion de la caminata aleatoria y se almacenan sus resultados
    tabla <- caminata_aleatoria(n_pasos = n_pasos,
                       probabilidad =probabilidad,
                       origen =origen,
                       graficar = F)
    #se identifica cuando la caminata regreso al punto de origen
    #nota: la primera fila no se tendra en cuenta pues esta representa el punto de origen
    #se usa la funcion apply para evaluar en cada fila si el punto de la caminata es igual al punto de origen
    #en cuyo caso el resultado sera TRUE, lo que indicaria que la caminata retorno al punto de origen
    regreso_origen <- apply(tabla[,c(1,2)] == origen, 1, all)
    #se almacena en el vector resultado un valor de TRUE si la caminata almenos una ves regreso al punto de origen

  resultados[i] <- any(regreso_origen[-1])
  }
  
  return(cat(paste0("La probabilidad de regreso a la posicion de origen es: ", sum(resultados)/n_sim )))
  
}

#se esejecuta la funcion con los parametros establecidos
simulacion(n_sim= 1000, # numero de simulaciones
           n_pasos= 100,# numero de pasos para la simulacion,
           probabilidad = c(0.3,0.2,0.3,0.2), # vector de dos valores c(p, 1-p) con las probabilidades
           origen= c(0.0))

