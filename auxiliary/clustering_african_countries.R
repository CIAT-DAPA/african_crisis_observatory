######################################################################
### Script to generate conflict cluster for Africa continent
#' @author: Andres Mendez and Benson Kenduiywo
#' year: 2023
#' 

#'importar paquetes
rm(list=ls(all=TRUE))
require(pacman)
pacman::p_load(terra, sf, exactextractr, tidyverse, FactoMineR)

region <- 'East_and_Horn_of_Africa_zscored'
#' define global vars and load inputs
root <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/"
out <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/UNCHR/"

wrld_shp <- terra::vect(paste0(root, "/_global/world_shapefile/all_country/all_countries.shp"))

af_shp <- wrld_shp[wrld_shp$CONTINENT == "Africa"] 

countries_lst <- c('Burundi','Djibouti', 'Eritrea', 'Ethiopia', 'Kenya','Rwanda','Sudan','Somalia', 'South Sudan', 'Tanzania', 'Uganda')
countries_iso <- c("BDI", "DJI", "ERI", "ETH", "KEN", "RWA", "SDN", "SOM", "SSD", "TZA", "UGA" )

#' consider downloading lowest level because the world shapefile is missing some countries like SSD
conflic_raw <- readxl::read_excel(paste0(root, "/_global/conflict/Africa_1997-2022_Jul08.xlsx"))
cols <- c("EVENT_DATE", "YEAR", "EVENT_TYPE","SUB_EVENT_TYPE", "ACTOR1", "ACTOR2","FATALITIES", "LATITUDE", "LONGITUDE", "COUNTRY" ) 

#load population density rasters 
fls <- list.files(paste0(root, "/", countries_iso, "/", "population_density/" ), 
                  pattern = "medn_popd.tif$",
                  full.names = T)

rsts <- lapply(fls, terra::rast)

popd_af <- purrr::reduce(rsts, terra::mosaic, fun = "mean")

terra::writeRaster(popd_af, paste0(out,"pop_af.tif"))

#load and process conflict data (uptadted to JULY 2022)
conflict_clean <- conflic_raw[conflic_raw$COUNTRY %in% countries_lst, cols ] 
#' define starting year for time series
mn <- by(conflict_clean$YEAR, conflict_clean$COUNTRY, function(i){to_ret <- list(); to_ret$min = min(i, na.rm = T); to_ret$max = max(i, na.rm =T);return(to_ret)}, simplify = T)
start_year <- max(sapply(mn, function(i){i$min}))
end_yead <- max(sapply(mn, function(i){i$max}))
rm(mn)

conflict_clean <- conflict_clean[conflict_clean$YEAR >=  start_year, ]
conflict_sf <- sf::st_as_sf(conflict_clean, coords = c("LONGITUDE", "LATITUDE"))

#' generate grid for africa continent
grd <- sf::st_make_grid(st_bbox(terra::ext(af_shp)+2), cellsize = 0.2, square =  T) %>% 
  st_as_sf(.) %>%
  dplyr::mutate(id = 1:nrow(.))
#' filter megapixels where conflict was observed
grd<- grd[!is.na(sapply(st_intersects(grd, conflict_sf) , function(i){if(length(i)>0){return(1)}else{return(NA)}})), ]

#' Add to conflict sf the grid id
conflict_sf$grd_id <- sapply(st_intersects(conflict_sf, grd), function(i){ifelse(length(i)==0,NA, unlist(grd[i, "id"]))})
conflict_sf$pop_dens <- terra::extract(popd_af, conflict_sf)$layer
#remove ACLED location with pop density NA
conflict_sf <- conflict_sf[!is.na(conflict_sf$pop_dens), ]
#' min Max normalization
minMax <- function(x){
  n <- x-min(x, na.rm = T)
  d <- max(x, na.rm=T)-min(x, na.rm=T)
  return(n/d)
}

#' Z normalization

zscore <- function(y, log=FALSE){
  # +.05 to avoid NA from log(x) where x <= 0
  if (log) y <- log(y+.05)
  return((y - mean(y, na.rm=TRUE) ) / (sd(y, na.rm=TRUE)))
}

#' by grid id, aggregate variables  and joint in a data.frame
aa <- data.frame(id = sort(unique(conflict_sf$grd_id)),
                 EVENTS =  (as.vector(by(conflict_sf$EVENT_TYPE, conflict_sf$grd_id, length))),#median(EVENTS, na.rm = T),
                 TYPE_RICHNESS = (as.vector(by(conflict_sf$EVENT_TYPE, conflict_sf$grd_id, function(i){length(unique(i))}))),#max(TYPE_RICHNESS, na.rm = T),
                 SUBTYPE_RICHNESS = (as.vector(by(conflict_sf$SUB_EVENT_TYPE, conflict_sf$grd_id, function(i){length(unique(i))}))),#max(SUBTYPE_RICHNESS, na.rm = T),
                 ACTOR1_RICHNESS = (as.vector(by(conflict_sf$ACTOR1, conflict_sf$grd_id, function(i){length(unique(na.omit(i)))}))),#max(ACTOR1_RICHNESS, na.rm = T),
                 ACTOR2_RICHNESS = (as.vector(by(conflict_sf$ACTOR2, conflict_sf$grd_id, function(i){length(unique(na.omit(i)))}))),#,
                 FATALITIES = (as.vector(by(conflict_sf$FATALITIES, conflict_sf$grd_id, sum))),
                 pop_dens = (as.vector(by(conflict_sf$pop_dens, conflict_sf$grd_id, mean, na.rm = T)))
)

aa$EVENTS <- aa$EVENTS/aa$pop_dens
aa$FATALITIES <- aa$FATALITIES/aa$pop_dens
aa <- aa[!is.nan(aa$FATALITIES) | !is.infinite(aa$EVENTS) , ]
aa <- aa[!is.infinite(aa$FATALITIES),]
aa$pop_dens<- NULL



bb <- data.frame(id = sort(unique(conflict_sf$grd_id)),
                 EVENTS =  zscore(as.vector(by(conflict_sf$EVENT_TYPE, conflict_sf$grd_id, length))),#median(EVENTS, na.rm = T),
                 TYPE_RICHNESS = zscore(as.vector(by(conflict_sf$EVENT_TYPE, conflict_sf$grd_id, function(i){length(unique(i))}))),#max(TYPE_RICHNESS, na.rm = T),
                 SUBTYPE_RICHNESS = zscore(as.vector(by(conflict_sf$SUB_EVENT_TYPE, conflict_sf$grd_id, function(i){length(unique(i))}))),#max(SUBTYPE_RICHNESS, na.rm = T),
                 ACTOR1_RICHNESS = zscore(as.vector(by(conflict_sf$ACTOR1, conflict_sf$grd_id, function(i){length(unique(na.omit(i)))}))),#max(ACTOR1_RICHNESS, na.rm = T),
                 ACTOR2_RICHNESS = zscore(as.vector(by(conflict_sf$ACTOR2, conflict_sf$grd_id, function(i){length(unique(na.omit(i)))}))),#,
                 FATALITIES = zscore(as.vector(by(conflict_sf$FATALITIES, conflict_sf$grd_id, sum)))
)



#' PCA 


pca1 <- FactoMineR::PCA(aa[, -1], scale.unit = T, ncp = 2)

plot(pca1)


cords <- rbind(pca1$ind$coord) %>% 
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
pca_w <-  FactoMineR::PCA(aa[, -1], scale.unit = T, ncp =2, graph = T , row.w = wh)

plot(pca_w, choix = "varcor")
plot(pca_w)
x_wh <- rbind(pca_w$ind$coord)%>% 
  dist(., method = "euclidean")

set.seed(1000) #https://towardsdatascience.com/three-versions-of-k-means-cf939b65f4ea
km_clust <- kmeans(x_wh, centers = 3, iter.max = 1000, nstart = 20, algorithm="MacQueen")#stats::kmeans(x_wh, 3)
table(km_clust$cluster)


#km_clust <- kmeans(dist(pca1$ind$coord[, c(1,2)], method = "euclidean"), centers = 3)
#table(km_clust$cluster)

plot(pca1$ind$coord[,1], pca1$ind$coord[,2], col = km_clust$cluster)

aa$km_cluster <- km_clust$cluster

aa %>% 
  dplyr::select(EVENTS:FATALITIES, starts_with("km_")) %>% 
  dplyr::mutate(km_cluster = as.factor(km_cluster)) %>% 
  tidyr::pivot_longer(-km_cluster, names_to = "var", values_to = "vals") %>% 
  ggplot(aes(y= vals, x = km_cluster))+
  geom_violin(trim=FALSE, fill="gray")+
  geom_boxplot(width=0.1, fill="white")+
  facet_wrap(~var, scales = "free")+
  xlab("")+
  theme_bw(base_size = 10)+ 
  scale_y_log10()

#cluster labeling

conf_cluts_vals <- aa %>% 
  dplyr::select(EVENTS:FATALITIES, starts_with("km_cluster")) %>% 
  group_by(km_cluster) %>% 
  dplyr::summarise(m_events = mean(EVENTS)) %>% 
  arrange(desc(m_events)) %>% 
  dplyr::mutate(label = c("High conflict", "Moderate conflict", "Limited conflict"),
                short_label = c("High", "Moderate", "Limited")) %>% 
  dplyr::mutate(across(everything(.), as.character),
                label = factor(label, levels = c("High conflict","Moderate conflict",  "Limited conflict")))

aa <- base::merge(aa, conf_cluts_vals , by = "km_cluster", all.x = T ) 

#' merge grid sf with conflict megapixel aggregated variables to consolidate the to_cluster data.frame
grd_conflict <- base::merge(grd, aa , by = "id", all.x = T ) 

#save results
#sf::st_write(grd_conflict, paste0(root, "/conflict_clusters_af.geojson"))
writeVector(vect(grd_conflict), filename = paste0(out,region,'.shp'))










