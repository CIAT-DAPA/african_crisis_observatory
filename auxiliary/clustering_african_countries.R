######################################################################
### Script to generate conflict cluster for Africa continent
#' @author: Andres Mendez and Benson Kenduiywo
#' year: 2023
#' 

#'importar paquetes
rm(list=ls(all=TRUE))
require(pacman)
pacman::p_load(terra, sf, exactextractr, tidyverse, FactoMineR)


#' define global vars and load inputs
root <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/"

wrld_shp <- terra::vect(paste0(root, "/_global/world_shapefile/all_country/all_countries.shp"))

af_shp <- wrld_shp[wrld_shp$CONTINENT == "Africa"] 

countries_lst <- c('Burundi','Djibouti', 'Eritrea', 'Ethiopia', 'Kenya','Rwanda','Sudan','Somalia', 'South Sudan', 'Tanzania', 'Uganda')
#' consider downloading lowest level because the world shapefile is missing some countries like SSD
conflic_raw <- readxl::read_excel(paste0(root, "/_global/conflict/Africa_1997-2022_Jul08.xlsx"))
cols <- c("EVENT_DATE", "YEAR", "EVENT_TYPE","SUB_EVENT_TYPE", "ACTOR1", "ACTOR2","FATALITIES", "LATITUDE", "LONGITUDE", "COUNTRY" ) 

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
           EVENTS =  zscore(as.vector(by(conflict_sf$EVENT_TYPE, conflict_sf$grd_id, length))),#median(EVENTS, na.rm = T),
           TYPE_RICHNESS = zscore(as.vector(by(conflict_sf$EVENT_TYPE, conflict_sf$grd_id, function(i){length(unique(i))}))),#max(TYPE_RICHNESS, na.rm = T),
           SUBTYPE_RICHNESS = zscore(as.vector(by(conflict_sf$SUB_EVENT_TYPE, conflict_sf$grd_id, function(i){length(unique(i))}))),#max(SUBTYPE_RICHNESS, na.rm = T),
           ACTOR1_RICHNESS = zscore(as.vector(by(conflict_sf$ACTOR1, conflict_sf$grd_id, function(i){length(unique(na.omit(i)))}))),#max(ACTOR1_RICHNESS, na.rm = T),
           ACTOR2_RICHNESS = zscore(as.vector(by(conflict_sf$ACTOR2, conflict_sf$grd_id, function(i){length(unique(na.omit(i)))}))),#,
           FATALITIES = zscore(as.vector(by(conflict_sf$FATALITIES, conflict_sf$grd_id, sum)))
)


#' merge grid sf with conflict megapixel aggregated variables to consolidate the to_cluster data.frame
grd_conflict <- base::merge(grd, aa , by = "id", all.x = T ) 

#rm(aa)


