######################################################################
### Script to generate conflict cluster for Africa continent
#' @author: Andres Mendez and Benson Kenduiywo
#' year: 2023
#' #NOTEs K-Means, Hierchical clustering and on PCA https://rpubs.com/williamsurles/310847

#'importar paquetes
rm(list=ls(all=TRUE))
require(pacman)
pacman::p_load(terra, sf, exactextractr, tidyverse, FactoMineR, fpc, dbscan)

#' ===================================================================================
#' 1). define global variables 
#' ===================================================================================
root <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/"
out <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/UNCHR/"

#' Define regions and corresponding ISO3 codes
MENA <- c("Iraq", "Israel", "Jordan", "Kuwait",	"Lebanon", "Saudi Arabia", "Syria",
          "United Arab Emirates",	"Yemen", "Algeria", "Egypt", "Libya", "Mauritania",
          "Morocco", "Tunisia",	"Western Sahara",	"Bahrain", "Oman", "Palestinian Territory",	"Qatar")
menaISO <- c("IRQ", "ISR", "JOR", "KWT", "LBN", "SAU", "SYR", "ARE", "YEM",	"DZA", 
             "EGY",	"LBY", "MRT",	"MAR", "TUN", "ESH", "BHR",	"OMN", "PSE",	"QAT")

EHorn <- c('Burundi','Djibouti', 'Eritrea', 'Ethiopia', 'Kenya','Rwanda','Sudan','Somalia', 'South Sudan', 'Tanzania', 'Uganda')
eastAsiaPacific <- c("KAZ", "KGZ", "TJK", "TKM", "UZB", "AUS", "CHN",	"HKG", "MAC",	"JPN",	
                     "KOR",	"NZL", "PNG",	"TWN", "IND",	"NPL", "LKA",	"BGD", "KHM", "IDN",
                     "LAO", "MYS", "MMR",	"PHL", "THA",	"VNM", "AFG",	"IRN", "PAK",	"ASM",	
                     "BTN",	"BRN", "CXR",	"CCK", "COK",	"FJI", "GUM",	"KIR", "PRK",  "MDV",
                     "MHL",	"FSM", "MNG", "NRU", "NCL", "NIU", "MNP", "PLW", "WSM", "SGP",
                     "SLB",	"TLS", "TKL",	"TON", "TUV",	"UMI", "VUT",	"WLF")

Americas_ISO <- c("BLZ", "BRA", "CHL", "COL",	"CRI", "CUB", "ECU", "SLV",	"GTM", "HND",	"MEX", "PAN",	"PRY", "PER",	"URY", "VEN", "ABW", "CAN",	"DOM", "GUY",	"HTI", "TTO", "USA", "AIA",	"ATG",	"BHS",	"BRB",	"BMU",	"BOL",	"VGB",	"CYM",	"DMA", "FLK",	"PYF", "GRD",	"GLP",	"JAM",	"MTQ",	"NIC",	"PRI",	"KNA",	"LCA",	"SPM",	"VCT",	"SGS",	"SUR",	"TCA",	"VIR")
Americas_lst <- c("Argentina", "Belize", "Brazil",	"Chile",	"Colombia",	"Costa Rica",	"Cuba",	"Ecuador",	"El Salvador", "Guatemala",	"Honduras",	"Mexico",	"Panama",	"Paraguay",	"Peru",	"Uruguay",	"Venezuela",	"Aruba",	"Canada",	"Dominican Republic",	"Guyana",	"Haiti",	"Trinidad and Tobago",	"United States",	"Anguilla",	"Antigua and Barbuda",	"Bahamas",	"Barbados",	"Bermuda",	"Bolivia",	"British Virgin Islands",	"Cayman Islands",	"Dominica",	"Falkland Islands",	"French Polynesia",	"Grenada",	"Guadeloupe",	"Jamaica",	"Martinique",	"Nicaragua",	"Puerto Rico",	"Saint Kitts and Nevis",	"Saint Lucia",	"Saint Pierre and Miquelon",	"Saint Vincent and Grenadines",	"South Georgia and the South Sandwich Islands",	"Suriname",	"Turks and Caicos Islands",	"Virgin Islands, US")

## DEFINE REGION< ACLED FILE AND ASIGN VAR to corresponding iso and conuntry names
region        <- 'Americas'
ACLEDFilename  <- "Americas_2023-05-09.xlsx"

countries_iso <-  Americas_ISO
countries_lst <- Americas_lst

#' ===================================================================================
#' 2. Load data
#' ===================================================================================
wrld_shp <- terra::vect(paste0(root, "/_global/world_shapefile/all_country/all_countries.shp"))

shp <- subset(wrld_shp, wrld_shp$ISO3 %in% countries_iso) #test <- wrld_shp[wrld_shp$ISO3 %in% countries_iso] 

#' consider downloading lowest level because the world shapefile is missing some countries like SSD
conflic_raw <- readxl::read_excel(paste0(root, "/_global/conflict/", ACLEDFilename))
cols <- c("EVENT_DATE", "YEAR", "EVENT_TYPE","SUB_EVENT_TYPE", "ACTOR1", "ACTOR2","FATALITIES", "LATITUDE", "LONGITUDE", "COUNTRY" ) 

#load population density rasters 
fls <- list.files(paste0(root, "/", countries_iso, "/", "population_density/" ), 
                  pattern = "medn_popd.tif$",
                  full.names = T)

rsts <- lapply(fls, terra::rast)

popd_af <- terra::rast(paste0(root, '_global/population_density/world_medn_popn.tif'))

#load and process conflict data (uptadted to JULY 2022)
conflict_clean <- conflic_raw[conflic_raw$COUNTRY %in% countries_lst, cols ] 
#' define starting year for time series
mn <- by(conflict_clean$YEAR, conflict_clean$COUNTRY, function(i){to_ret <- list(); to_ret$min = min(i, na.rm = T); to_ret$max = max(i, na.rm =T);return(to_ret)}, simplify = T)
start_year <- max(sapply(mn, function(i){i$min}))
end_yead <- max(sapply(mn, function(i){i$max}))
cat('Start year is ', start_year,'\n')
cat('End year is ', end_yead,'\n')
rm(mn)

conflict_clean <- conflict_clean[conflict_clean$YEAR >=  start_year, ]
conflict_sf <- sf::st_as_sf(conflict_clean, coords = c("LONGITUDE", "LATITUDE"))

#' ===================================================================================
#' 3. Genral Mega pixels grids
#' ===================================================================================

#' generate grid for africa continent
grd <- sf::st_make_grid(st_bbox(terra::ext(shp)+2), cellsize = 0.2, square =  T) %>% 
  st_as_sf(.) %>%
  dplyr::mutate(id = 1:nrow(.))
#' filter megapixels where conflict was observed
grd <- grd[!is.na(sapply(st_intersects(grd, conflict_sf) , function(i){if(length(i)>0){return(1)}else{return(NA)}})), ]
plot(st_geometry(grd))
#' Add to conflict sf the grid id
conflict_sf$grd_id <- sapply(st_intersects(conflict_sf, grd), function(i){ifelse(length(i)==0,NA, unlist(grd[i, "id"]))})
conflict_sf$pop_dens <- terra::extract(popd_af, conflict_sf)$layer
#remove ACLED location with pop density NA
conflict_sf <- conflict_sf[!is.na(conflict_sf$pop_dens), ]

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

aa_original <- aa
aa$EVENTS <- aa$EVENTS/aa$pop_dens
aa$FATALITIES <- aa$FATALITIES/aa$pop_dens
aa_original <- aa_original[!is.nan(aa$FATALITIES) | !is.infinite(aa$EVENTS) ,]
aa <- aa[!is.nan(aa$FATALITIES) | !is.infinite(aa$EVENTS) , ]
aa_original <- aa_original[!is.infinite(aa$FATALITIES),]
aa <- aa[!is.infinite(aa$FATALITIES),]
aa$pop_dens<- NULL

nrow(aa)
nrow(aa_original)

#' ===================================================================================
#' 4. Dimensionality reduction with unit variance of data (Z-scoer normalization)
#' ===================================================================================
#' PCA 

pca1 <- FactoMineR::PCA(aa[, -1], scale.unit = T, ncp = 2)

plot(pca1)

#' Eliminate outliers that 1.5 times IQR
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
pca_w <-  FactoMineR::PCA(aa[, -1], scale.unit = T, ncp =5, graph = T , row.w = wh)

plot(pca_w, choix = "varcor")
plot(pca_w)
x_wh <- rbind(pca_w$ind$coord[, 1:2])%>% 
  dist(., method = "euclidean")

#' ===================================================================================
#' 5. K-Means clustering
#' ===================================================================================
set.seed(1000) #https://towardsdatascience.com/three-versions-of-k-means-cf939b65f4ea
km_clust <- kmeans(x_wh, centers = 3, iter.max = 1000, nstart = 20, algorithm="MacQueen")#stats::kmeans(x_wh, 3)
table(km_clust$cluster)


#km_clust <- kmeans(dist(pca1$ind$coord[, c(1,2)], method = "euclidean"), centers = 3)
#table(km_clust$cluster)

plot(pca1$ind$coord[,1], pca1$ind$coord[,2], col = km_clust$cluster)

aa$km_cluster <- km_clust$cluster

aa_original %>% 
  dplyr::mutate(km_cluster = km_clust$cluster) %>% 
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

#' ===================================================================================
#' 6. Cluster labelling
#' ===================================================================================
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
#sf::st_write(grd_conflict, paste0(root,region, "/conflict_clusters.geojson"))

labs <- c('High', 'Moderate', 'Limited')

for( i in 1:length(labs)){
  grd_conflict$km_cluster[grd_conflict$short_label==labs[i]] <- i
}


regionShp <- wrld_shp[wrld_shp$ISO3 %in% countries_iso,]
regionShp <- sf::st_as_sf(regionShp)

library(tmap)
library(mapview)
tmap_mode("view")
map <- tm_shape(regionShp)+
  tm_borders()+
  tm_shape(grd_conflict)+
  tm_fill(col= "km_cluster", palette="-YlOrRd", title='Conflict', labels = labs) +
  tm_compass(type = "8star", position = c("right", "top")) +#c('red','yellow','orange')
  tm_scale_bar(breaks = c(0, 50, 100), text.size = 1, width=1,
               position = c("left", "bottom"))+
  tm_layout(legend.outside=F, 
            legend.text.size = 1.1,
            legend.title.size= 1.3,
            legend.frame=F,
            legend.just = c("right", "bottom"), 
            #legend.width= 1,
            #legend.height= -0.2 
  )#+ 
#tm_add_legend(type = "text", text = "labels")

#tm_format("World")
map

tmap_save(map,  dpi= 600,  height=8, width=10, units="in",
          filename=paste0(out,'/', region,'.png'))

to_save <- grd_conflict[,  c("label", "short_label", "km_cluster")]
writeVector(vect(to_save), filename = paste0(out, region,'.shp'), overwrite=TRUE)


