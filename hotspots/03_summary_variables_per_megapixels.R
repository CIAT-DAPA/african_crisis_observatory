# Produce summarized metrics per mega-pixel - master code
# African Crisis Observatory
# By: H. Achicanoy
# Alliance Bioversity-CIAT, 2021

# R options
g <- gc(reset = T); rm(list = ls()) # Emptying the garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, terra, exactextractr, vroom))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

iso <- 'KEN'
# Country shapefile
shp <- raster::shapefile(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
# Conflict clusters
clt <- raster::shapefile(paste0(root,'/data/',iso,'/_results/cluster_results/conflict/conflict_regular_clust.shp'))
if(sum(is.na(clt@data$clust)) > 0){
  clt <- clt[!is.na(clt@data$clust),]
}
# Raster template
tmp <- raster::raster(paste0(root,'/data/_global/masks/mask_world_1km.tif'))
tmp <- tmp %>% raster::crop(raster::extent(shp)) %>% raster::mask(shp)

vrs <- list.files(path = paste0(root,'/data/',iso), pattern = '*.tif$', recursive = T, full.names = T)
vrs <- vrs[-grep(pattern = '_results', x = vrs)]
vrs <- vrs[-grep(pattern = '_mask', x = vrs)]

rst <- vrs %>% purrr::map(.f = function(vr){
  r <- raster::raster(x = vr)
  r <- r %>% terra::resample(tmp)
  return(r)
}) %>% raster::stack()
vls <- exactextractr::exact_extract(x = rst, y = clt, full_colnames = T)

dfm <- 1:length(vls) %>% purrr::map(.f = function(i){
    df <- vls[[i]]
    df <- cbind(clt@data[i,], df)
    return(df)
  }) %>% dplyr::bind_rows()

names(dfm)[1:9] <- c('id','EVENTS',
                     'TYPE_RICHNESS','SUBTYPE_RICHNESS',
                     'ACTOR1_RICHNESS','ACTOR2_RICHNESS',
                     'FATALITIES','KERNEL_DENSITY','CONFLICT_CLUSTER')
dfm <- dfm %>%
  dplyr::select(id:CONFLICT_CLUSTER, acess:trnd_wasting,
                deforest:trnd_male_edu, irrigation:lvst_tlu,
                cvar_migration:trnd_popd, cvar_piped_water:soil_pH,
                KEN_AWE,KEN_rwi,
                cvar_cwdf:THI_3_trnd, cvar_aet:flood, cvar_prec:trnd_prec, cvar_tmax:trnd_tmax)

pca <- dfm %>%
  dplyr::select(-id,-CONFLICT_CLUSTER) %>%
  FactoMineR::PCA(X = ., scale.unit = T, ncp = ncol(.), graph = F)

corrplot::corrplot(corr = pca$var$cos2[,1:20] %>%
                     round(digits = 2) %>%
                     base::as.data.frame() %>%
                     tidyr::drop_na() %>%
                     base::as.matrix(), is.corr = FALSE)

# dfm %>%
#   ggplot2::ggplot(aes(x = KEN_AWE, y = medn_migration)) +
#   ggplot2::geom_point() +
#   ggplot2::theme_bw() +
#   ggplot2::facet_wrap(~CONFLICT_CLUSTER)

vroom::vroom_write(x = dfm, file = paste0('D:/',iso,'_stats.csv'), delim = ',')
dfm <- vroom::vroom(paste0('D:/',iso,'_stats.csv'), delim = ',')
