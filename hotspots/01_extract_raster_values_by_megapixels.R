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

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'

#primero correr la creaci'on de los cluster de conflicto

iso <- 'MLI'
# Country shapefile
shp <- terra::vect(paste0(root,'/data/', iso,'/_shps/', iso,'.shp'))
# Conflict clusters
clt <- terra::vect(paste0(root,'/data/',iso,'/_results/cluster_results/conflict/conflict_regular_clust.shp'))
if(sum(is.na(clt$clust)) > 0){
  clt <- clt[!is.na(clt$clust),]
}
# Raster template
tmp <- terra::rast(paste0(root,'/data/_global/masks/mask_world_1km.tif'))
tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)

fls <- readxl::read_excel(paste0(root,"/Hostpots_data_dictionary.xlsx")) %>% 
  dplyr::select(Variable, to_use, Code) %>% 
  dplyr::filter(as.logical(to_use)) %>% 
  dplyr::pull(Code) %>% 
  stringr::str_replace(., "\\{iso\\}", iso) %>% 
  paste0(.,".tif")

vrs <- data.frame(full_pth = list.files(path = paste0(root,'/data/',iso), pattern = '*.tif$', recursive = T, full.names = T),
                  file_names = list.files(path = paste0(root,'/data/',iso), pattern = '*.tif$', recursive = T, full.names = F)) %>%
  dplyr::filter(!str_detect(pattern = '_results', string = full_pth)) %>% 
  dplyr::filter(!str_detect(pattern = '_mask', string = full_pth)) %>% 
  dplyr::filter(!str_detect(pattern = 'old', string = full_pth)) %>%
  dplyr::mutate(file_names = stringr::str_extract(file_names, "[a-zA-Z0-9_]+.tif$")) %>% 
  dplyr::filter(file_names %in% fls) %>% 
  dplyr::pull(full_pth)


rst <- vrs %>% purrr::map(.f = function(vr){
  r <- terra::rast(x = vr)
  r <- r %>% terra::resample(tmp)
  names(r) <-  str_extract(vr, "[A-Za-z0-9_]+.tif") %>% str_replace(., ".tif", "")
  return(r)
}) %>% terra::rast()

vls <- exactextractr::exact_extract(x = rst, y = sf::st_as_sf(clt), full_colnames = T)

dfm <- 1:length(vls) %>% purrr::map(.f = function(i){
    df <- vls[[i]]
    df <- cbind(sf::st_as_sf(clt[i,]) %>% sf::st_drop_geometry(.) , df)
    return(df)
  }) %>% 
  dplyr::bind_rows() %>% 
  dplyr::select(-matches("coverage_fraction"))


names(dfm)[1:9] <- c('id','EVENTS',
                     'TYPE_RICHNESS','SUBTYPE_RICHNESS',
                     'ACTOR1_RICHNESS','ACTOR2_RICHNESS',
                     'FATALITIES','KERNEL_DENSITY','CONFLICT_CLUSTER')
# dfm <- dfm %>%
#   dplyr::select(id:CONFLICT_CLUSTER, acess:trnd_wasting,
#                 deforest:trnd_male_edu, irrigation:lvst_tlu,
#                 cvar_migration:trnd_popd, cvar_piped_water:soil_pH,
#                 KEN_AWE,KEN_rwi,
#                 cvar_cwdf:THI_3_trnd, cvar_aet:flood, cvar_prec:trnd_prec, cvar_tmax:trnd_tmax)

# pca <- dfm %>%
#   dplyr::select(-id,-CONFLICT_CLUSTER) %>%
#   FactoMineR::PCA(X = ., scale.unit = T, ncp = ncol(.), graph = F)
# 
# corrplot::corrplot(corr = pca$var$cos2[,1:20] %>%
#                      round(digits = 2) %>%
#                      base::as.data.frame() %>%
#                      tidyr::drop_na() %>%
#                      base::as.matrix(), is.corr = FALSE)

# dfm %>%
#   ggplot2::ggplot(aes(x = KEN_AWE, y = medn_migration)) +
#   ggplot2::geom_point() +
#   ggplot2::theme_bw() +
#   ggplot2::facet_wrap(~CONFLICT_CLUSTER)

to_compare <- read_csv("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data_extracted/KEN_stats.csv")


names(to_compare) == names(dfm)

vroom::vroom_write(x = dfm, file = paste0(root,'/data_extracted/',iso,'_stats.csv'), delim = ',')
#dfm <- vroom::vroom(paste0(root,'/data_extracted/',iso,'_stats.csv'), delim = ',')
