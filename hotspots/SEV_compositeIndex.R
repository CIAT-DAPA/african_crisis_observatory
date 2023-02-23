
g <- gc(reset = T); rm(list = ls())

pacman::p_load(tidyverse, terra, dplyr)


root <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO"

iso <- 'KEN'
pth <- paste0(root, '/data/',iso)
countryName <- 'Kenya'
ip <- 'ip_all'
n_vars <- 10
#SEV in the selected country
ip_var_list <- read.csv(paste0(root, "/data/", iso, "/_results/hotspots/soc_eco_all_variables.csv"))
ip_var_list$Code <- tolower(ifelse(grepl("\\{iso\\}", 
                                         ip_var_list$Code), gsub("\\{iso\\}", iso, ip_var_list$Code), 
                                   ip_var_list$Code))
#List of selected SEV
fl <- list.files(paste0(root, "/data/",iso,  "/_results/hotspots/"), pattern = paste0("_sorted_",countryName, '_', ip, ".xlsx"), full.names = T)
msk <- terra::rast(paste0(root,'/data/_global/masks/mask_world_1km.tif'))
# Country shapefile
shp <- terra::vect(paste0(pth,'/_shps/',iso,'.shp'))
stmp <- raster::shapefile(paste0(pth,'/_shps/',iso,'.shp'))

#tb      <- readxl::read_excel(fl)
#tb$Code <- tolower(tb$Code) 
#tb      <- merge(x = tb, y = ip_var_list, all.x = TRUE)

tb <- readxl::read_excel(fl) %>% 
  dplyr::mutate(Code = tolower(Code)) %>% 
  dplyr::left_join(., ip_var_list %>% 
                     dplyr::filter(IP_id == ip ) %>% 
                     dplyr::select(-Variable, - Classification), by = c("Code" = "Code")) %>% 
  dplyr::mutate(Code = ifelse(grepl("_awe", Code), paste0(iso,"_AWE"), Code),
                Code = ifelse(grepl("_rwi", Code), paste0(iso, "_rwi"), Code)) %>% 
  dplyr::slice(1:n_vars)


regions <- stringr::str_split(unique(tb$Region_value), ";") %>% 
  unlist %>%
  stringr::str_trim()

var_name <- unique(tb$Region_key)


shp_region <- stmp[stmp@data %>% pull(!!var_name) %in% regions, ] 
shp_region <- as(shp_region, 'SpatVector')

tmp  <- msk %>% 
  terra::crop(x = ., y = terra::ext(shp_region)) %>% 
  terra::mask(mask = shp_region)

# Load raster variables
htp <- 1:length(tb$Code) %>%
  purrr::map(.f = function(i){
    
    var <- tb$Code[i]
    
    cat(">>> Processing: ", var, "\n")
    # Raster template
    r <- terra::rast(list.files(path = pth, pattern = paste0(var,'.tif$'), full.names = T, recursive = T))
    
    r <- r %>%
      terra::crop(x = ., y = terra::ext(tmp)) %>%
      terra::resample(x = ., y = tmp) %>%
      terra::mask(mask = shp_region)
    names(r) <- tb$Code[i]
    return(r)
    
  })


#' list with raster files 
img <- rast(htp)
names(img)
summary(img)

##Use for other variables
minMax <- function(x){
  n <- x - minmax(x)[1]
  d <- minmax(x)[2] - minmax(x)[1]
  return(n/d)
}

#' Function for normalizing variables that have negative and positive
#' scales, like migration. We wan to capture extremes in such cases.
#'
negMaxNorm <- function(x){
  temp <- abs(x)
  return(temp/minmax(temp)[2])
}

ivars <- names(img)
#' Min max normalization


for(i in 1:length(ivars)){
  if(ivars[i] != "p90_migration"){
    img[[ivars[i]]] <- minMax(img[[ivars[i]]])
  }else{
    img[[ivars[i]]] <- negMaxNorm(img[[ivars[i]]])
  }
}

#' Change direction of variables where high values indicate high vulnearbility
vars  <- c("medn_male_edu", paste0(iso, "_AWE"), paste0(iso, "_rwi"), "medn_female_edu", "acess")
for(i in 1:length(vars)){
    img[[vars[i]]] <- 1-img[[ivars[i]]]
}

#' Compute a mean composite of a vulnerabilities

library(RColorBrewer)
c_img <- median(img, na.rm=T)

png(paste0(iso, "_Vulnerability_Composite.png"), units="px", width=2480, height=2508, res=300)
plot(c_img, col=rev(brewer.pal(nlyr(img),'RdYlGn')))
legend("bottomleft", names(img), ncol=2, bty='n')
#legend("bottomleft", legend = paste("Group", 1:3), col = 1:3, pch = 19, bty = "n")
dev.off()
