# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Obtain Absolute Wealth Index at 5 km resolution
# Original data source: https://data.humdata.org/dataset/relative-wealth-index
# Steps:
# 1. Download manually the relative wealth index files for each country of interest:
#    [Country]_relative_wealth_index.csv
# 2. Execute this script to obtain:
#    [ISO3]_AWE.tif
# Author: Andres Mendez
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(raster, sp, sf, tidyverse, VGAM))


ICDF <- function(x, theta, gini, gdppc){
  
  alpha = (1+gini)/(2*gini)
  thr <- (1-(1/alpha))*gdppc
  
  dens_pareto <-  VGAM::qpareto(p = x , 
                                scale = thr, 
                                shape = alpha)  #(x/ (theta+y))^alpha inverse cdf  #(alpha*theta*x^(alpha - 1))/((theta+x)^(alpha+1)) - density
  
  sd_log <- sqrt(2)*VGAM::probitlink(theta = (gini+1)/2)
  mean_log <-  log(gdppc)-((sd_log)^2/2)
  
  dens_logN <- qlnorm(p = x, 
                      meanlog = mean_log, 
                      sdlog = sd_log)
  
  return(dens_pareto^(0.32)*dens_logN^(1-0.32))
  
}

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'

gdp  <- readr::read_csv("D:/african observatory/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488.csv") # Change this path properly
gini <- readr::read_csv("D:/african observatory/API_SI.POV.GINI_DS2_en_csv_v2_2445276/API_SI.POV.GINI_DS2_en_csv_v2_2445276.csv") # Change this path properly
mask <- raster::raster(paste0(root,'/data/_global/masks/mask_world_1km.tif'))

ISO3 <- "SDN"

wealth_dir <- paste0(root,'/data/_global/wealth_index')
rwi_out_dir <-  paste0(root,'/data/', ISO3,'/wealth_index')

if(!dir.exists(rwi_out_dir)){dir.create(rwi_out_dir, F, T)}

country_gini <- gini %>% 
  dplyr::filter(`Country Code` == ISO3) %>% 
  dplyr::select(tidyselect::where(~!all(is.na(.)))) %>% 
  dplyr::pull(ncol(.))
country_gini <- country_gini/100

country_gdp <- gdp %>% 
  filter(`Country Code` == ISO3) %>% 
  dplyr::pull(`2009`)

shp_country <- raster::shapefile(paste0(root,'/data/',ISO3,'/_shps/',ISO3,'.shp'))

country_rwi <- read_csv(paste0(wealth_dir,'/',ISO3,'_relative_wealth_index.csv')) %>% 
  dplyr::select(rwi, longitude, latitude)

country_rwi$rank <- rank(country_rwi$rwi, ties.method = "random")/(max(rank(country_rwi$rwi, ties.method = "random"))+1)
country_rwi$icdf <- ICDF(x     = country_rwi$rank, 
                         theta = 1, 
                         gini  = country_gini, 
                         gdppc = country_gdp) # 41986 - OCED mean

country_rwi$AWE <- country_rwi$icdf*(country_gdp/(mean(country_rwi$icdf)) )
cap <- quantile(country_rwi$AWE, probs = 0.98)
country_rwi <- country_rwi %>% 
  dplyr::mutate(AWE = ifelse(AWE > cap, cap, AWE))

coordinates(country_rwi) <-  ~longitude+latitude
crs(country_rwi) <- "+proj=longlat +datum=WGS84 +no_defs"

r <- raster(resolution = 0.04166667, ext = extent(shp_country))
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
r <- raster::mask(r, shp_country)
r_f <- rasterize(country_rwi, r, fun = mean, field = "AWE")
r_f <- raster::resample(r_f, mask %>% raster::crop(., extent(shp_country)) )

raster::writeRaster(r_f, paste0(rwi_out_dir,"/",ISO3,"_AWE.tif"), overwrite = T)

## AWE for all Africa

world_shp  <- raster::shapefile(paste0(root,'/data/_global/world_shapefile/all_country/all_countries.shp'))
wealth_dir <- paste0(root,'/data/_global/wealth_index')
mask       <- raster::raster(paste0(root,'/data/_global/masks/mask_world_1km.tif'))

wealth_iso3 <- list.files(wealth_dir) %>% substring(., 1, 3)
avalaible_iso3 <- world_shp@data$ISO3[(world_shp@data$ISO3 %in% wealth_iso3) & (world_shp@data$CONTINENT == "Africa")]

gdp  <- read_csv("D:/african observatory/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488.csv") # Change this path properly
gini <- read_csv("D:/african observatory/API_SI.POV.GINI_DS2_en_csv_v2_2445276/API_SI.POV.GINI_DS2_en_csv_v2_2445276.csv") # Change this path properly

get_gini <- function(gini, iso3){
  x <- gini%>% 
    dplyr::filter(`Country Code` == iso3) %>% 
    dplyr::select(where(~!all(is.na(.)))) %>%
    dplyr::select(ncol(.)) 
  
  return(data.frame(year = as.numeric(names(x)),
              gini = as.numeric(x)/100))
}
get_gdp  <- function(gdp, iso3, year){
  
  x <- gdp %>% 
    dplyr::filter(`Country Code` == iso3) %>% 
    dplyr::pull(as.character(year)) %>% 
    as.numeric()
  return(x)
  
}
ICDF     <- function(x, theta, gini, gdppc){
  
  alpha = (1+gini)/(2*gini)
  thr <- (1-(1/alpha))*gdppc
  
  dens_pareto <-  VGAM::qpareto(p = x , 
                                scale = thr, 
                                shape = alpha)  #(x/ (theta+y))^alpha inverse cdf  #(alpha*theta*x^(alpha - 1))/((theta+x)^(alpha+1)) - density
  
  sd_log <- sqrt(2)*VGAM::probitlink(theta = (gini+1)/2)
  mean_log <-  log(gdppc)-((sd_log)^2/2)
  
  dens_logN <- qlnorm(p = x, 
                      meanlog = mean_log, 
                      sdlog = sd_log)
  
  return(dens_pareto^(0.32)*dens_logN^(1-0.32))
  
}
calc_AWE <- function(rwi_pth, country_gini, country_gdp){
  
  country_rwi <- readr::read_csv(rwi_pth) %>% 
    dplyr::select(rwi, longitude, latitude)
  
  country_rwi$rank <- rank(country_rwi$rwi, ties.method = "random")/(max(rank(country_rwi$rwi, ties.method = "random"))+1)
  country_rwi$icdf <- ICDF(x     = country_rwi$rank, 
                           theta = 1, 
                           gini  = country_gini, 
                           gdppc = country_gdp) #41986 - ocde mean
  
  country_rwi$AWE <- country_rwi$icdf*(country_gdp/(mean(country_rwi$icdf)) )
  cap <- quantile(country_rwi$AWE, probs = 0.98)
  country_rwi <- country_rwi %>% dplyr::mutate(AWE = ifelse(AWE > cap, cap, AWE))
  
  return(country_rwi)
  
}

f <- tibble(files = list.files(wealth_dir, full.names = T)[which(wealth_iso3 %in%avalaible_iso3 )],
            iso3 = wealth_iso3[which(wealth_iso3 %in%avalaible_iso3 )]) %>% 
  bind_cols(., lapply(.$iso3, get_gini, gini = gini ) %>% bind_rows) %>% 
  drop_na %>% 
  mutate(gdp =unlist(purrr::map2(.x = iso3, .y = year, get_gdp, gdp = gdp)),
         AWE = purrr::pmap(.l=list(files, gini, gdp), .f = calc_AWE))

country_rwi <- f$AWE %>% 
  bind_rows() %>% 
  dplyr::select(longitude, latitude, AWE)

shp_country <- world_shp[(world_shp@data$CONTINENT == "Africa"), ]

coordinates(country_rwi) <-  ~longitude+latitude
crs(country_rwi) <- "+proj=longlat +datum=WGS84 +no_defs"

r <- raster(resolution = 0.04166667, ext = extent(shp_country))
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
r <- raster::mask(r, shp_country)
r_f <- rasterize(country_rwi, r, fun = mean, field = "AWE")
r_f <- raster::resample(r_f, mask %>% raster::crop(., extent(shp_country)) )

writeRaster(r_f, paste0(wealth_dir, "/AWE_africa.tif"), overwrite = T)
