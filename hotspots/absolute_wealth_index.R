### Absolute wealth index

pacman::p_load(raster, sp, sf, tidyverse, VGAM)


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

gdp <- read_csv("D:/african observatory/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488.csv")
gini <- read_csv("D:/african observatory/API_SI.POV.GINI_DS2_en_csv_v2_2445276/API_SI.POV.GINI_DS2_en_csv_v2_2445276.csv")
mask <- raster("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/_global/masks/mask_world_1km.tif")


ISO3 <- "SDN"

wealth_dir <- "D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/_global/wealth_index"

rwi_out_dir <-  paste0("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/", ISO3,"/wealth_index" )


if(!dir.exists(rwi_out_dir)){dir.create(rwi_out_dir)}


country_gini <- gini %>% 
  filter(`Country Code` == ISO3) %>% 
  dplyr::select(where(~!all(is.na(.)))) %>% 
  dplyr::pull(ncol(.))
country_gini <- country_gini/100


country_gdp <- gdp %>% 
  filter(`Country Code` == ISO3) %>% 
  dplyr::pull(`2009`)



shp_country <- raster::shapefile(paste0("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/", ISO3, "/_shps/",ISO3, ".shp"))
 

country_rwi <- read_csv(paste0("D:/african observatory/relative-wealth-index-april-2021/", ISO3, "_relative_wealth_index.csv")) %>% 
  dplyr::select(rwi, longitude, latitude)


country_rwi$rank <- rank(country_rwi$rwi, ties.method = "random")/(max(rank(country_rwi$rwi, ties.method = "random"))+1)



country_rwi$icdf <- ICDF(x = country_rwi$rank, 
                  theta = 1, 
                  gini = country_gini, 
                  gdppc = country_gdp) #41986 - ocde mean

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

writeRaster(r_f, paste0(rwi_out_dir, "/",ISO3, "_AWE.tif"), overwrite = T)

#######################################################
####### AWE para todo africa        ##################
#####################################################


world_shp <- raster::shapefile("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/_global/world_shapefile/all_country/all_countries.shp")
wealth_dir <- "D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/_global/wealth_index"
mask <- raster("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/_global/masks/mask_world_1km.tif")

wealth_iso3 <- list.files(wealth_dir) %>% 
  substring(., 1, 3)

avalaible_iso3 <- world_shp@data$ISO3[(world_shp@data$ISO3 %in% wealth_iso3) &  (world_shp@data$CONTINENT == "Africa")]


gdp <- read_csv("D:/african observatory/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2531488.csv")
gini <- read_csv("D:/african observatory/API_SI.POV.GINI_DS2_en_csv_v2_2445276/API_SI.POV.GINI_DS2_en_csv_v2_2445276.csv")

get_gini <- function(gini, iso3){
  x <- gini%>% 
    dplyr::filter(`Country Code` == iso3) %>% 
    dplyr::select(where(~!all(is.na(.)))) %>%
    dplyr::select(ncol(.)) 
  
  return(data.frame(year = as.numeric(names(x)),
              gini = as.numeric(x)/100))
}


get_gdp <- function(gdp, iso3, year){
  
  x <- gdp %>% 
    dplyr::filter(`Country Code` == iso3) %>% 
    dplyr::pull(as.character(year)) %>% 
    as.numeric()
  return(x)
  
}

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



calc_AWE<- function(rwi_pth, country_gini, country_gdp){
  
  
  country_rwi <- read_csv(rwi_pth) %>% 
    dplyr::select(rwi, longitude, latitude)
  
  
  country_rwi$rank <- rank(country_rwi$rwi, ties.method = "random")/(max(rank(country_rwi$rwi, ties.method = "random"))+1)
  
  
  
  country_rwi$icdf <- ICDF(x = country_rwi$rank, 
                           theta = 1, 
                           gini = country_gini, 
                           gdppc = country_gdp) #41986 - ocde mean
  
  country_rwi$AWE <- country_rwi$icdf*(country_gdp/(mean(country_rwi$icdf)) )
  cap <- quantile(country_rwi$AWE, probs = 0.98)
  country_rwi <- country_rwi %>% 
    dplyr::mutate(AWE = ifelse(AWE > cap, cap, AWE))
 
  
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

