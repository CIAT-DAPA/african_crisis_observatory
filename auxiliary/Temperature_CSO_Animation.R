#temperature animation
install.packages("gifski")
library(gifski)
library(terra)
library(magick)
library(sf)
library(tmap)
library(tmaptools)
library(raster)
library(rlist)
#set working directory
dir <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/CMIP6/Africa/'
setwd(dir)

#find matching files
matching_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", 'tas_', ".*"))
matching_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", 'tasmax', ".*"))
matching_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", 'tasmin', ".*"))
matching_files
#extract models
model1 <- terra::rast(matching_files[[1]])
model2 <- terra::rast(matching_files[[2]])
model3 <- terra::rast(matching_files[[3]])
model4 <- terra::rast(matching_files[[4]])
plot(model4[[1]])
model2
#resample to fine resolution
ref <- terra::rast(ncols=83, nrows=89, xmin=-23.9, xmax=59.5, ymin=-37.4, ymax=40.2, nlyrs=1, res=0.5)
ref
model1_resampled <- terra::resample(model1, ref, method="ngb")
model1_resampled
plot(model1_resampled[[1]])
model2_resampled <- terra::resample(model2, ref)
model3_resampled <- terra::resample(model3, ref)
model4_resampled <- terra::resample(model4, ref)
#crop to kenya
kenya <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/shps/KEN_GIT/kenya.shp'
ken <- read_sf(kenya)
crop <- function(x){
  model1_kenya <- terra::mask(x,ken)
  model1_kenya <- terra::crop(model1_kenya,ext(c(33,42.7,-6,6)))
  return(model1_kenya)
}
model1_kenya <- crop(model1_resampled)
model2_kenya <- crop(model2_resampled)
model3_kenya <- crop(model3_resampled)
model4_kenya <- crop(model4_resampled)
plot(model1_kenya[[1]])
#extract mean for each year
mean_2030 <- terra::mean(terra::mean(model1_kenya["2030"],model2_kenya["2030"],model3_kenya["2030"],model4_kenya["2030"]))
mean_2030
names(mean_2030) <- "2030"
mean_2031 <- terra::mean(terra::mean(model1_kenya["2031"],model2_kenya["2031"],model3_kenya["2031"],model4_kenya["2031"]))
names(mean_2031) <- "2031"
mean_2032 <- terra::mean(terra::mean(model1_kenya["2032"],model2_kenya["2032"],model3_kenya["2032"],model4_kenya["2032"]))
names(mean_2032) <- "2032"
mean_2033 <- terra::mean(terra::mean(model1_kenya["2033"],model2_kenya["2033"],model3_kenya["2033"],model4_kenya["2033"]))
names(mean_2033) <- "2033"
mean_2034 <- terra::mean(terra::mean(model1_kenya["2034"],model2_kenya["2034"],model3_kenya["2034"],model4_kenya["2034"]))
names(mean_2034) <- "2034"
mean_2035 <- terra::mean(terra::mean(model1_kenya["2035"],model2_kenya["2035"],model3_kenya["2035"],model4_kenya["2035"]))
names(mean_2035) <- "2035"
mean_2036 <- terra::mean(terra::mean(model1_kenya["2036"],model2_kenya["2036"],model3_kenya["2036"],model4_kenya["2036"]))
names(mean_2036) <- "2036"
mean_2037 <- terra::mean(terra::mean(model1_kenya["2037"],model2_kenya["2037"],model3_kenya["2037"],model4_kenya["2037"]))
names(mean_2037) <- "2037"
mean_2038 <- terra::mean(terra::mean(model1_kenya["2038"],model2_kenya["2038"],model3_kenya["2038"],model4_kenya["2038"]))
names(mean_2038) <- "2038"
mean_2039 <- terra::mean(terra::mean(model1_kenya["2039"],model2_kenya["2039"],model3_kenya["2039"],model4_kenya["2039"]))
names(mean_2039) <- "2039"
mean_2040 <- terra::mean(terra::mean(model1_kenya["2040"],model2_kenya["2040"],model3_kenya["2040"],model4_kenya["2040"]))
names(mean_2040) <- "2040"
mean_2041 <- terra::mean(terra::mean(model1_kenya["2041"],model2_kenya["2041"],model3_kenya["2041"],model4_kenya["2041"]))
names(mean_2041) <- "2041"
mean_2042 <- terra::mean(terra::mean(model1_kenya["2042"],model2_kenya["2042"],model3_kenya["2042"],model4_kenya["2042"]))
names(mean_2042) <- "2042"
mean_2043 <- terra::mean(terra::mean(model1_kenya["2043"],model2_kenya["2043"],model3_kenya["2043"],model4_kenya["2043"]))
names(mean_2043) <- "2043"
mean_2044 <- terra::mean(terra::mean(model1_kenya["2044"],model2_kenya["2044"],model3_kenya["2044"],model4_kenya["2044"]))
names(mean_2044) <- "2044"
mean_2045 <- terra::mean(terra::mean(model1_kenya["2045"],model2_kenya["2045"],model3_kenya["2045"],model4_kenya["2045"]))
names(mean_2045) <- "2045"
mean_2046 <- terra::mean(terra::mean(model1_kenya["2046"],model2_kenya["2046"],model3_kenya["2046"],model4_kenya["2046"]))
names(mean_2046) <- "2046"
mean_2047 <- terra::mean(terra::mean(model1_kenya["2047"],model2_kenya["2047"],model3_kenya["2047"],model4_kenya["2047"]))
names(mean_2047) <- "2047"
mean_2048 <- terra::mean(terra::mean(model1_kenya["2048"],model2_kenya["2048"],model3_kenya["2048"],model4_kenya["2048"]))
names(mean_2048) <- "2048"
mean_2049 <- terra::mean(terra::mean(model1_kenya["2049"],model2_kenya["2049"],model3_kenya["2049"],model4_kenya["2049"]))
names(mean_2049) <- "2049"
mean_2050 <- terra::mean(terra::mean(model1_kenya["2050"],model2_kenya["2050"],model3_kenya["2050"],model4_kenya["2050"]))
names(mean_2050) <- "2050"
plot(mean_2050)

#animation
means <- c(mean_2030,mean_2031,mean_2032,mean_2033,mean_2034,mean_2035,mean_2036,mean_2037,mean_2038,mean_2039,
           mean_2040,mean_2041,mean_2042,mean_2043,mean_2044,mean_2045,mean_2046,mean_2047,mean_2048,mean_2049,mean_2050)
means
kenya_adm0 <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/shps/KEN_GIT/Kenya_adm0.shp'
kenya_adm0 <- read_sf(kenya_adm0)
plot(kenya_adm0)
temp_plot <- tm_shape(means)+
  tm_raster(style='cont', palette=get_brewer_pal('-RdYlGn',n=9, plot=FALSE), 
            title='Daily Temp (°C)')+
  tm_shape(kenya_adm0) +
  tm_borders(col="black", lwd=1)+
  tm_layout(legend.position = c('left','bottom'),legend.text.size=1.3, 
           legend.title.size = 1.1, legend.title.fontface =2,
           legend.hist.width = 20,title.size = 10,main.title.size = 5,
           legend.width = 1) +
  tm_facets(nrow=1,ncol=1)
temp_plot
animation <- tmap_animation(temp_plot, dpi=400, 'Daily_temp.gif',asp = 0)





