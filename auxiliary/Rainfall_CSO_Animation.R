install.packages('rlist')
#load required packages
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
matching_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", 'pr', ".*"))
matching_files
length(matching_files)
hist<- terra::rast(paste0(dir,"CMIP6_ACCESS-ESM1-5_historical_r1i1p1f1_pr_Africa_daily.tif"))
hist[[1]]
hist["2015-01-01 12:00:00 GMT"]
#extract models
model1 <- terra::rast(matching_files[[1]])
model2 <- terra::rast(matching_files[[2]])
model3 <- terra::rast(matching_files[[3]])
model4 <- terra::rast(matching_files[[4]])
plot(model1[[1]])

#resampling to have same resolution 
first_resampled <- terra::resample(model1,model3)
first_resampled
second_resampled <- terra::resample(model2,model3)
second_resampled
fourth_resampled <- terra::resample(model4,model3)
fourth_resampled

#crop to have same extents, 
#find object with the least extent to use as reference for cropping
first_area <- terra::expanse(first_resampled[[1]],unit='km')
first_area
second_area <- terra::expanse(second_resampled[[1]],unit='km')
second_area
third_area <- terra::expanse(model3[[1]],unit='km')
third_area
fourth_area <- terra::expanse(fourth_resampled[[1]],unit='km')
fourth_area
min(first_area,second_area,third_area, fourth_area)
min_extent <- fourth_resampled

#cropping
first_cropped <- terra::crop(first_resampled,fourth_resampled)
first_cropped
second_cropped <- terra::crop(second_resampled,fourth_resampled)
second_cropped
third_cropped <- terra::crop(model3,fourth_resampled)
third_cropped
first <- first_cropped
second <- second_cropped
third <- third_cropped
fourth <- fourth_resampled
crs(first)

#find mean between models
model_mean <- terra::mean(first,second,third, fourth)
model_mean
plot(model_mean[[1]])
crs(model_mean)

#get yearly means

mean_2030 <- terra::mean(model_mean["2030"])
mean_2031 <- terra::mean(model_mean["2031"])
mean_2032 <- terra::mean(model_mean["2032"])
mean_2033 <- terra::mean(model_mean["2033"])
mean_2034 <- terra::mean(model_mean["2034"])
mean_2035 <- terra::mean(model_mean["2035"])
plot(mean_2032)
names(mean_2030) <- "2030"
mean_2030
#write yearly means to files
output <- paste0(dir,'results/kenya/ssp585/')
output
mean_2030_f <-  terra::writeRaster(mean_2030,paste0(output, 'mean_2030.tif'))
mean_2031_f <- terra::writeRaster(mean_2031,paste0(output, 'mean_2031.tif'))
mean_2032_f <- terra::writeRaster(mean_2032,paste0(output, 'mean_2032.tif'))
mean_2033_f <- terra::writeRaster(mean_2033,paste0(output, 'mean_2033.tif'))
mean_2034_f <- terra::writeRaster(mean_2034,paste0(output, 'mean_2034.tif'))
mean_2035_f <- terra::writeRaster(mean_2035,paste0(output, 'mean_2035.tif'))

#resampling and cropping to kenya
kenya <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/shps/KEN_GIT/kenya.shp'
ken <- read_sf(kenya)
plot(ken)
reference <- raster(ncol=83, nrow=89, xmn=-23.9, xmx=59.5, ymn=-37.4, ymx=40.2)
reference
res(reference) <- 0.5
res(reference)
reference
resample_function <- function(img){
  data <- raster::raster(img)
  data
  resampled <- raster::resample(data,reference, method="ngb")
  masked <- raster::mask(resampled, ken)
  cropped <- raster::crop(masked,extent(c(31,44,-8,8)))
  return (cropped)
}
ssp585_2030 <- resample_function(paste0(output, 'mean_2030.tif'))
names(ssp585_2030) <- "2030 Precipitation"
plot(ssp585_2030)
ssp585_2030
ssp585_2031 <- resample_function(paste0(output, 'mean_2031.tif'))
plot(ssp585_2031)
names(ssp585_2031) <- '2031'
ssp585_2032 <- resample_function(paste0(output, 'mean_2032.tif'))
plot(ssp585_2032)
names(ssp585_2032) <- '2032'
ssp585_2032
ssp585_2033 <- resample_function(paste0(output, 'mean_2033.tif'))
names(ssp585_2033) <- ' 2033'
ssp585_2033
ssp585_2034 <- resample_function(paste0(output, 'mean_2034.tif'))
ssp585_2034
ssp585_2035 <- resample_function(paste0(output, 'mean_2035.tif'))
ssp585_2035
ssp585_2030
st <- raster::stack(ssp585_2030,ssp585_2031)
plot(st)

#visualizing data
stack <- raster::stack(ssp585_2030,ssp585_2031,ssp585_2032,ssp585_2033,ssp585_2034,ssp585_2035)
plot(stack)
kenya_adm0 <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/shps/KEN_GIT/Kenya_adm0.shp'
kenya_adm0 <- read_sf(kenya_adm0)
kenya_adm0
plot(kenya_adm0)
tmraster <- tm_shape(stack) + 
  tm_raster(style='cont', palette=get_brewer_pal('RdYlGn',n=9), title='Precipitation flux (mm/s)') +
  tm_shape(kenya_adm0) +
  tm_borders(col="black", lwd=1)+
  tm_compass(position = c("right", "top"))+
  tm_scale_bar(position = c("right", "bottom"))+
  tm_layout(legend.position = c('left','bottom')) +
  tm_facets(nrow=1,ncol=1)
tmraster
# Create an animation using tm_animation
animation <- tmap_animation(tmraster, dpi=400, 'ssp585_Precipitation.gif')
# Display the animation
animation
# Calculate the maximum and minimum values from the raster 
max_value <- max(stack)
max_value
min_value <- min(stack)
min_value
# Create a common legend
legend <- tm_legend(title = "Common Legend",labels = c(min_value, max_value),at = c(min_value, max_value),
          labels.labels = c("Min Value", "Max Value"), col = "viridis",auto.palette.mapping = FALSE)

map1 <- tm_shape(kenya_2035)+
  tm_raster(style='cont', palette=get_brewer_pal("Blues", n=4), title = 'Precipitation')+
  tm_shape(ken) +
  tm_borders(col="black", lwd=1)+
  tm_compass(position = c("right", "top"))+
  tm_scale_bar(position = c("right", "bottom"))+
  tm_layout(main.title = "2035 Precipitation SSP585", 
            title.size = 1.5, title.position = c("left", "top"),legend.position =  c('left','bottom'))
map1
tmap_save(map1, filename="kenya_2031.jpg", height=8.5, width=11, units="in", dpi=300)


#animation using magick 
first <- image_read(paste0(dir,'kenya_2030.jpeg'))
second <- image_read(paste0(dir,'kenya_2031.jpeg'))
third <- image_read(paste0(dir,'kenya_2032.jpeg'))
fourth <- image_read(paste0(dir,'kenya_2033.jpeg'))
fifth <- image_read(paste0(dir,'kenya_2034.jpeg'))
sixth <- image_read(paste0(dir,'kenya_2035.jpeg'))
plot(first)
combined <- c(first,second, third, fourth, fifth,sixth)
combined
image_scale(combined,"500x500")
image_info(combined)
animation <- image_animate(combined,fps=1)
animation
output <- paste0(dir,'animation.gif')
image_write(animation,output)


#temperature animation

