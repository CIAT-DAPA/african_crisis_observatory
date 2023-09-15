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

#extract models
model1 <- terra::rast(matching_files[[1]])
model1
model2 <- terra::rast(matching_files[[2]])
model2
model3 <- terra::rast(matching_files[[3]])
model3
model4 <- terra::rast(matching_files[[4]])
model4

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
means <- list(mean_2030, mean_2031, mean_2032, mean_2033, mean_2034, mean_2035)
means


#write yearly means to files
output <- paste0(dir,'results/kenya/ssp585/')
output
mean_2030_f <-  terra::writeRaster(mean_2030,paste0(output, 'mean_2030.tif'))
mean_2031_f <- terra::writeRaster(mean_2031,paste0(output, 'mean_2031.tif'))
mean_2032_f <- terra::writeRaster(mean_2032,paste0(output, 'mean_2032.tif'))
mean_2033_f <- terra::writeRaster(mean_2033,paste0(output, 'mean_2033.tif'))
mean_2034_f <- terra::writeRaster(mean_2034,paste0(output, 'mean_2034.tif'))
mean_2035_f <- terra::writeRaster(mean_2035,paste0(output, 'mean_2035.tif'))


#cropping to kenya
kenya <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr1/1.Data/shps/KEN_GIT/kenya.shp'
ken <- read_sf(kenya)
plot(ken)

crop_fun <- function(img){
  data <- raster(img)
  masked <- raster::mask(data,ken)
  cropped <- raster::crop(masked,extent(c(31,44,-8,8)))
  return (cropped)
}
resample_fun <- function(img,res){
  temp <- rast(xmin=ext(img)$xmin, xmax=ext(img)$xmax, ymin=ext(img)$ymin, ymax=ext(img)$ymax, res=res, crs=crs(img))
  resampled_img <- terra::resample(img,temp)
  return(resampled_img)
}
resample_kenya <- resample_fun(raster(paste0(output, 'mean_2030.tif')),0.1)
kenya_2030 <- crop_fun(paste0(output, 'mean_2030.tif'))
kenya_2031 <- crop_fun(paste0(output, 'mean_2031.tif'))
kenya_2032 <- crop_fun(paste0(output, 'mean_2032.tif'))
kenya_2033 <- crop_fun(paste0(output, 'mean_2033.tif'))
kenya_2034 <- crop_fun(paste0(output, 'mean_2034.tif'))
kenya_2035 <- crop_fun(paste0(output, 'mean_2035.tif'))
plot(kenya_2030)
kenya_2031
#visualizing data

map1 <- tm_shape(kenya_2031)+
  tm_raster(style='cont', palette=get_brewer_pal("Blues", n=4), title = 'Precipitation')+
  tm_shape(ken) +
  tm_borders(col="black", lwd=1)+
  tm_compass(position = c("right", "top"))+
  tm_scale_bar(position = c("right", "bottom"))+
  tm_layout(main.title = "2031 Precipitation SSP585", 
            title.size = 1.5, title.position = c("left", "top"),legend.position =  c('left','bottom'))
map1
tmap_save(map1, filename="kenya_2031.jpg", height=8.5, width=11, units="in", dpi=300)

#animation
raster_layers <- lapply(1:nlyr(ssp_kenya), function(i) ssp_kenya[[i]])
tm_animation <- tm_shape(raster_layers[[1]]) +
  tm_raster() +
  tm_borders() +
  tm_text("Layer {frame}", x = 0.1, y = 0.05) +
  tmap_animation(frames = raster_layers, width = 800, height = 600)
tmraster <- tm_shape(ssp_kenya) +
  tm_raster()
tmraster
# Create an animation using tm_animation
animation <- tmap_animation(tmraster, width = 400, height = 400, duration = 2)
path <- paste0(dir,'results/kenya/')
dir.create(path)
dir.create(paste0(path,'ssp245/'))
name <- names(ssp_kenya[[1]])
name
new <- image_read(path)
writeRaster(ssp_245mean,paste0(path, 'real.tif'))
no <- nlyr(ssp_245mean)
for (x in 1:no){
  name <- names(ssp_245mean[[x]])
  name <- substr(name,1,10)
  writeRaster(ssp_245mean[[x]],paste0(path, name,'.tiff'))
}
first <- image_read(paste0(dir,'new1.tif'))
second <- image_read(paste0(dir,'new2.tif'))
plot(first)
combined <- c(first,second)
combined
image_scale(combined,"500x500")
image_info(combined)
aniamtion <- image_animate(combined,fps=1)
aniamtion
output <- paste0(path,'animation.gif')
image_write(aniamtion,output)
tiff_files <- list.files(path, pattern = ".tiff$", full.names = TRUE)
tiff_files
images <- image_read(tiff_files)

image_animate(images)
