#Compute monthly averages for CMIP6 data for precipitation and temperature
#Author:Brenda Chepngetich 2024

#load libraries
library(raster)
library(terra)
#set working directory
wd <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/Data/CMIP6/Monthly_1991_2050/"
setwd(wd)

#regex for hitorical data
hist_pr <- list.files(pattern=paste0(".*", 'historical', ".*", 'pr', ".*"))
hist_tasmax <- list.files(pattern=paste0(".*", 'historical', ".*", "tasmax", ".*"))
hist_tasmin <- list.files(pattern=paste0(".*", 'historical', ".*", "tasmin", ".*"))
hist_tasmax
#regex for future files
pr_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", 'pr', ".*"))
tasmax_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", "tasmax", ".*"))
tasmin_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", "tasmin", ".*"))

#find raster with finest resolution to use as reference for resampling between models
find_finest <- function(files){
  finest <- terra::rast(files[[1]])
  smallest_area <- prod(res(finest))
  for (f in 2:length(files)) {
    current_finest <- terra::rast(files[[f]])
    current_smallest <- prod(res(current_finest))
    if (current_smallest < smallest_area) {
      finest <- current_finest
      smallest_area <- current_smallest
    }
  }
  finest
}
#future data
pr_ref <- find_finest(pr_files)
tasmax_ref <- find_finest(tasmax_files)
tasmin_ref <- find_finest(tasmin_files)
tasmin_ref
#historical data
hist_pr_ref <- find_finest(hist_pr)
hist_tasmax_ref <- find_finest(hist_tasmax)
hist_tasmin_ref <- find_finest(hist_tasmin)
#create raster stacks for each variable
create_stack <- function (files, reference){
  stck <- terra::rast()
  for (x in files) {
    model <- terra::rast(x)
    model <- terra::resample(model, reference, method="near")
    stck <- terra::`add<-`(stck, model)
  }
  stck
}
pr_stack <- create_stack(pr_files, pr_ref)
tasmin_stack <- create_stack(tasmin_files, tasmin_ref)
tasmax_stack <- create_stack(tasmax_files,tasmax_ref)
hist_pr_stack <- create_stack(hist_pr, hist_pr_ref)
hist_tasmin_stack <- create_stack(hist_tasmin, hist_tasmin_ref)
hist_tasmax_stack <- create_stack(hist_tasmax,hist_tasmax_ref)
#functions to compute averages
monthly_avg <- function(raster_stack, month) {
  # Extract the rasters corresponding to the given month
  monthly_rasters <- raster_stack[[seq(month, nlyr(raster_stack), 12)]]
  # Calculate the mean raster for the given month
  mean(monthly_rasters)
}

compute_avg <- function(stack_) {
  monthly_stack <- terra::rast()
  for (month in 1:12){
    cat("calculating monthly averages..")
    month_avg <- monthly_avg(stack_,month)
    monthly_stack <- terra::`add<-`(monthly_stack, month_avg)
    }
  monthly_stack  
}
#precipitation
prep_avg <- compute_avg(pr_stack)
hist_prep_avg <- compute_avg(hist_pr_stack)
prep_avg
names(precipitation) <- (month.name)
prep <- terra::rast()
prep <- terra::`add<-`(prep,hist_prep_avg)
prep <- terra::`add<-`(prep,prep_avg)
precipitation <- compute_avg(prep)
plot(prep_mm[[1]])
#convert to mm
prep_mm <- precipitation * 86400
#tasmin
tasmin_avg <- compute_avg(tasmin_stack)



#tasmax
tasmax_avg <- compute_avg(tasmax_stack)
names(tasmax) <- (month.name)
hist_tasmax_avg <- compute_avg(hist_tasmax_stack)
hist_tasmax_avg
tasmax <- terra::rast()
tasmax <- terra::`add<-`(tasmax, hist_tasmax_avg)
tasmax <- terra::`add<-`(tasmax, tasmax_avg)
tasmax <- compute_avg(tasmax)
tasmax

#write outputs to file
terra::writeRaster(prep_mm, 
                   '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/Data/CMIP6/Monthly_Average_1991_2050/Precipitation_mean_monthly_1991_2050.tif')
