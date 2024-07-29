#Compute monthly averages for CMIP6 data for precipitation and temperature
#Author:Brenda Chepngetich
#set working directory
library(raster)
wd <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/Data/CMIP6/monthly/"
setwd(wd)

#regex for files
pr_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", 'pr', ".*"))
tasmax_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", "tasmax", ".*"))
tasmin_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", "tasmin", ".*"))

#define the reference extent
ref <- terra::rast(ncols=83, nrows=89, xmin=-23.9, xmax=59.5, ymin=-37.4, ymax=40.2, nlyrs=1, res=0.5)

#create raster stacks for each variable
create_stack <- function (files){
  stck <- terra::rast()
  for (x in files) {
    pr_model <- terra::rast(x)
    pr_model <- terra::resample(pr_model, ref, method="near")
    stck <- terra::`add<-`(stck, pr_model)
  }
  stck
}
pr_stack <- create_stack(pr_files)
tasmin_stack <- create_stack(tasmin_files)
tasmax_stack <- create_stack(tasmax_files)

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
names(prep_avg) <- (month.name)
#convert to mm
prep_mm <- prep_avg * 86400
#tasmin
tasmin_avg <- compute_avg(tasmin_stack)
names(tasmin_avg) <- (month.name)
#tasmax
tasmax_avg <- compute_avg(tasmax_stack)
names(tasmax_avg) <- (month.name)

#write outputs to file
terra::writeRaster(prep_mm, 
                   '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/Data/CMIP6/Monthly Average/Precipitation_mean_monthly.tif')

