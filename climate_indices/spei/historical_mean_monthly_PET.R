################################################################################
#*Compute the monthly historical mean of potential Evapotranspiration(PET)
#* 
#* 
#* Author: Victor Korir

################################################################################
#Read data
monthly_PET_dir <- ('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/
                    Data/intermediate/spei/monthly_evapotranspiration')
monthly_PET_files <- list.files(monthly_PET_dir,pattern = '.tif', full.names = T)
monthly_PET <- rast(lapply(monthly_PET_files,function(x) terra::rast(x)))
names(monthly_PET) <- sources(monthly_PET)

#Function to separate the monthly PET and average them
monthly_average <- function(raster_stack, month) {
  # Extract the rasters corresponding to the given month
  monthly_rasters <- raster_stack[[seq(month, nlyr(raster_stack), 12)]]
  # Calculate the mean raster for the given month
  mean(monthly_rasters)
}

#Compute the monthly averages and add them to a stack
monthly_averages <- terra::rast()
for (month in 1:12) {
  print(paste("Processing month:", month))
  monthly_avg <- monthly_average(monthly_PET, month)
  plot(monthly_avg)
  monthly_averages <- terra::`add<-`(monthly_averages, monthly_avg)
}
# Rename the layers in the stack and write them to a file
names(monthly_averages) <- (month.name)
terra::writeRaster(monthly_averages, '//alliancedfs.alliance.cgiar.org/
                   WS18_Afrca_K_N_ACO2/FCM/Data/intermediate/spei/historical_mean_monthly_PET.tif')

