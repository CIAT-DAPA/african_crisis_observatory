#Todo: Compute the mean of climate indices over the years

#Author: Victor Korir
########################################################################################
library(terra)
#Read all indices from the folder
dir <-'//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/climate_indices/KEN/annual'
indices_list <- list.files(dir, pattern="*.tif", full.names=TRUE)

#Define function to compute the average of each raster file
index_average <- function(file){
  index_raster <- terra::rast(file)
  mean <- mean(index_raster)
  return(mean)
}

# Compute the mean raster for each file and stack them together
averaged_indices <- do.call(c, lapply(indices_list, index_average))

#get layer names from the folder and rename the stacked means
layer_names <- gsub(".tif$", "", basename(indices_list))
names(averaged_indices) <- layer_names

# Export the stacked
terra::writeRaster(averaged_indices, paste0(dir, '/average_indices.tiff') )
