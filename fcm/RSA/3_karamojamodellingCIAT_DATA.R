# Install necessary packages (if not already installed)
install.packages(c("Recocrop", "tmap", "terra", "sf", "dplyr"))

# Load necessary libraries
library(Recocrop)
library(tmap)
library(terra)
library(sf)
library(dplyr)

# Define working directory
work_dir <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES"
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

# Define the path to the AOI shapefile
aoi_shapefile <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/igad_cluster_1_1/igad_cluster_1_1.shp"

# Load AOI shapefile
aoi <- vect(aoi_shapefile)

# Define the path to the in-house climate data (1991-2050)
climate_data_dir <- "D:/OneDrive - CGIAR/SA_Team/Data/suitability/Monthly_Average_1991_2050"

# Load the climate data
tasmin <- rast(file.path(climate_data_dir, "Tasmin_mean_monthly_1991_2050.tif"))
tasmax <- rast(file.path(climate_data_dir, "Tasmax_mean_monthly_1991_2050.tif"))
precipitation <- rast(file.path(climate_data_dir, "Precipitation_mean_monthly_1991_2050.tif"))

# Rename layers to match WorldClim format
names(tasmin) <- paste0("KEN_w~avg_", 1:12)
names(tasmax) <- paste0("KEN_w~avg_", 1:12)
names(precipitation) <- paste0("KEN_w~prec_", 1:12)

# Calculate average temperature (tavg)
tavg <- (tasmin + tasmax) / 2
names(tavg) <- paste0("KEN_w~avg_", 1:12)

# Align raster layers to AOI
tavg <- project(tavg, crs(aoi))
rain <- project(precipitation, crs(aoi))

# Mask and crop climate data to AOI
tavg <- crop(mask(tavg, aoi), aoi)
rain <- crop(mask(rain, aoi), aoi)

# Grass types to model
grass_types <- c("Cenchrus ciliaris L.", "Chloris gayana", "Eragrostis superba",
                 "Pennisetum clandestinum", "Eragrostis tef", 
                 "Dichanthium aristatum", "Desmodium intortum",
                 "Mucuna pruriens", "Acacia senegal",
                 "Balanites aegyptiaca")

# Function to run the EcoCrop model and save plots
run_and_save_ecocrop <- function(grass_name, tavg, rain, aoi) {
  # Create an Ecocrop model
  crop <- Recocrop::ecocropPars(grass_name)
  model <- Recocrop::ecocrop(crop)
  
  # Use the model to make predictions
  plant <- predict(model, prec = rain, tavg = tavg)
  
  # Process predictions to suitability index
  p <- classify(plant > 0, cbind(0, NA)) * 1:12
  pm <- median(p, na.rm = TRUE)
  hv <- pm + model$duration
  hv <- ifel(hv > 12, hv - 12, hv)
  
  # Normalize suitability values to 0-100 range
  hv <- (pm / 12) * 100
  hv_clipped <- crop(hv, aoi)
  
  # Filter out NA values
  hv_clipped <- mask(hv_clipped, aoi)
  
  # Ensure suitability map has non-NA values
  if (all(is.na(hv_clipped[]))) {
    print(paste("No valid suitability values for", grass_name, "- skipping."))
    return(NULL)
  }
  
  # Create a subdirectory for each grass type
  grass_dir <- file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability", gsub(" ", "_", grass_name))
  dir.create(grass_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Convert the AOI to sf for plotting
  aoi_sf <- st_as_sf(aoi)
  
  # Create the thematic map
  them_map <- tm_shape(hv_clipped) +
    tm_raster(style = "cont", palette = "-RdYlGn", title = "Suitability Index") +
    tm_shape(aoi_sf) +
    tm_borders(col = "black", lwd = 2) +
    tm_layout(
      title = paste("Suitability Map for", grass_name),
      title.size = 2.0,
      title.position = c("center", "top"),
      legend.outside = TRUE,
      legend.outside.position = "right",
      legend.title.size = 1.2,
      legend.text.size = 0.8,
      inner.margins = c(0.02, 0.02, 0.02, 0.02)  # Adjust margins
    ) +
    tm_scale_bar(width = 0.15, position = c("right", "bottom")) +
    tm_compass(type = "8star", position = c("right", "top"), size = 1) +
    tm_graticules(n.x = 6, n.y = 6, lines = TRUE, labels.size = 0.6)
  
  # Save the suitability map as a JPEG file
  tmap_save(them_map, file.path(grass_dir, paste0("suitability_map_", gsub(" ", "_", grass_name), ".jpg")), width = 1200, height = 800, dpi = 300)
  
  # Return the raster object for further processing
  return(hv_clipped)
}


# Run and save EcoCrop model outputs for all grass types
suitability_all <- list()
for (grass_type in grass_types) {
  suitability_all[[grass_type]] <- run_and_save_ecocrop(grass_type, tavg, rain, aoi)
}

# Calculate and plot median suitability across all grass types
median_suitability <- median(stack(suitability_all), na.rm = TRUE)

# Create a thematic map for the median suitability
median_map <- tm_shape(median_suitability) +
  tm_raster(style = "cont", palette = "-RdYlGn", title = "Median Suitability Index") +
  tm_shape(st_as_sf(aoi)) +
  tm_borders(col = "black", lwd = 2) +
  tm_layout(
    title = "Median Suitability Map",
    title.size = 2.0,
    title.position = c("center", "top"),
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.title.size = 1.2,
    legend.text.size = 0.8,
    inner.margins = c(0.02, 0.02, 0.02, 0.02)  # Adjust margins
  ) +
  tm_scale_bar(width = 0.15, position = c("right", "bottom")) +
  tm_compass(type = "8star", position = c("right", "top"), size = 1) +
  tm_graticules(n.x = 6, n.y = 6, lines = TRUE, labels.size = 0.6)

# Save the median suitability map as a JPEG file
tmap_save(median_map, file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability", "median_suitability_map.jpg"), width = 1200, height = 800, dpi = 300)
