# Install necessary packages (if not already installed)
install.packages(c("geodata", "Recocrop", "tmap", "terra", "sf", "dplyr", "raster"))

# Load necessary libraries
library(geodata)
library(Recocrop)
library(tmap)
library(terra)
library(sf)
library(dplyr)
library(raster)

# Define working directory
work_dir <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES"
data_dir <- file.path(work_dir, "3_Data/1_Raw/Historical_Climate_data_1991_2020/")
output_dir <- file.path(work_dir, "3_Data/3_Outputs/suitability_maps")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

# Define the path to the AOI shapefile
aoi_shapefile <- "D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/igad_cluster_1_1/igad_cluster_1_1.shp"
aoi <- vect(aoi_shapefile)

# List Raster Files for tmax, tmin, and prec for the period 1991-2020
# Adjust patterns to match your filenames
tmax_files <- list.files(file.path(data_dir, "tmax_1991_2020"), pattern = "tmax", full.names = TRUE)
tmin_files <- list.files(file.path(data_dir, "tmin_1991_2020"), pattern = "tmin", full.names = TRUE)
prec_files <- list.files(file.path(data_dir, "prec_1991_2020"), pattern = "prec", full.names = TRUE)

# Load and stack rasters using terra
tmax_stack <- rast(tmax_files)
tmin_stack <- rast(tmin_files)
prec_stack <- rast(prec_files)

# Compute monthly averages across all years
mean_tmax <- tapp(tmax_stack, index = 1:12, fun = mean)
mean_tmin <- tapp(tmin_stack, index = 1:12, fun = mean)
mean_prec <- tapp(prec_stack, index = 1:12, fun = mean)

# Generate tavg from the mean tmax and tmin
mean_tavg <- (mean_tmax + mean_tmin) / 2

# Mask and crop to AOI
tavg <- mask(crop(mean_tavg, aoi), aoi)
rain <- mask(crop(mean_prec, aoi), aoi)

# SRename layers to match WorldClim format
names(tavg) <- paste0("KEN_w~avg_", 1:12)
names(rain) <- paste0("KEN_w~prec_", 1:12)


# Grass types to model
grass_types <- c("Cenchrus ciliaris L.", "Chloris gayana", "Eragrostis superba",
                 "Pennisetum clandestinum", "Eragrostis tef", 
                 "Dichanthium aristatum",
                 "Desmodium intortum",
                 "Mucuna pruriens", "Acacia senegal",
                 "Balanites aegyptiaca")

# Load place names for Kenya, Uganda, and South Sudan
ken_setl <- st_read("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/Place_Names/Kenya/hotosm_ken_populated_places_points_shp.shp") %>%
  st_set_crs(st_crs(aoi)) %>%
  st_as_sf() %>%
  st_intersection(st_as_sf(aoi)) %>%
  filter(place == 'town')

uga_setl <- st_read("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/Place_Names/Uganda/hotosm_uga_populated_places_points_shp.shp") %>%
  st_set_crs(st_crs(aoi)) %>%
  st_as_sf() %>%
  st_intersection(st_as_sf(aoi)) %>%
  filter(place == 'town')

ssd_setl <- st_read("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw/Place_Names/South_Sudan/ssd_pppls_ocha_20221216.shp") %>%
  st_set_crs(st_crs(aoi)) %>%
  st_as_sf() %>%
  st_intersection(st_as_sf(aoi))

# Step 7: Function to run the EcoCrop model and save plots
run_and_save_ecocrop <- function(grass_name, rain, tavg, aoi, ken_setl, uga_setl, ssd_setl) {
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
  
  # Clip the suitability index to the AOI
  hv_clipped <- crop(hv, aoi)
  
  # Create a subdirectory for each grass type
  grass_dir <- file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability_1991_2020", gsub(" ", "_", grass_name))
  dir.create(grass_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Convert the AOI and place names to sf for plotting
  aoi_sf <- st_as_sf(aoi)
  
  # Create the thematic map
  them_map <- tm_shape(hv_clipped) +
    tm_raster(style = "cat", palette = "-RdYlGn", title = "Suitability Index") +
    tm_shape(aoi_sf) +
    tm_borders(col = "red", lwd = 2) +
    tm_shape(ken_setl) + tm_text("name", size = 0.3, col = "black", remove.overlap = T) +
    tm_shape(uga_setl) + tm_text("name", size = 0.3, col = "black", remove.overlap = T) +
    tm_shape(ssd_setl) + tm_text("featureNam", size = 0.3, col = "black", remove.overlap = T) +
    tm_layout(
      title = paste("Suitability Map for", grass_name, "in Karamoja"),
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
  
  # Save the suitability map as JPEG
  tmap_save(them_map, file.path(grass_dir, paste0("suitability_map_", gsub(" ", "_", grass_name), ".jpg")), width = 1200, height = 800, dpi = 300)
  
  # Return the raster object for further processing
  return(hv_clipped)
}

# Step 8: Initialize a list to store suitability indices
suitability_list <- list()

# Loop through each grass type and run the model
for (grass_name in grass_types) {
  tryCatch({
    suitability_index <- run_and_save_ecocrop(grass_name, rain, tavg, aoi, ken_setl, uga_setl, ssd_setl)
    suitability_list[[grass_name]] <- suitability_index
  }, error = function(e) {
    message("Error in processing ", grass_name, ": ", e$message)
  })
}

# Stack all suitability indices
suitability_stack <- rast(suitability_list)

# Calculate the average suitability
avg_suitability <- app(suitability_stack, fun = mean, na.rm = TRUE)

# Calculate the median suitability
median_suitability <- app(suitability_stack, fun = median, na.rm = TRUE)

# Save the average suitability map
avg_suitability_path <- file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability_1991_2020/models", "average_suitability.tif")
writeRaster(avg_suitability, avg_suitability_path, overwrite = TRUE)

# Save the median suitability map
median_suitability_path <- file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability_1991_2020/models", "median_suitability.tif")
writeRaster(median_suitability, median_suitability_path, overwrite = TRUE)

# Plot the median suitability with additional details
median_suitability_jpg <- file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability_1991_2020/models", "median_suitability.jpg")

# Create an sf object for plotting
aoi_sf <- st_as_sf(aoi)

# Plot the median suitability map with added details
median_map <- tm_shape(median_suitability) +
  tm_raster(style = "cont", palette = "-RdYlGn", title = "Median Suitability") +
  tm_shape(aoi_sf) +
  tm_borders(col = "black", lwd = 2) +
  tm_shape(ken_setl) + tm_text("name", size = 0.3, col = "black", remove.overlap = TRUE) +
  tm_shape(uga_setl) + tm_text("name", size = 0.3, col = "black", remove.overlap = TRUE) +
  tm_shape(ssd_setl) + tm_text("featureNam", size = 0.3, col = "black", remove.overlap = TRUE) +
  tm_layout(
    title = "Karamoja Pasture Suitability Index 1991-2020",
    title.size = 3.6,
    title.fontface = "bold",  # 
    title.position = c("center", "top"),
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.title.size = 3.8,
    legend.text.size = 0.8,
    inner.margins = c(0.02, 0.02, 0.02, 0.02)
  ) +
  tm_legend(outside = TRUE, outside.position = "right") +
  tm_scale_bar(width = 0.15, position = c("left", "bottom")) +
  tm_compass(type = "8star", position = c("right", "top"), size = 1) +
  tm_graticules(n.x = 6, n.y = 6, lines = TRUE, labels.size = 0.6)


# Calculate median suitability values for each grass type
suitability_values <- sapply(suitability_list, function(r) median(values(r), na.rm = TRUE))

# Save the median suitability map as JPEG
tmap_save(median_map, median_suitability_jpg, width = 1200, height = 800, dpi = 300)

# Ensure the map is displayed
median_map

