#*Pasture Suitability Modelling to inform resource sharing agreement in Karamoja
#* 
#* 
#* Author:: Ogero Derrick
#* 
########################################################################
# Install necessary  packages
install.packages("geodata")
install.packages("Recocrop")
install.packages("tmap")
install.packages("terra")
install.packages("sf")
install.packages("dplyr")

# Load necessary libraries
library(geodata)
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

# Define ISO codes for the Karamoja cluster (Kenya, Uganda, South Sudan, Ethiopia)
iso_codes <- c("KEN", "UGA", "SSD", "ETH")

# Download climate data for the Karamoja cluster
rain_list <- lapply(iso_codes, function(iso) 
  geodata::worldclim_country(iso, var="prec", path="D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw"))
tavg_list <- lapply(iso_codes, function(iso) 
  geodata::worldclim_country(iso, var="tavg", path="D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/1_Raw"))

# Merge the climate data
rain <- do.call(merge, rain_list)
tavg <- do.call(merge, tavg_list)

# Mask the climate data to the AOI
rain_aoi <- mask(rain, aoi)
tavg_aoi <- mask(tavg, aoi)

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

# Function to run the EcoCrop model and save plots
run_and_save_ecocrop <- function(grass_name, rain, tavg, aoi, work_dir, ken_setl, uga_setl, ssd_setl) {
  # Create an Ecocrop model
  crop <- Recocrop::ecocropPars(grass_name)
  model <- Recocrop::ecocrop(crop)
  
  # Use the model to make predictions
  plant <- predict(model, prec=rain, tavg=tavg)
  
  # Process predictions to suitability index
  p <- classify(plant > 0, cbind(0, NA)) * 1:12
  pm <- median(p, na.rm=TRUE)
  hv <- pm + model$duration
  hv <- ifel(hv > 12, hv - 12, hv)
  
  # Normalize suitability values to 0-100 range
  hv <- (pm / 12) * 100
  print(hv)
  
  # Clip the suitability index to the AOI
  hv_clipped <- crop(hv, aoi)
  print(hv_clipped)
  # Create a subdirectory for each grass type
  grass_dir <- file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs","suitability", gsub(" ", "_", grass_name))
  dir.create(grass_dir, showWarnings = FALSE, recursive = TRUE)
  print(grass_dir)
  
  # Convert the AOI and place names to sf for plotting
  aoi_sf <- st_as_sf(aoi)
  
  # Create the thematic map
  them_map <- tm_shape(hv_clipped) +
    tm_raster(style = "cat", palette = "-RdYlGn", title = "Suitability Index",legend.show = TRUE) +
    tm_shape(aoi_sf) +
    tm_borders(col = "red", lwd = 2) +
    tm_shape(ken_setl) + tm_text("name", size = 0.3, col = "black", remove.overlap = T) +
    tm_shape(uga_setl) + tm_text("name", size = 0.3, col = "black", remove.overlap = T) +
    tm_shape(ssd_setl) + tm_text("featureNam", size = 0.3, col = "black",remove.overlap = T) +
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
  them_map
  # Save the suitability map as JPEG
  tmap_save(them_map, file.path(grass_dir,paste0("suitability_map_", gsub(" ", "_", grass_name), ".jpg")), width = 1200, height = 800, dpi = 300)
  
  # Return the raster object for further processing
  return(hv_clipped)
}

# Initialize a list to store suitability indices
suitability_list <- list()

# Loop through each grass type and run the model
for (grass_name in grass_types) {
  tryCatch({
    suitability_index <- run_and_save_ecocrop(grass_name, rain_aoi, tavg_aoi, aoi, 
                                              work_dir, ken_setl, uga_setl, ssd_setl)
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
avg_suitability_path <- file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability/models", "average_suitability.tif")
writeRaster(avg_suitability, avg_suitability_path, overwrite = TRUE)

# Save the median suitability map
median_suitability_path <- file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability/models", "median_suitability.tif")
writeRaster(median_suitability, median_suitability_path, overwrite = TRUE)

# Plot the median suitability with additional details
jpeg(file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability/models", "median_suitability.jpg"))

# Create an sf object for plotting
aoi_sf <- st_as_sf(aoi)

# Plot the median suitability map with added details
median_map <- tm_shape(median_suitability) +
  tm_raster(style = "cont", palette = "-RdYlGn", title = "Median Suitability") +
  tm_shape(aoi_sf) +
  tm_borders(col = "black", lwd = 2) +
  tm_shape(ken_setl) + tm_text("name", size = 0.3, col = "black", remove.overlap = T) +
  tm_shape(uga_setl) + tm_text("name", size = 0.3, col = "black", remove.overlap = T) +
  tm_shape(ssd_setl) + tm_text("featureNam", size = 0.3, col = "black",remove.overlap = T) +
  tm_layout(
    title = "Karamoja Pasture Suitability Index",
    title.size = 2.6,
    title.fontface = "bold",  # Make the title bold
    title.position = c("center", "top"),
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.title.size = 1.8,
    legend.text.size = 0.8,
    inner.margins = c(0.02, 0.02, 0.02, 0.02)
  ) +
  tm_legend(outside = TRUE, outside.position = "right") +
  tm_scale_bar(width = 0.15, position = c("left", "bottom")) +
  tm_compass(type = "8star", position = c("right", "top"), size = 1) +
  tm_graticules(n.x = 6, n.y = 6, lines = TRUE, labels.size = 0.6)

# Define a function to create a legend for the grass types
create_legend <- function(grass_types, suitability_values) {
  colors <- tmaptools::get_brewer_pal("RdYlGn", n = 7)
  bins <- seq(40, 75, length.out = 7)
  labels <- paste(grass_types, round(suitability_values, 1), sep = " - ")
  tm_add_legend(type = "fill", labels = labels, col = colors, title = "Grass Types with Median Suitability")
}

# Calculate median suitability values for each grass type
suitability_values <- sapply(suitability_list, function(r) median(values(r), na.rm = TRUE))

# Add the legend to the map
median_map <- median_map + create_legend(grass_types, suitability_values)

# Save the median suitability map as JPEG
tmap_save(median_map, file.path("D:/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/suitability/models", "median_suitability.jpg"), width = 1200, height = 800, dpi = 300)

median_map
