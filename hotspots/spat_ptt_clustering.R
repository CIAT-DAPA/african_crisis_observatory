options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, raster, trend))

root <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory'
source('https://raw.githubusercontent.com/CIAT-DAPA/african_crisis_observatory/main/base__lowest_gadm.R') # Get lowest administrative level per country

iso <- 'KEN'
fls <- list.files(path = paste0(root,'/data/',iso,'/conflict'), pattern = '*.tif$', full.names = T)
vrs <- gsub(pattern = '.tif', replacement = '', x = basename(fls), fixed = T)

testInteger <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if(test == TRUE){ return(TRUE) }
  else { return(FALSE) }
}

# Function to reclassify float rasters into integer rasters
reclass <- function(file = fl){
  
  r <- raster::raster(file) # Load raster
  if(testInteger(as.numeric(names(table(r[]))))){
    # Integer raster
    fnrast <- r %>% stars::st_as_stars()
    return(fnrast)
  } else {
    # Float raster
    qtl    <- raster::quantile(r, seq(.2,1,.2)) # Obtain quantiles
    rclmat <- data.frame(x1 = c(min(r[], na.rm = T),qtl[-length(qtl)]), x2 = qtl, y = 1:5) %>% as.matrix() # Reclassify matrix
    fnrast <- raster::raster(file) %>% raster::reclassify(., rclmat, include.lowest = T) %>% stars::st_as_stars(); rm(r,qtl,rclmat) # Final raster
    return(fnrast)
  }
  
}

# Reclassify conflict rasters
rst <- fls %>% purrr::map(.f = function(fl){reclass(fl)})

eco_data <- purrr::reduce(.x = rst, .f = stars:::c.stars, try_hard = T)

eco_signature <- motif::lsp_signature(eco_data, type = "incove", window = 30)
eco_dist      <- motif::lsp_to_dist(eco_signature, dist_fun = "jensen-shannon")

eco_hclust <- hclust(eco_dist, method = "ward.D2")
plot(eco_hclust, hang = -1, cex = 0.6)

clusters <- cutree(eco_hclust, k = 4)

eco_grid_sf <- motif::lsp_add_clusters(eco_signature, clusters)
tm_clu <- tmap::tm_shape(eco_grid_sf) +
  tmap::tm_polygons("clust", style = "cat", palette = "Set2", title = "Cluster:") +
  tmap::tm_layout(legend.position = c("LEFT", "BOTTOM"))
tm_clu

pacman::p_load(maptree)

optcl <- maptree::kgs(cluster = eco_hclust, diss = eco_dist, maxclust = 10)
plot(names(optcl), optcl, xlab = "# clusters", ylab = "penalty")
optcl[which(optcl == min(optcl))]

# eco_grid_sf %>%
#   sf::write_sf(paste0(root,'/data/',iso,'/_results/conflict_clusters.shp'))

# Irregular approach
shp <- raster::shapefile(paste0(root,'/data/',iso,'/_shps/',iso,'.shp'))
adm <- grep(pattern = '^NAME_', x = names(shp), value = T)
shp@data$key <- tolower(do.call(paste, c(shp@data[,adm], sep="-")))
shp@data$id  <- as.character(1:nrow(shp@data))

shp <- sf::st_as_sf(shp)

eco_signature2 <- motif::lsp_signature(eco_data, type = "incove", window = shp["id"])
eco_dist2      <- motif::lsp_to_dist(eco_signature2, dist_fun = "jensen-shannon")

eco_hclust2 <- hclust(eco_dist2, method = "ward.D2")
plot(eco_hclust2, hang = -1, cex = 0.6)

clusters2 <- cutree(eco_hclust2, k = 3)

eco_grid_sf2 <- motif::lsp_add_clusters(eco_signature2, clusters2, window = shp["id"])
tm_clu2 <- tmap::tm_shape(eco_grid_sf2) +
  tmap::tm_polygons("clust", style = "cat", palette = "Set2", title = "Cluster:") +
  tmap::tm_layout(legend.position = c("LEFT", "BOTTOM"))
tm_clu2

# eco_grid_sf2 %>%
#   sf::write_sf(paste0(root,'/data/',iso,'/_results/conflict_clusters_poly.shp'))

pacman::p_load(maptree)

optcl2 <- maptree::kgs(cluster = eco_hclust2, diss = eco_dist2, maxclust = 10)
plot(names(optcl2), optcl2, xlab = "# clusters", ylab = "penalty")
optcl2[which(optcl2 == min(optcl2))]
