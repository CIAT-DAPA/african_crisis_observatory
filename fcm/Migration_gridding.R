library(geodata)
library(sf)
library(terra)
library(tmap)
library(geodata)
library(tidyverse)
library(raster)
library(data.table)
#Shapefiles of individual countries in the region
DRC <- geodata::gadm(country = 'Democratic Republic of the Congo', level = 0, path = tempdir())
Uganda <- geodata::gadm(country = 'Uganda', level = 0, path = tempdir())
Ethiopia <- geodata::gadm(country = 'Ethiopia', level = 0, path = tempdir())
South_Sudan <- geodata::gadm(country = 'SSD', level = 0, path = tempdir())
Tanzania <- geodata::gadm(country = 'Tanzania ', level = 0, path = tempdir())
Kenya <- geodata::gadm(country = 'Kenya ', level = 0, path = tempdir())
Rwanda <- geodata::gadm(country = 'Rwanda ', level = 0, path = tempdir())
Burundi <- geodata::gadm(country = 'Burundi ', level = 0, path = tempdir())
Somalia <- geodata::gadm(country = 'Somalia ', level = 0, path = tempdir())
Djibouti <- geodata::gadm(country = 'Djibouti ', level = 0, path = tempdir())
Eritrea <- geodata::gadm(country = 'Eritrea ', level = 0, path = tempdir())


#merging the individual country boundaries
region <- rbind(DRC, Uganda, Ethiopia, South_Sudan, Tanzania, Kenya, Rwanda, Burundi, Somalia)
region <-st_as_sf(region)
region$COUNTRY[region$COUNTRY=="Dem.Rep.Congo"] <- "D.R.Congo"


# Create a grid covering the region
EA_region <- rast(ext = ext(region), res = 1.6)
EA_region <- terra::crop(EA_region, region)
values(EA_region) <- 1:ncell(EA_region)
EA_region <- terra::mask(EA_region, region)


# Define function to read and process geocoded data
read_and_process_geocode <- function(year) {
  file_path <- paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/Migration/geo_", year, ".csv")
  origins_geocode <- fread(file_path)
  
  mig_points <- vect(origins_geocode, geom=c('lon', 'lat'), crs = crs(EA_region))
  
  # Remove unnecessary objects from memory
  rm(origins_geocode)
  
  mig_raster <- terra::rasterize(mig_points, EA_region, fun = 'count')
  
  # Remove unnecessary objects from memory
  rm(mig_points)
  
  return(mig_raster)
}

# Read and process geocoded data for each year
years <- 2018:2023
mig_rasters <- lapply(years, read_and_process_geocode)

# Create a list of raster names
mig_raster_names <- paste("mig_rasters_", years, sep = "")

# Combine the raster objects into a single list
mig_rasters <- setNames(mig_rasters, mig_raster_names)

# Optionally, you can assign the individual raster objects to the global environment
list2env(mig_rasters, envir = .GlobalEnv)


#Conflict data

# Define a function to read and filter conflict data
read_and_filter_conflict <- function(country_code) {
  file_path <- sprintf("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/%s/conflict/%s_conflict.csv", country_code, country_code)
  confl_data <- read.csv(file_path)
  confl_data <- confl_data[confl_data$YEAR >= 2018,]
  return(confl_data)
}

# List of country codes
#, "SSD", "TZA", "UGA", 'COD'
country_codes <- c("ETH", "SSD", "TZA", "UGA", 'COD','KEN','SOM')

# Read and filter conflict data for each country
conflict_data_list <- lapply(country_codes, read_and_filter_conflict)

# Combine the data into a single dataframe
merged_conf <- bind_rows(conflict_data_list)

# Group by YEAR and nest the data
conf_grouped <- merged_conf %>% group_by(YEAR) %>% nest()

# Define a function to rasterize the data
rasterize_data <- function(data) {
  vect_data <- vect(data, geom = c('LONGITUDE', 'LATITUDE'), crs = crs(EA_region))
  terra::rasterize(vect_data, EA_region, fun = 'count')
}

# Apply the rasterization function to each nested dataset
conf_rasters <- map(conf_grouped$data, rasterize_data)
names(conf_rasters) <- paste0("conf_raster_", years)

#Climate data

ETH_TR <- subset(rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/ETH/climatic_indexes/season_type_1/TR.tif"), c( "2018" ,"2019", '2020', '2021', '2022'))
SSD_TR <- subset(rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/SSD/climatic_indexes/season_type_1/TR.tif"), c( "2018" ,"2019", '2020', '2021', '2022'))
UGA_TR <- subset(rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/UGA/climatic_indexes/season_type_1/TR.tif"), c( "2018" ,"2019", '2020', '2021', '2022'))
TZA_TR <- subset(rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/TZA/climatic_indexes/season_type_1/TR.tif"), c( "2018" ,"2019", '2020', '2021', '2022'))
COD_TR <- subset(rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/COD/climatic_indexes/seanon_type_1/TR.tif"), c( "2018" ,"2019", '2020', '2021', '2022'))

TR <- terra::merge(ETH_TR, SSD_TR, UGA_TR, TZA_TR, COD_TR)
rm(ETH_TR, SSD_TR, UGA_TR, TZA_TR, COD_TR)


ETH_AT <- subset(rast('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/ETH/climatic_indexes/season_type_1/AT.tif'), c( "2018" ,"2019", '2020', '2021', '2022'))
SSD_AT <- subset(rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/SSD/climatic_indexes/season_type_1/AT.tif"), c( "2018" ,"2019", '2020', '2021', '2022'))
UGA_AT <- subset(rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/UGA/climatic_indexes/season_type_1/AT.tif"), c( "2018" ,"2019", '2020', '2021', '2022'))
TZA_AT <- subset(rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/TZA/climatic_indexes/season_type_1/AT.tif"), c( "2018" ,"2019", '2020', '2021', '2022'))
COD_AT <- subset(rast("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/IOM/COD/climatic_indexes/seanon_type_1/AT.tif"), c( "2018" ,"2019", '2020', '2021', '2022'))

AT <- terra::merge(ETH_AT, SSD_AT, UGA_AT, TZA_AT, COD_AT)
rm(ETH_AT, SSD_AT, UGA_AT, TZA_AT, COD_AT)

AT <- terra::resample(AT, EA_region, method = 'bilinear')
TR <- terra::resample(TR, EA_region, method = 'bilinear')

#Preparing data for correlation

conf_rasters <- c(conf_rasters$conf_raster_2018,conf_rasters$conf_raster_2019,
                  conf_rasters$conf_raster_2020, conf_rasters$conf_raster_2021, conf_rasters$conf_raster_2022)
mig_rasters<- c(mig_rasters_2018, mig_rasters_2019, mig_rasters_2020, 
                mig_rasters_2021, mig_rasters_2022)


AT_mig <- stack(as(AT, 'Raster'), as(mig_rasters, 'Raster'))
TR_mig <- stack(as(TR, 'Raster'), as(mig_rasters, 'Raster'))
conf_mig <- stack(as(conf_rasters, 'Raster'), as(mig_rasters, 'Raster'))

#Pixelwise correlation fuction
gridcorts <- function(rasterstack, method, type=c("corel","pval","both")){
  # Values for (layers, ncell, ncol, nrow, method, crs, extent) come straight from the input raster stack
  # e.g. nlayers(rasterstack), ncell(rasterstack)... etc.
  print(paste("Start Gridcorts:",Sys.time()))
  print("Loading parameters")
  layers=nlayers(rasterstack);ncell=ncell(rasterstack);
  ncol=ncol(rasterstack);nrow=nrow(rasterstack);crs=crs(rasterstack);
  extent=extent(rasterstack);pb = txtProgressBar(min = 0, max = ncell, initial = 0)
  print("Done loading parameters")
  mtrx <- as.matrix(rasterstack,ncol=layers)
  empt <- matrix(nrow=ncell, ncol=2)
  print("Initiating loop operation")
  if (type == "corel"){
    for (i in 1:ncell){
      setTxtProgressBar(pb,i)
      if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
        empt[i,1] <- NA 
      } else 
        if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
          empt[i,1] <- NA 
        } else 
          empt[i,1] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$estimate)
    }
    print("Creating empty raster")
    corel <- raster(nrows=nrow,ncols=ncol,crs=crs)
    extent(corel) <- extent
    print("Populating correlation raster")
    values(corel) <- empt[,1]
    print(paste("Ending Gridcorts on",Sys.time()))
    corel
  } 
  else
    if (type == "pval"){
      for (i in 1:ncell){
        setTxtProgressBar(pb,i)
        if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
          empt[i,2] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
            empt[i,2] <- NA 
          } else 
            empt[i,2] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$p.value)
      }
      pval <- raster(nrows=nrow,ncols=ncol,crs=crs)
      extent(pval) <- extent
      print("Populating significance raster")
      values(pval) <- empt[,2]
      print(paste("Ending Gridcorts on",Sys.time()))
      pval
    }
  else
    if (type == "both"){
      for (i in 1:ncell){
        setTxtProgressBar(pb,i)
        if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
          empt[i,] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
            empt[i,] <- NA 
          } else {
            empt[i,1] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$estimate) 
            empt[i,2] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$p.value)
          }
      }
      c <- raster(nrows=nrow,ncols=ncol,crs=crs)
      p <- raster(nrows=nrow,ncols=ncol,crs=crs)
      print("Populating raster brick")
      values(c) <- empt[,1]
      values(p) <- empt[,2]
      brk <- brick(c,p)
      extent(brk) <- extent
      names(brk) <- c("Correlation","Pvalue")
      print(paste("Ending Gridcorts on",Sys.time()))
      brk
    }
}

#Correlations
mig_conf_corr <- gridcorts(conf_mig, method="pearson", type="both")
mig_AT_corr <- gridcorts(AT_mig, method="pearson", type="both")
mig_TR_corr <- gridcorts(TR_mig, method="pearson", type="both")



###################################################
#plotting 
tm_shape(mig_TR_corr)+
  tm_raster(palette="RdYlGn", title = '', style="jenks", colorNA = 'grey', legend.show=T,legend.reverse = T)+
  tm_shape(region)+
  tm_borders()+
  tm_graticules(ticks = T, lines = F)+
  tm_compass(type="arrow", position=c('left', 'top'))+
  tm_facets( free.scales = TRUE) +
  tm_layout(legend.position= c("right", "bottom"),
            legend.outside = F,
            legend.width = 2, 
            legend.text.size = 0.6)


#plot the conf rasters and migration rasters(grided data)
migrasters1 <- stack(as(mig_rasters_2018, 'Raster'), 
                     as(mig_rasters_2019, 'Raster'),
                     as(mig_rasters_2020, 'Raster'),
                     as(mig_rasters_2021, 'Raster'),
                     as(mig_rasters_2022, 'Raster'),
                     as(mig_rasters_2023, 'Raster'))
region_cast <- st_cast(region)
names(mig_rasters1) <- c("2018", "2019", "2020", "2021", "2022", "2023")
names(conf_rasters) <- c("2018", "2019", "2020", "2021", "2022", "2023")
mig_rasters <- mig_rasters/1000
migration_plot <- tm_shape(conf_rasters)+
  tm_raster(palette= 'YlOrBr',style = 'jenks', colorNA = 'white', title = 'Conflict Events')+
  tm_shape(region_cast)+
  tm_borders(col = "black")+
  tm_text('COUNTRY', size = 0.4, remove.overlap = T, just = 'top')+
  tm_compass( position=c('left', 'top'))+
  tm_graticules(ticks = T, lines = F)+
  tm_scale_bar(breaks = c(0, 250, 500), text.size = 0.5, 
               position = c("left", "bottom"))+
  tm_layout(legend.outside=T,
            legend.show = T,
            legend.outside.position = 'right',
            legend.outside.size = 0.15,
            #legend.title.size= 1,
            legend.frame=F,
            legend.just = c("left", "top"), 
            legend.position  = c("left", "bottom"),
            inner.margins = c(0.02, 0.02, 0.01, 0.02)
  )

out <- "C:/Users/vkorir/Documents/CS_CGIAR/"
tmap_save(migration_plot,  dpi= 600,  height=4, width=8, units="in",
          filename=paste0(out, "conf.png"))


#plot the correlation results()
plot_list <- list()
# Iterate over layers of the raster stack
for (i in 1:nlayers(mig_AT_corr)) {
  
  # Create the tm_shape object for the current layer
  current_plot <- tm_shape(mig_AT_corr[[i]]) +
    tm_raster(palette = "RdYlGn", style = "jenks", colorNA = 'grey', legend.show = TRUE, legend.reverse = TRUE) +
    tm_shape(region) +
    tm_borders() +
    tm_graticules(ticks = TRUE, lines = FALSE) +
    tm_compass(type = "arrow", position = c('left', 'top')) +
    tm_layout(legend.position = c("right", "bottom"),
              legend.outside = FALSE,
              legend.width = 1,
              legend.title.size = 0.6,
              legend.frame = T,
              legend.text.size = 0.5)
  
  # Add the current plot to the list
  plot_list[[i]] <- current_plot
}
coorp <- tmap_arrange(plot_list, ncol = 2)
tmap_save(coorp, dpi=600, height=3.11,width = 8,units = 'in', filename=paste0(out, "migATcorrplot.png"))

region<-tm_shape(region_cast)+
  tm_borders(col = 'brown')+
  tm_text('COUNTRY',size= 0.8, remove.overlap = T)+
  tm_graticules(ticks = TRUE, lines = FALSE) 
tmap_save(region, dpi = 600, height=4, width=8, units="in",
          filename=paste0(out, "region.png"))