
#* TODO: Grid number of migrants originating from every pixel
#*  and perform pixelwise correlation with climate and conflict data
#* 
#* Author: Victor Korir and Benson Kenduiywo
################################################################################
rm(list=ls(all=TRUE))
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
list.of.packages <- c("tidyverse","sf","terra","geodata", "tmap", "raster", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/'


countries <- c('Democratic Republic of the Congo', 'Uganda', 'Ethiopia', 'South Sudan', 'Tanzania', 'Kenya',
               'Rwanda', 'Burundi', 'Somalia', 'Djibouti', 'Eritrea')
region <- gadm(country = countries , level = 0, path = paste0(root,'raw/admin'), version="latest")
region <-st_as_sf(region)

#'  Create vector square grids of approximately 51 km2 (0.46 degrees)  covering the region

grd <- st_make_grid(st_bbox(extent(region)+2), cellsize = 0.46, square = TRUE)
grd <- st_as_sf(grd)
grd$id <- 1:nrow(grd)
st_crs(grd) <- st_crs(region)
plot(st_geometry(grd))
plot(st_geometry(region), add = TRUE)
 
#' =======================================================================================
#' 1.0 MIGRATION PRESENCE in GRIDS (MEGA-PIXELS)
#' =======================================================================================
#' Load migration geocoded points

load_mig <- function(year){
  filename <- paste0(root,"intermediate/fms_geocoded/geo_", year, ".csv")
  temp <- read.csv(filename, header = T)
  names(temp) <- c('Index', 'YEAR', 'Reason', 'Forcibly_displaced', 'District', 'City', 'ID', 'Long', 'Lat', 'Address' )
  return(temp)
}

years <- 2018:2023
mig <- lapply(years, load_mig)
mig <- do.call(rbind, mig)
#' Extract Forced migrationn only
mig$Forcibly_displaced <-toupper(mig$Forcibly_displaced)
mig <- mig[mig$Forcibly_displaced =='YES', ]
mig <- na.omit(mig)

coordinates(mig) <- ~Long+Lat
proj4string(mig) <-  CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") 
mig <- st_as_sf(mig)
st_crs(mig) <- st_crs(grd)
plot(st_geometry(grd))
plot(st_geometry(region), add=T)
plot(st_geometry(mig), pch=16,col="red", cex=0.5, add=T)

#' Count number of migration points in a grid


df <- sapply(st_intersects(grd, mig), 
             function(i){if(length(i)>0){ st_drop_geometry( mig[i,]) %>% group_by(YEAR) %>% tally() }else{NA} }
             )

for(id in 1:length(df)){
  if(length(df[[id]]) > 1){df[[id]] <- data.frame(id = id, df[[id]]) }
  }
df_m <- do.call(rbind, df)
df_m <- na.omit(df_m)
names(df_m)[names(df_m)=='n'] <-'Migration'


#' =======================================================================================
#' 2.0 CONFLICT PRESENCE in GRIDS (MEGA-PIXELS)
#' =======================================================================================

#Conflict data
filename <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/_global/conflict/Africa.xlsx"
conf <- readxl::read_excel(filename, sheet = 1)
country <- c("Democratic Republic of Congo", 'Uganda', 'Ethiopia', 'South Sudan', 'Tanzania', 'Kenya',
               'Rwanda', 'Burundi', 'Somalia', 'Djibouti', 'Eritrea')
conf <- conf[conf$YEAR >= 2018 & conf$COUNTRY==country,]

coordinates(conf) <- ~LONGITUDE+LATITUDE
proj4string(conf) <-  CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") 
conf <- st_as_sf(conf)
st_crs(conf) <- st_crs(grd)

#' Count number of migration points in a grid

temp <- sapply(st_intersects(grd, conf), 
             function(i){if(length(i)>0){ st_drop_geometry(conf[i,]) %>% group_by(YEAR) %>% tally() }else{NA} }
)

for(id in 1:length(temp)){
  if(length(temp[[id]]) > 1){temp[[id]] <- data.frame(id = id, temp[[id]]) }
}
df_c <- do.call(rbind, temp)
df_c <- na.omit(df_c)
names(df_c)[names(df_c)=='n'] <-'Conflict'

#' Merge Conflict & Migration into Megapixels
dff <- merge(df_m, df_c, by =c("id","YEAR"))

#' Merge @dff with megapixels

dff <- merge(grd, dff, by ='id')







#' OLD CODE xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#' 
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



EA_region <- terra::crop(EA_region, region)
values(EA_region) <- 1:ncell(EA_region)
EA_region <- terra::mask(EA_region, region)


# Define function to read and process geocoded data
read_and_process_geocode <- function(year) {
  file_path <- paste0(root,"intermediate/fms_geocoded/geo_", year, ".csv")
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