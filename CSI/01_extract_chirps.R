#' Extract Chirps Rainfall per admin level 1
#' By: Benson Kenduiywo
#' March 2023
#' R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
list.of.packages <- c('terra', 'geodata', 'doParallel', 'tictoc')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
#' Set up libraries
library(terra)
library(geodata)
library(tictoc)
#' variable declarations
#' @param iso is the country name ISO code.
#' @root root base directory for saving the results
#' @data_path directory path for chirps data
iso <- 'KEN'
root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSI/'
data_path  <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr/1.Data/Chirps/'

bdy1 <- gadm(country=iso, level=1, path=paste0(root,'_global/countries/'))
knitr::kable(head(bdy1[['NAME_1']]))

files <- list.files(paste0(data_path),pattern = (".tif$"), recursive = TRUE, full.names = TRUE)
knitr::kable(head(files), caption = 'CHIRPS raster files')
#df <- data.frame(NAME_1=bdy1[['NAME_1']])

#' Extract spatial CHIRPS rainfall by aggregating salinity in a given year in each admin 4 boundary.

tic("Load files..")
img <- terra::rast(files, win=ext(bdy1))
NAflag(img) <- -9999
toc()

#x11()
plot(img[[1]])
plot(bdy1, add=T, type='l') 

tic("Extract rainfall..")
temp <- exactextractr::exact_extract(x = img, y = sf::st_as_sf(bdy1), 'mean', full_colnames = T, append_cols="NAME_1") 
exectime <- toc()
exectime <- (exectime$toc - exectime$tic)/3600
exectime
colnames(temp)[-1] <- gsub("mean.chirps-v2.0.","",names(temp)[-1])
df <- reshape2::melt(temp, variable.name = "date", value.name = "chirps_rainfall",id.vars = "NAME_1")
df$date <- as.Date(df$date, format = '%Y.%m.%d')
df$year <- format(df$date, format="%Y")
knitr::kable(head(df, n=5))
filename <- paste0(root,"chirps/",iso,"/",iso,"_daily_chirps_rainfall.csv")
dir.create(paste0(root,"chirps/",iso))
write.csv(df, filename)
#' Aggregate rainfall per admin level 1 by year
#temp <- aggregate(df[,"chirps_rainfall"], df[,c("NAME_1","year"), drop=FALSE], mean, na.rm=T)
temp <- aggregate(chirps_rainfall~NAME_1+year, data=df, sum, na.rm=T)
filename <- paste0(root,"chirps/",iso,"/",iso,"_annual_chirps_rainfall.csv")
write.csv(temp, filename)


