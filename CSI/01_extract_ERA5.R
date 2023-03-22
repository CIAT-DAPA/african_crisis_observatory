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
iso <- 'UGA'
root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSI/'
data_path  <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr/1.Data/ERA5/'
variable <- '2m_temperature-24_hour_maximum'
bdy1 <- gadm(country=iso, level=1, path=paste0(root,'_global/countries/'))
knitr::kable(head(bdy1[['NAME_1']]))

files <- list.files(paste0(data_path,variable),pattern = (".nc$"), recursive = TRUE, full.names = TRUE)
knitr::kable(head(files), caption = 'CHIRPS raster files')
#df <- data.frame(NAME_1=bdy1[['NAME_1']])

#' Extract spatial CHIRPS rainfall by aggregating salinity in a given year in each admin 4 boundary.

tic("Load files..")
img <- terra::rast(files, win=ext(bdy1))
toc()

#x11()
plot(img[[1]])
plot(bdy1, add=T, type='l') 

names(img) <- gsub('.*AgERA5_|_final.*', '',files)
head(names(img))

terra::gdalCache(20000)
tic("Extract temperature..")
temp <- exactextractr::exact_extract(x = img, y = sf::st_as_sf(bdy1), 'mean', full_colnames = T, append_cols="NAME_1") 
exectime <- toc()
exectime <- (exectime$toc - exectime$tic)/3600
exectime

colnames(temp)[-1] <- gsub("mean.","",names(temp)[-1])
df <- reshape2::melt(temp, variable.name = "date", value.name = "ERA5_Tmax",id.vars = "NAME_1")
df$date <- as.Date(df$date, format = '%Y%m%d')
df$year <- format(df$date, format="%Y")
df$ERA5_Tmax <- df$ERA5_Tmax - 273.15
knitr::kable(head(df, n=5))
filename <- paste0(root,"ERA5/",iso,"/",iso,"_daily_ERA5_tmax_rainfall.csv")
dir.create(paste0(root,"ERA5/",iso), recursive = T)
write.csv(df, filename)
#' Aggregate rainfall per admin level 1 by year
#temp <- aggregate(df[,"chirps_rainfall"], df[,c("NAME_1","year"), drop=FALSE], mean, na.rm=T)
temp <- aggregate(ERA5_Tmax~NAME_1+year, data=df, mean, na.rm=T)
filename <- paste0(root,"ERA5/",iso,"/",iso,"_annual_ERA5_tmax_rainfall.csv")
write.csv(temp, filename)


