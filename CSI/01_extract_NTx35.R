#' Extract CV Rainfall per admin level 1
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
#' @season season type
#'@data_path directory path for chirps data
#'@variab Variable being extracted

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSI/' 
data_path <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/'
iso <- 'SOM'
season <- 'season_type_1'
variab <- 'NTx35'
img <- terra::rast(paste0(data_path,iso, '/climatic_indexes/temp/season_type_1/',variab,'.tif'))

bdy1 <- gadm(country=iso, level=1, path=paste0(root,'_global/countries/'))
knitr::kable(head(bdy1[['NAME_1']]))

tic("Extract average number days where temp was above 30 degress rainfall..")
temp <- exactextractr::exact_extract(x = img, y = sf::st_as_sf(bdy1), 'mean', full_colnames = T, append_cols="NAME_1") 
exectime <- toc()
exectime <- (exectime$toc - exectime$tic)/3600
exectime

colnames(temp)[-1] <- gsub("mean.","",names(temp)[-1])
df <- reshape2::melt(temp, variable.name = "date", value.name = paste0(variab,"_",season),id.vars = "NAME_1")
knitr::kable(head(df, n=5))

filename <- paste0(root,variab,"/",iso,"/",iso,"_", variab,"_",season,".csv")
dir.create(paste0(root,variab,"/",iso), recursive = T)
write.csv(df, filename)
