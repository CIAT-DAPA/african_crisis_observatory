##Download MODIS fpar
#' ############################################################################################
#' Author: Benson Kenduiywo
#' ############################################################################################

#https://viirsland.gsfc.nasa.gov/Products/NASA/LAI_FparESDR.html
#https://lpdaac.usgs.gov/products/myd15a2hv006/
root <-  '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/raw/modis/'
dir.create(root, recursive = T)

list.of.packages <- c("luna","terra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)

prod <- getProducts("^MYD15A2H|^VNP15A2H|^VJ115A2H|")
product <- 'MYD15A2H'
#productInfo(product)
#login information

#saveRDS(df, paste0(root,'earthdata.rds'))
df <- readRDS(paste0(root,'earthdata.rds'))
iso  <- 'KEN'
start <- "2001-01-01"
end <- "2023-07-31"
aoi <- geodata::gadm(iso, level = 0, version = 'latest', path = tempdir())
path <- paste0(root, iso, '/')
dir.create(path, recursive = T)
mf <- luna::getNASA(product, start, end, aoi=aoi, download=TRUE,
                     path=path, username=df$user, password=df$pwd)

rlist <- purrr::map(mf, rast)
rlist <- terra::sprc(rlist)
im <- merge(rlist)
plot(im['Fpar_500m'])