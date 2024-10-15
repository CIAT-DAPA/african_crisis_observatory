# ***********
# From monthly SPEI rasters, extract values at coordinates indicated in a specified csv file
# at quartely range i.e., every three months and compute the average
# store the output values in a csv file
# ***********
# to do
# read the folder, list the tif files using full paths and create a raster stack of all tifs
# extract spei at csv coordinates

rm(list=ls(all=TRUE))
library(terra)
library(sf)

monthly_spei <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/intermediate/spei/monthly_spei"
#coords_csv <- "C:/Users/bchepngetich/Downloads/KEN_2024_coordinates.csv"

kenya <- vect("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/raw/HH/KEN_2024_coordinates.shp")

# read csv containing coordinates list
#coords <- read.csv(coords_csv)
#coords <- coords[, c("longnum","latnum")]
#spei_coords <- coords
# list spei raster files 
spei <- list.files(monthly_spei, pattern = "\\.tif$", full.names = TRUE)

years <- 2001:2023 #c("2003", "2008", "2009", "2014", "2022")
years <- as.character(years)
#year_files <- grep(years, spei, value = TRUE)
selectFiles <- function(years){
  grep(years, spei, value = TRUE)
}
files <- lapply(years, selectFiles)
img <- rast(unlist(files))
df <- terra::extract(img, kenya, bind=T)
df <- as.data.frame(df)
temp <- df[, -4]

dff <- reshape(temp,
              direction = "long",
              varying = list(names(temp)[-c(1:5)]),
              v.names = "Spei",
              idvar = c("hv001", "year", "quarter", "lat", "long"),
              timevar = "Year_Month",
              times= names(temp)[-c(1:5)],
              )
lr <- reshape(LR, direction="long", varying = 2:ncol(LR), v.names="NDVI", timevar="year", times=2000:2015)

aa=melt(df, id.vars=c("hv001", "year", "quarter", "FID", "lat", "long"), measure.vars="SPEI")
# for (year in years) {
#   yr <- as.character(year)
#   print(yr)
#   year_files <- grep(yr, spei, value = TRUE)
#   print(year_files)
#   year_files <- lapply(year_files, rast)
#   print("rast done")
#   year_files <- do.call(c, year_files)
#   year_files <- terra::mask(year_files, kenya)
#   print("done mask")
#   x <- 1
#   #for (i in seq(1,nlyr(year_files), by=3)){
#   for (i in seq(1,nlyr(year_files), by=3)){
#     print(x)
#     quarter <- year_files[[i:(i+2)]]
#     quarter <- app(quarter, mean)
#     names(quarter) <- paste0(yr,"-",x)
#     spei_values <- extract(quarter, coords)
#     #spei_coords <- cbind(spei_coords, spei_values)
#     x <- x + 1
#     print("done")
#   }
# }


#id_cols <- grep("ID", colnames(spei_coords))
#final <- spei_coords[,-id_cols]
write.csv(df, "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/raw/HH/SPEI_Values.csv")

