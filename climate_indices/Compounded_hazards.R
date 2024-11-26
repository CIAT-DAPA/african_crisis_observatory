#********
#*This script combines NWLD, heat stress on humans, livestock, SPI
#*Author:Brenda Chepngetich, 2024
#*******
rm(list=ls(all=TRUE))
library(terra)

wd <- "C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Brenda/Hazards/"
#Baseline hazards
NWLD <- terra::rast(paste0(wd,"Floods/Floods/Karamoja_NWLD.tif"))
SPI <- terra::rast("C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/drought_spi/gee_spi/baseline_spi_raster_average.tif")
SPI_ <- terra::resample(SPI, NWLD)
HSH <- terra::rast(paste0(wd,"Heat stress/mean_Human_HSI.tif"))
HSH_ <- terra::resample(HSH, NWLD)
HSC <- terra::rast(paste0(wd,"Heat stress/HSC.tif"))
HSC_ <- terra::resample(HSC, NWLD)

normalization <- function(x){
  y <- (x - min(x, na.rm=TRUE))/(max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  return(y)
}
NWLD_ <- (terra::app(NWLD, normalization)) * 100
HSC_ <- (terra::app(HSC_, normalization)) * 100
HSH_ <- (terra::app(HSH_, normalization)) * 100
SPI_ <- (terra::app(SPI_, normalization)) * 100
SPI_Inverted <- 100 - SPI_

combined <- (NWLD_ + HSC_ + HSH_ + SPI_Inverted) / 4
terra::writeRaster(combined, filename = paste0(wd,"compounded_hazards.tif"))

#overall suitability
hazards <- terra::rast(paste0(wd,"compounded_hazards.tif"))
hazards <- (terra::app(hazards, normalization) * 100)
agnes <- "C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Projects/AGNES/3_Data/3_Outputs/"
locust <- terra::rast(paste0(agnes,"Locust_suitability_Index/Outputs/risk_maps/tif/final_baseline__locust_risk.tif"))
locust <- terra::resample(locust, hazards)
locust <- (terra::app(locust, normalization) * 100)
species <- terra::rast("C:/Users/bchepngetich/OneDrive - CGIAR/SA_Team/Brenda/Invasive_Species/data/karamoja_occurence.tif")
species <- terra::resample(species, hazards)
species <- (terra::app(species,normalization) * 100)
pasture <- terra::rast(paste0(agnes,"suitability_1991_2020/models/average_suitability_baseline.tif"))
pasture <- terra::resample(pasture, hazards)
pasture <- (terra::app(pasture, normalization)) * 100
pasture  <- 100 - pasture

overall <- (hazards + locust + species + pasture ) / 4
terra::writeRaster(overall,filename = paste0(wd,"overall_suitability.tif"))
