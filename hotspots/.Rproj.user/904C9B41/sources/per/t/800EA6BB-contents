require(pacman)
pacman::p_load(raster, tidyverse, readxl, sf, sp )


ISO3 <- "UGA"

wealth_dir <- "D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/_global/wealth_index"

rwi_out_dir <-  paste0("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/", ISO3,"/wealth_index" )

if(!dir.exists(rwi_out_dir)){dir.create(rwi_out_dir)}

wealth_df <- read_csv(list.files(wealth_dir, pattern = paste0("^",ISO3), full.names = T))

mask <- raster("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/_global/masks/mask_world_1km.tif")


shp_country <- raster::shapefile(paste0("D:/CGIAR/Achicanoy Estrella, Harold Armando (Alliance Bioversity-CIAT) - African_Crisis_Observatory/data/", ISO3, "/_shps/",ISO3, ".shp"))


coordinates(wealth_df) <-  ~longitude+latitude
crs(wealth_df) <- "+proj=longlat +datum=WGS84 +no_defs"

r <- raster(resolution = 0.04166667, ext = extent(shp_country))
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
r <- raster::mask(r, shp_country)
r_f <- rasterize(wealth_df, r, fun = mean, field = "rwi")
r_f <- raster::resample(r_f, mask %>% raster::crop(., extent(shp_country)) )

writeRaster(r_f, paste0(rwi_out_dir, "/",ISO3, "_rwi.tif"), overwrite = T)

