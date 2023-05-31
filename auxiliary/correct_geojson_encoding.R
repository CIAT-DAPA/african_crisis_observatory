pacman::p_load(sf, tidyverse, stringi)


gj <- sf::st_read("D:/OneDrive - CGIAR/Documents/prueba_1/geojsons/SEN_megapixels.geojson")


gj_original <- sf::st_read("Z:/1.Data/Palmira/CSO/data/SEN/_results/clim_conflict_ips_overlays.geojson",  options = "ENCODING=WINDOWS-1252")


gj$NAME_1 <- iconv(gj_original$NAME_1, from = "UTF-8", to = "MS-ANSI")
gj$NAME_2 <- iconv(gj_original$NAME_2, from = "UTF-8", to = "MS-ANSI")
gj$NAME_3 <- iconv(gj_original$NAME_3, from = "UTF-8", to = "MS-ANSI")


st_write(gj, "D:/OneDrive - CGIAR/Documents/prueba_1/geojsons/SEN_megapixels.geojson", delete_dsn = T)

#####---------
gj <- sf::st_read("D:/OneDrive - CGIAR/Documents/prueba_1/geojsons/ETH_megapixels.geojson")

gj %>% 
  filter(!is.na(clim_cluster_short_label)) %>% 
  st_write(., "D:/OneDrive - CGIAR/Documents/prueba_1/geojsons/ETH_megapixels.geojson", delete_dsn = T)


