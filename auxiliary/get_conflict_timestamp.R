require(tidyverse, terra)

af_conf <- readxl::read_excel("Z:/1.Data/Palmira/CSO/data/_global/conflict/Africa_1997-2022_Jul08.xlsx")
as_conf <- readxl::read_excel("Z:/1.Data/Palmira/CSO/data/_global/conflict/East-Asia-Pacific_2023-05-09.xlsx")
lt_conf <- readxl::read_excel("Z:/1.Data/Palmira/CSO/data/_global/conflict/LatinAmerica_2018-2022_Jul08.xlsx")

ISOS <- list.dirs("Z:/1.Data/Palmira/CSO/data", recursive = F) %>% 
  .[!grepl("_[a-zA-Z]+", .)] %>% 
  basename()

shp <- sf::st_read("Z:/1.Data/Palmira/CSO/data/_global/world_shapefile/all_country/all_countries.shp")


c_names <- sapply(ISOS, function(i){
 
   c_names <- shp %>% 
    dplyr::filter(!duplicated(ISO3)) %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(ISO3, NAME) %>% 
    dplyr::filter(ISO3 == i) %>% 
    pull(NAME)
  return(c_names)
  
}, simplify = T) %>% unlist




full_conf <- bind_rows(af_conf %>% select(EVENT_DATE, COUNTRY  )  ,
                       lt_conf %>% select(EVENT_DATE, COUNTRY  ),
                       as_conf %>% select(EVENT_DATE, COUNTRY  )) 

dts_conf_tbl <- lapply(c_names, function(cty){
  
  min_date <- full_conf %>% 
    filter(COUNTRY == cty) %>% 
    pull(EVENT_DATE) %>% 
    min()
  
  max_date <- full_conf %>% 
    filter(COUNTRY == cty) %>% 
    pull(EVENT_DATE) %>% 
    max()
  
  return(data.frame(start_date= min_date, end_date = max_date))
}) %>% bind_rows


ISOS <- ISOS[!grepl("SSD", ISOS)] 

dts_conf_tbl %>% 
  dplyr::mutate(country = c_names, ISO = ISOS) %>% 
  writexl::write_xlsx(., "Z:/1.Data/Palmira/CSO/data/country_conflic_timeframes.xlsx")















