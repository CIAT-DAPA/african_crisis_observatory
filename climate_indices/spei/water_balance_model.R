
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')


iso <- 'KEN'
soil_cp_pth =  paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/temp/soilcp.tif")
soil_sat_pth = paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/temp/soilsat.tif")
# Raster template
tmp <- terra::rast('//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif')

#Loading data
soil_cp_pth =  paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/temp/soilcp.tif")
soil_sat_pth = paste0("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/", iso, "/climatic_indexes/temp/soilsat.tif")
rain <- ppt <- terra::rast('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/Data/CMIP6/monthly/CMIP6_EC-Earth3_ssp585_r1i1p1f1_pr_Africa_monthly.tif')
etmax <- terra::rast('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/intermediate/spei/historical_mean_monthly_PET.tif')
rain <- rain[[c(49:nlyr(rain))]]
sat <- terra::rast(soil_sat_pth)
scp <- terra::rast(soil_cp_pth)

shp <- geodata::gadm(country = 'Kenya', level = 0, path = tempdir())

#Cropping the data
rain <- terra::crop(rain, shp) %>% terra:: mask(., shp)
sat <- terra::crop(sat, shp) %>% terra:: mask(., shp)%>% terra::resample(.,rain, 'bilinear')
scp <- terra::crop(scp, shp) %>% terra:: mask(., shp)%>% terra::resample(.,rain, 'bilinear')
etmax <- terra::crop(etmax, shp) %>% terra:: mask(., shp)%>% terra::resample(.,rain, 'bilinear')
etmax <- rast(replicate(nlyr(rain), etmax, simplify = FALSE))
tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp) %>% terra::resample(.,rain, 'bilinear')
tmp[!is.na(tmp)] <- 1
# Compute water balance model
AVAIL <<- tmp
AVAIL[!is.na(AVAIL)] <- 0

water_bal <- eabyep_calc(soilcp  = scp,
                         soilsat = sat,
                         avail   = AVAIL,
                         rain    = rain,
                         evap    = etmax)

ERATIO  <- watbal %>% purrr::map('Eratio') %>% terra::rast()
LOGGING <- watbal %>% purrr::map('Logging') %>% terra::rast()
IRR     <- etmax - rain
GDAY    <- terra::lapp(x = terra::sds(tav, ERATIO), fun = function(TAV, ERATIO){ifelse(TAV >= 6 & ERATIO >= 0.35, 1, 0)})

NDWS    <- terra::app(x = ERATIO, fun = function(ERATIO){ifelse(ERATIO < 0.5, 1, 0)}) %>% sum()
NWLD    <- terra::app(x = LOGGING, fun = function(LOGGING){ifelse(LOGGING > 0, 1, 0)}) %>% sum()
NWLD50  <- sum(LOGGING > (sat*0.5))
NWLD90  <- sum(LOGGING >= sat)
IRR     <- sum(IRR)

