#*TODO: Compute Future SPEI based on CMIP6 SSP585 Tmax, TMin, Pr data 
#*
#*Author; Victor Korir

################################################################################
library(terra)


#Load projected Tmax, Tmin and Precipitation sourced from CMIP6
ppt <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/1.Data/CMIP6/intermediate_africa/monthly_annual/CMIP6_ACCESS-ESM1-5_ssp245_r1i1p1f1_pr_Africa_monthly_annual_2030_2050.tif')
tmin <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/1.Data/CMIP6/intermediate_africa/monthly_annual/CMIP6_ACCESS-ESM1-5_ssp245_r1i1p1f1_tasmin_Africa_monthly_annual_2030_2050.tif')
tmax <- terra::rast('//CATALOGUE/WFP_ClimateRiskPr1/1.Data/CMIP6/intermediate_africa/monthly_annual/CMIP6_ACCESS-ESM1-5_ssp245_r1i1p1f1_tasmax_Africa_monthly_annual_2030_2050.tif')

#Future PET using Hargreaves equation
# compute_pet <- function(tmin, tmax) {
#   tmean <- (tmin + tmax) / 2
#   pet <- 0.0023 * (tmean + 17.8) * (tmax - tmin) ^ 0.5 * 0.408  # Hargreaves equation
#   return(pet)
# }

# # Compute future PET
# pet <- compute_pet(tmin, tmax)
# pet <- terra::app(tmin, tmax, fun = compute_pet)

# Read historical mean of PET per month and replicate for the number of years
KEN_shp <- geodata::gadm(country = 'Kenya', level = 0, path = tempdir())
years <- 21
historical_PET <- rast('//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/intermediate/spei/historical_mean_monthly_PET.tif')
historical_PET_replic <- rast(replicate(years, historical_PET, simplify = FALSE))
historical_PET_replic <- terra::crop(historical_PET_replic, KEN_shp)
historical_PET_replic <- terra::mask(historical_PET_replic, KEN_shp)
ppt <- terra::crop(ppt, KEN_shp)
ppt <- terra::mask(ppt, KEN_shp)

#Resampling PET
historical_PET_replic <-resample(historical_PET_replic, ppt, method='bilinear')
# Compute future water balance (Precipitation - PET)
water_balance_future <- ppt - historical_PET_replic

#output_path <- "path/to/output/water_balance_future.tif"
#writeRaster(water_balance_future, output_path, overwrite = TRUE)

calc_spei <- function(x){
  x <- ts(data = x, start = c(2030,1), frequency = 12)
  Spei <- SPEI::spei(data = x, scale = 12, na.rm = T)$fitted |> as.numeric()
  return(Spei)
}

# Calculate SPEI
Spei <- terra::app(x = water_balance_future, fun = calc_spei, cores = 20)
#Spei <- terra::app(x = water_balance_future, fun = function(i, ff) ff(i), cores = 20, ff = calc_spei)
names(Spei)<-names(water_balance_future)
plot(mean(Spei, na.rm=T))
