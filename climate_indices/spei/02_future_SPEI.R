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
compute_pet <- function(tmin, tmax) {
  tmean <- (tmin + tmax) / 2
  pet <- 0.0023 * (tmean + 17.8) * (tmax - tmin) ^ 0.5 * 0.408  # Hargreaves equation
  return(pet)
}

# Compute future PET
pet <- compute_pet(tmin, tmax)
pet <- terra::app(tmin, tmax, fun = compute_pet)

# Compute future water balance (Precipitation - PET)
water_balance_future <- ppt - pet

#output_path <- "path/to/output/water_balance_future.tif"
#writeRaster(water_balance_future, output_path, overwrite = TRUE)

calc_spei <- function(x){
  x <- ts(data = x, start = c(1979,1), frequency = 12)
  Spei <- SPEI::spei(data = x, scale = 6, na.rm = T)$fitted |> as.numeric()
  return(Spei)
}

# Calculate SPEI
# Spei <- terra::app(x = blc, fun = calc_spei, cores = 20)
Spei <- terra::app(x = water_balance_future, fun = function(i, ff) ff(i), cores = 20, ff = calc_spei)
names(Spei)<-names(water_balance_future)
