#Author:Brenda Chepngetich

wd <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/Data/CMIP6/monthly/"
setwd(wd)
first <- paste0(wd,"CMIP6_ACCESS-ESM1-5_ssp585_r1i1p1f1_pr_Africa_monthly.tif")
x <- terra::rast(first)
x
names(x)[24]
pr_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", 'pr', ".*"))

tasmax_files
pr <- terra::rast(paste0(wd,"CMIP6_ACCESS-ESM1-5_ssp585_r1i1p1f1_pr_Africa_monthly.tif"))
pr
temp <- terra::rast(paste0(wd,"CMIP6_ACCESS-ESM1-5_ssp585_r1i1p1f1_tasmax_Africa_monthly.tif"))
temp
jan_rasters <- temp[[seq(1, nlyr(temp), 12)]]
jan_rasters
mean_jan <- mean(jan_rasters)
mean_jan
plot(mean_jan)
#get tas max monthly means
tasmax_files <- list.files(pattern=paste0(".*", 'ssp585', ".*", "tasmax", ".*"))