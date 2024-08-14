#download monthly aggregated CMIP6 data from UDG mirror, using UDG tools
# August 2024

#java headspace
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
#install required packages
install.packages("rgdal")
install.packages("rJava")
install.packages('rgdal',repos="http://www.stats.ox.ac.uk/pub/RWin")
#load libraries
library(rJava)
library(devtools)
library(loadeR)
library(loadeR.java)
library(climate4R.UDG) 
library(transformeR)
library(downscaleR)
library(visualizeR)
library(rgdal)
library(tidyverse)
library(terra)

remove.packages("climate4R.UDG")
install_github("SantanderMetGroup/loadeR")
# # in case the above packages do not exist, install as follows
install_github(c("SantanderMetGroup/loadeR.java",
                 "SantanderMetGroup/climate4R.UDG",
                 "SantanderMetGroup/loadeR",
                 "SantanderMetGroup/transformeR",
                 "SantanderMetGroup/visualizeR",
                 "SantanderMetGroup/convertR",
                 "SantanderMetGroup/climate4R.indices",
                 "SantanderMetGroup/downscaleR"))


#working directory
wd <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/Data/CMIP6/Monthly_1991_2050"
wd
if (!file.exists(wd)) {dir.create(wd, recursive=TRUE)}

lons <- c(-23, 59)  # Africa
lats <- c(-37, 40)   # Africa
years.hist <- 1991:2014
years.rcp <- 2015:2050
#varname one of "tas","tasmin","tasmax","pr"
#rcp one of "ssp126", "ssp245", "ssp370", "ssp585"

#GCMs of interest: ACCESS-ESM1-5, EC-Earth3-Veg, INM-CM5-0, MPI-ESM1-2-HR, MRI-ESM2-0
dataset_list <- c("CMIP6_ACCESS-ESM1-5_scenario_r1i1p1f1",
                  "CMIP6_MPI-ESM1-2-HR_scenario_r1i1p1f1",
                  "CMIP6_EC-Earth3_scenario_r1i1p1f1",
                  "CMIP6_INM-CM5-0_scenario_r1i1p1f1",
                  "CMIP6_MRI-ESM2-0_scenario_r1i1p1f1")

#function to download CMIP6 data
downloadCMIP6 <- function(ds_name="CMIP6_ACCESS-ESM1-5_scenario_r1i1p1f1", rcp="ssp585", varname="pr", 
                          years.hist=1991:2014, years.rcp=2015:2050, lons=c(-23, 59), lats=c(-37, 40),
                          basedir) {
  #info
  cat("dataset=", ds_name, "/ rcp=", rcp, "/ variable=", varname, "\n")
  
  #names of datasets per CMIP6 Atlas (see https://github.com/SantanderMetGroup/ATLAS/)
  dataset.hist <- gsub(pattern="_scenario_", replacement="_historical_", x=ds_name)
  dataset.rcp <- gsub(pattern="_scenario_", replacement=paste0("_", rcp, "_"), x=ds_name)
  
  #function to load data
  load.data <- function (dset, years, var) loadGridData(dataset = dset, var = var,
                                                        years = years,
                                                        latLim = lats, lonLim = lons,
                                                        season = 1:12,aggr.m = "mean") 
  #file name, historical
  fname_his <- paste0(basedir, "/", dataset.hist,"_",varname,"_Africa_monthly.tif")
  if (!file.exists(fname_his)) {
    #loading mean temperature, historical
    cat("downloading historical data, please wait...\n")
    data_his <- load.data(dataset.hist, years.hist, var=varname)
    cat("downloaded")
    #convert 'grid' to sp object
    r_his <- grid2sp(data_his)
    r_his <- terra::rast(r_his)
    
    #write file
    terra::writeRaster(r_his, filename=fname_his)
    
    #clean up
    rm(data_rcp)
    gc(verbose=FALSE, full=TRUE)
  } else {
    cat("historical data already exists, loading\n")
    r_his <- terra::rast(fname_his)
  }
  
  #rcp data, file name
  fname_rcp <- paste0(basedir, "/", dataset.rcp,"_",varname,"_Africa_monthly.tif")
  if (!file.exists(fname_rcp)) {
    #loading mean temperature, rcp
    cat("downloading rcp data, please wait...\n")
    data_rcp <- load.data(dataset.rcp, years.rcp, var=varname)
    cat("downloaded rcp")
    #convert 'grid' to sp object
    r_rcp <- grid2sp(data_rcp)
    r_rcp <- terra::rast(r_rcp)
    
    #write file
    terra::writeRaster(r_rcp, filename=fname_rcp)
    
    #clean up
    rm(data_rcp)
    gc(verbose=FALSE, full=TRUE)
  } else {
    cat("rcp data already exists, loading\n")
    r_rcp <- terra::rast(fname_rcp)
  }
  return(list(his=r_his, rcp=r_rcp))
}

#run function
for (i in 1:length(dataset_list)) {
  for (scenario in c("ssp585")) {
    for (varname in c("tasmin", "tasmax", "pr")) {
      cmip6_data <- downloadCMIP6(ds_name=dataset_list[i], 
                                  rcp=scenario,  
                                  varname=varname, 
                                  years.hist=1991:2014, 
                                  years.rcp=2015:2050, 
                                  lons=c(-23, 59), 
                                  lats=c(-37, 40), 
                                  basedir=wd)
      rm(cmip6_data)
      gc(verbose=FALSE, full=TRUE)
      .jcall("java/lang/System", method = "gc")
    }
  }
}


