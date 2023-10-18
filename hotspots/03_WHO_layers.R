#' #################################
#' @author alex orestein
#' 2023
#' Script to get WHO layers for CSO
#' see https://github.com/oren-sa/cso_contextual
require(pacman)
pacman::p_load(raster, rgdal, sf, fasterize, dplyr, tidyr, tidyverse)

main_dir <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data"

#prep- create folders for "whole_country_rasters", "masked_rasters" and add "megapixels.geojson" to your desired output directory.
# This script should be run first, followed by the "megapixels_zonal" script afterwards.
# FEWS and LHZ data are not included in this script. They were processed separately through QGIS for inclusion into megapixels.
#  Megapixels.geojson must contain the megapixels you will use for analysis.
# Once you specify the directory paths, you can re-run this code, only changing the 3 letter country code.


COUNTRY<- "NER" #be sure to replace this with the 3 letter country code

setwd(paste0(main_dir, "/", COUNTRY))

megapixels<- sf::st_read("_results/clim_conflict_ips_overlays.geojson")
megapixels$id<- 1:nrow(megapixels) # format megapixels

#create AOI and a population mask- this will pull directly from the data dictionary
aoi<-raster("population_density/medn_popd.tif")
aoi_mask<- aoi>1 
x11();plot(aoi) ##run this if you want to test what the AOI data looks like
  
megapixels_raster<- fasterize((megapixels),aoi,field="id")
megapixels_mask <- function(x){x*(megapixels_raster>0)}


  #create a function to resample the data to match the population mask and remove data from unpopulated pixels
  resample_aoi <- function(x){ (aoi_mask*(raster::resample(x,aoi,method="bilinear")))*(megapixels_raster>0)   }

#Get Data dictionary layers and resample them
popd<- aoi
stunting<- resample_aoi(raster("child_growth_failure/medn_stunting.tif"))
wasting<- resample_aoi(raster("child_growth_failure/medn_wasting.tif"))
underweight<- resample_aoi(raster("child_growth_failure/medn_underweight.tif"))
female_ed<-resample_aoi(raster("education/medn_female_edu.tif"))
male_ed<- resample_aoi(raster("education/medn_male_edu.tif"))
diff_ed<- resample_aoi(raster("education/medn_difference_edu.tif"))
nightlights <- resample_aoi(raster("nightlights/medn_lght.tif"))
piped_water<- resample_aoi(raster("sanitation/medn_piped_water.tif"))
sanitation <- resample_aoi (raster("sanitation/medn_sanitation_facilities.tif"))
migration <- resample_aoi(raster("migration/rcnt_migration.tif"))
awe <- resample_aoi(raster(paste0("wealth_index/",COUNTRY,"_AWE.tif")))
rwi<- resample_aoi(raster(paste0("wealth_index/",COUNTRY,"_rwi.tif")))

# aa <- st_read("D:/OneDrive - CGIAR/Desktop/ZWE_mpx.geojson")
# 
# #add non data dictionary layers here. these are outside of the data dictionary so you will need to specify their paths below
dep_ratio <- resample_aoi(raster(paste0(main_dir, "/_global/dependency_ratio/AFR_2010_SubNat_DepRatio.tif"))) 
lhz <- resample_aoi(raster(paste0(main_dir, "/_global/livelihood_africa/LHZ_AFRICA.tif")))
# 
# #Write all data into rasters masked by the megapixels
# writeRaster(stack(masked_list), paste0("masked_rasters/",names(masked_list)), bylayer = TRUE, format='GTiff')
# 
# #Write all data into rasters without the mask
# writeRaster(stack(layers_list), paste0("whole_country_rasters/",names(layers_list)), bylayer = TRUE, format='GTiff')


################################################################
###############################################################
#######              SECOND PART              ################
#############################################################


#calculate population for each megapixels and add them into megapixels
pop<- as.data.frame(zonal (popd,megapixels_raster,fun='sum'))
megapixels<-merge(megapixels,pop, by.x = 'id', by.y= 'zone')
megapixels <- rename(megapixels, pop = sum)
megapixels$pop<- round(megapixels$pop, digits = 0)


#create functions to weigh megapixel data by population
weight_zonal<- function(x){as.data.frame(zonal(x*popd,megapixels_raster,fun = 'sum'))} 
merge_megapixel<- function(y){merge(megapixels,y, by.x = 'id', by.y= 'zone')} 
value_megapixel<- function(x){round(x/megapixels$pop,digits=2)}  

##These lines apply the weights to each of the raster layers. If you wish to add more raster layers, copy and paste the code below accordingly.

zonal_female_ed<-weight_zonal(female_ed)
colnames(zonal_female_ed)[2] <- "female_ed_weighted"
megapixels<-merge_megapixel(zonal_female_ed)
megapixels$female_ed <- value_megapixel(megapixels$female_ed_weighted)
megapixels<- subset(megapixels, select = -(female_ed_weighted))

zonal_male_ed<-weight_zonal(male_ed)
colnames(zonal_male_ed)[2] <- "male_ed_weighted"
megapixels<-merge_megapixel(zonal_male_ed)
megapixels$male_ed <- value_megapixel(megapixels$male_ed_weighted)
megapixels<- subset(megapixels, select = -(male_ed_weighted))

zonal_stunting<-weight_zonal(stunting)
colnames(zonal_stunting)[2] <- "stunting_weighted"
megapixels<-merge_megapixel(zonal_stunting)
megapixels$stunting <- value_megapixel(megapixels$stunting_weighted)
megapixels<- subset(megapixels, select = -(stunting_weighted))

zonal_wasting<-weight_zonal(wasting)
colnames(zonal_wasting)[2] <- "wasting_weighted"
megapixels<-merge_megapixel(zonal_wasting)
megapixels$wasting <- value_megapixel(megapixels$wasting_weighted)
megapixels<- subset(megapixels, select = -(wasting_weighted))

zonal_underweight<-weight_zonal(underweight)
colnames(zonal_underweight)[2] <- "underweight_weighted"
megapixels<-merge_megapixel(zonal_underweight)
megapixels$underweight <- value_megapixel(megapixels$underweight_weighted)
megapixels<- subset(megapixels, select = -(underweight_weighted))

zonal_difference_ed<-weight_zonal(diff_ed)
colnames(zonal_difference_ed)[2] <- "difference_ed_weighted"
megapixels<-merge_megapixel(zonal_difference_ed)
megapixels$difference_ed <- value_megapixel(megapixels$difference_ed_weighted)
megapixels<- subset(megapixels, select = -(difference_ed_weighted))

zonal_nightlights<-weight_zonal(nightlights)
colnames(zonal_nightlights)[2] <- "nightlights_weighted"
megapixels<-merge_megapixel(zonal_nightlights)
megapixels$nightlights <- value_megapixel(megapixels$nightlights_weighted)
megapixels<- subset(megapixels, select = -(nightlights_weighted))

zonal_piped_water<-weight_zonal(piped_water)
colnames(zonal_piped_water)[2] <- "piped_water_weighted"
megapixels<-merge_megapixel(zonal_piped_water)
megapixels$piped_water <- value_megapixel(megapixels$piped_water_weighted)
megapixels<- subset(megapixels, select = -(piped_water_weighted))

zonal_sanitation<-weight_zonal(sanitation)
colnames(zonal_sanitation)[2] <- "sanitation_weighted"
megapixels<-merge_megapixel(zonal_sanitation)
megapixels$sanitation <- value_megapixel(megapixels$sanitation_weighted)
megapixels<- subset(megapixels, select = -(sanitation_weighted))

zonal_migration<-weight_zonal(migration)
colnames(zonal_migration)[2] <- "migration_weighted"
megapixels<-merge_megapixel(zonal_migration)
megapixels$migration <- value_megapixel(megapixels$migration_weighted)
megapixels<- subset(megapixels, select = -(migration_weighted))

zonal_awe<-weight_zonal(awe)
colnames(zonal_awe)[2] <- "awe_weighted"
megapixels<-merge_megapixel(zonal_awe)
megapixels$awe <- value_megapixel(megapixels$awe_weighted)
megapixels<- subset(megapixels, select = -(awe_weighted))

zonal_rwi<-weight_zonal(rwi)
colnames(zonal_rwi)[2] <- "rwi_weighted"
megapixels<-merge_megapixel(zonal_rwi)
megapixels$rwi <- value_megapixel(megapixels$rwi_weighted)
megapixels<- subset(megapixels, select = -(rwi_weighted))

zonal_dep_ratio<-weight_zonal(dep_ratio)
colnames(zonal_dep_ratio)[2] <- "dep_ratio_weighted"
megapixels<-merge_megapixel(zonal_dep_ratio)
megapixels$dep_ratio <- value_megapixel(megapixels$dep_ratio_weighted)
megapixels<- subset(megapixels, select = -(dep_ratio_weighted))

#write the megapixels as geojson
## you can change this to whichever format by altering the extension in file name and choosing a new driver
#writeOGR(megapixels,paste0("/media/alex/LaCie/GIS/CSO/sample_data/",COUNTRY,"/",COUNTRY,"_megapixels1.geojson"),driver= "GeoJSON", layer='megapixels', overwrite_layer = TRUE)
sf::st_write(megapixels, paste0(main_dir, "/", COUNTRY, "/_results/clim_conflict_ips_overlays_who_layers.geojson"), delete_dsn  = T)
writexl::write_xlsx(sf::st_drop_geometry(megapixels), paste0(main_dir, "/", COUNTRY, "/_results/clim_conflict_ips_overlays_who_layers.xlsx"))
#calculate national averages- creates a function to calculate averages of all raster values.
nat_avgs<-data.frame(COUNTRY)
nat_avg<- function(x){cellStats((x*popd),sum)/cellStats(popd,sum)}

nat_avgs$popd<-nat_avg(popd)
nat_avgs$stunting<-nat_avg(stunting)
nat_avgs$wasting<-nat_avg(wasting)
nat_avgs$female_ed<-nat_avg(female_ed)
nat_avgs$male_ed<- nat_avg(male_ed)
nat_avgs$diff_ed<- nat_avg(diff_ed)
nat_avgs$nightlights<- nat_avg(nightlights)
nat_avgs$piped_water <- nat_avg(piped_water)
nat_avgs$sanitation <- nat_avg(sanitation)
nat_avgs$migration <- nat_avg(migration)
nat_avgs$awe<- nat_avg(awe)
nat_avgs$rwi<- nat_avg(rwi)
nat_avgs$dep_ratio <- nat_avg(dep_ratio)
nat_avgs$female_ed<-cellStats((female_ed*popd),sum)/cellStats(popd,sum)
nat_avgs$male_ed<-cellStats((male_ed*popd),sum)/cellStats(popd,sum)

write.csv(nat_avgs,paste0(main_dir, "/", COUNTRY, "/_results/national_average.csv"))








