#############################################################
# script to create all necesary folders
# change ISO variable 
#



source_folders <- data.frame(folders =  list.dirs("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/KEN/", recursive = T))
source_folders <- base::subset(source_folders,!grepl("old|season_type",source_folders$folders) )
source_folders$folders <- gsub("KEN", "iso3", source_folders$folders) 

ISO <- "GTM"

for(i in source_folders$folders){
  i <- gsub("iso3", ISO, i)
  if(!file.exists(i)){dir.create(i, recursive = T)}
}