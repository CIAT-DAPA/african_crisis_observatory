# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Identify socio-economic indicators based on impact pathways
# Author: Harold Achicanoy, Andres Mendez
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# R options
g <- gc(reset = T); rm(list = ls()) # Emptying the garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))   # Loading R-packages
suppressMessages(pacman::p_load(tidyverse,readxl))

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'
iso  <- "KEN"
cntr <- 'Kenya'

select_eco_vars <- function(root, iso, cntr, ip){
  
  vr_tbl <- readxl::read_excel(path = paste0(root,'/Hostpots_data_dictionary.xlsx'), sheet = 1)
  # Impact pathways (IP) table
  ip_tbl <- readxl::read_excel(path = paste0(root,'/Africa Climate Security_Country Pathways.xlsx'), sheet = 2)
  ip_tbl <- ip_tbl %>% dplyr::filter(Country == cntr & IP_id == ip)
  
  ## ip_tbl$Values %>% purrr::map(.f = function(txt){base::strsplit(x = txt, split = "\\s{2,}")[[1]]})
  
  # Key climate arguments from IP
  key <- ip_tbl$Values[setdiff(1:nrow(ip_tbl),grep(pattern = '[cC]limate', x = ip_tbl$Dimension))]
  key <- key[!is.na(key)]
  key <- strsplit(x = key, split = "\\s{2,}") %>% unlist() %>% unique() %>% na.omit() %>% as.character() %>% tolower()
  
  region_key <- unique(ip_tbl$Region_key)
  region_val <- unique(ip_tbl$Region_value)
  
  # Climate variables from dictionary
  val <- vr_tbl$Category[vr_tbl$Classification != 'Climate']
  val <- strsplit(x = val, split = ';')
  
  # Function to do the match between key and values
  select_vars <- function(reference = key, values = val[[6]]){
    if(all(is.na(values))){ 
      out <- NA; return(NA)
    }else{
      grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')
      out <- unlist(grep2(pattern = values, x = key))
      if(length(out) > 0){out <- 1} else {out <- NA}
      return(out)
    }
  }
  
  selected <- val %>% 
    purrr::map(.f = function(x){select_vars(reference = key, values = x)}) %>% 
    unlist()
  
  result <- data.frame(Variable     = vr_tbl$Variable[vr_tbl$Classification != 'Climate'],
                       Code         = vr_tbl$Code[vr_tbl$Classification != 'Climate'],
                       Selected     = selected,
                       Region_key   = region_key,
                       Region_value = region_val)
  result <- result %>% tidyr::drop_na()
  result <- dplyr::left_join(x = result, y = vr_tbl %>% dplyr::select(Code,Threshold,Percentile,Classification), by = 'Code')
  result <- result[!is.na(result$Percentile),]
  thrs <- grep(pattern = ';', x = result$Percentile)
  if(length(thrs) > 0){
    ustckd <- purrr::map(.x = thrs, .f = function(th){
      prcs <- strsplit(x = result[th,'Percentile'], split = ';')[[1]]
      rptd <- result[rep(th, times = length(prcs)),]
      rptd$Percentile <- prcs
      return(rptd)
    }) %>%
      dplyr::bind_rows()
    result <- result[-thrs,]
    result <- rbind(result, ustckd); rownames(result) <- 1:nrow(result); rm(ustckd)
    result$Percentile <- as.numeric(result$Percentile)
  } else {
    result$Percentile <- as.numeric(result$Percentile)
  }
  
  return(result)
  
}
dest_dir <- paste0(root ,"/data/",iso,"/_results/hotspots/soc_eco_all_variables.csv")
if(!file.exists(dest_dir)){
  ip_tbl <- readxl::read_excel(path = paste0(root,'/Africa Climate Security_Country Pathways.xlsx'), sheet = 2)
  ip_tbl <- ip_tbl %>% dplyr::filter(Country == cntr)
  ips <- unique(ip_tbl$IP_id); rm(ip_tbl)
  ips %>%
    purrr::map(.f = function(ip){
      vrs <- select_eco_vars(root, iso, cntr, ip = ip)
      vrs$IP_id <- ip
      return(vrs)
    }) %>%
    dplyr::bind_rows() %>%
    write.csv(x = ., file = dest_dir, row.names = F)
} else {
  read.csv(file = dest_dir)
}

# View(result)
# Execute climate clusters code ...
