# ----------------------------------------------------------------------------------- #
# Climate Security Observatory
# Identify climate indicators based on impact pathways
# Author: Harold Achicanoy, Andres Mendez
# Alliance Bioversity International - CIAT, 2022
# ----------------------------------------------------------------------------------- #

# # R options
# g <- gc(reset = T); rm(list = ls()) # Emptying the garbage collector
# .rs.restartR()                      # Restart R session
# options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
# suppressMessages(library(pacman))   # Loading R-packages
# suppressMessages(pacman::p_load(tidyverse,readxl))

# root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'
# iso <- "KEN"
# cntr <- 'Kenya'

select_clim_vars <- function(root, iso, cntr){
  # Data dictionary
  
  dest_dir <- paste0(root ,"/data/",iso,"/_results/cluster_results/climate/climate_selected_variables.csv")
  
  if(!file.exists(dest_dir)){
    
    vr_tbl <- readxl::read_excel(path = paste0(root,'/Hostpots_data_dictionary.xlsx'), sheet = 1)
    # Impact pathways (IP) table
    ip_tbl <- readxl::read_excel(path = paste0(root,'/Africa Climate Security_Country Pathways.xlsx'), sheet = 2)
    ip_tbl <- ip_tbl %>% dplyr::filter(Country == cntr)
    
    ## ip_tbl$Values %>% purrr::map(.f = function(txt){base::strsplit(x = txt, split = "\\s{2,}")[[1]]})
    
    # Key climate arguments from IP
    key <- ip_tbl$Values[grep(pattern = '[cC]limate', x = ip_tbl$Dimension)]
    key <- strsplit(x = key, split = "\\s{2,}") %>% unlist() %>% unique() %>% na.omit() %>% as.character() %>% tolower()
    
    # Climate variables from dictionary
    val <- vr_tbl$Category[vr_tbl$Classification == 'Climate']
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
    
    result <- data.frame(Variable = vr_tbl$Variable[vr_tbl$Classification == 'Climate'],
                         Code     = vr_tbl$Code[vr_tbl$Classification == 'Climate'],
                         Selected = selected)
    result <- result %>% tidyr::drop_na() 
    write.csv(result, dest_dir, row.names = F)
    
  } else {
    result <- read.csv(dest_dir, header = T)
  }
  
  return(result)
  
}
# View(result)
# Execute climate clusters code ...
