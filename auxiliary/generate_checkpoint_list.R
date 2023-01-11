pacman::p_load(tidyverse, stringr)

baseDir <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data"


base_files <- readxl::read_excel("//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/Hostpots_data_dictionary.xlsx") %>% 
  dplyr::select(Variable, to_use, Code,Component) %>% 
  dplyr::filter(as.logical(to_use)) %>% 
  dplyr::filter(Component != "Conflict") %>% 
  drop_na(Code) %>% 
  dplyr::mutate()



avaliable_files <- lapply(list.files(baseDir, recursive = F, pattern = "[A-Z]{3}$", full.names = T)
       , function(x){
  av_files <- tibble(fl_pth_av = list.files(x, recursive = T, full.names = T) %>% 
                       .[!grepl("/_",.)] %>% 
                       .[!grepl("old|tmp", .)],
                     fl_name_av = basename(fl_pth_av) %>% 
                       gsub(pattern = ".tif", replacement = "", .) %>% 
                       gsub(pattern = "[A-Z]{3}_rwi", "{iso}_rwi", .) %>% 
                       gsub(pattern = "[A-Z]{3}_AWE", "{iso}_AWE", .)) 
  
  
  res <- base_files %>% 
    dplyr::left_join(., av_files, by = c("Code" = "fl_name_av")) %>% 
    dplyr::select(fl_pth_av)
  
  names(res) <- basename(x)
  return(res)
}) %>% 
  dplyr::bind_cols()



bind_cols(base_files, avaliable_files) %>% View


 # writexl::write_xlsx(., "D:/OneDrive - CGIAR/Documents/AFO_progress.xlsx")




