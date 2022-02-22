# 00

require(pacman)
pacman::p_load(tidyverse, raster, terra,sf, stars, fst, stringi, stringr, lubridate, furrr, purrr, future, ncdf4, trend)

source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/_main_functions.R')         # Main functions
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/_get_soil_data.R')          # Get soil data
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/_calc_indices.R')           # Calculating agro-indices
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/_calc_indices2.R')          # Calculating agro-indices
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/_calc_spi_drought.R')       # SPI calculation
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/time_series_plot.R')        # Time series graphs
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/maps.R')                    # Maps
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/time_series_plot_region.R') # Time series graphs by region
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/climatology_plot.R')        # Climatology graph. 
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/elv_map.R')                 # Elevation map. 
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/summary.R')                 # summary indices (mean, median...)
source('https://raw.githubusercontent.com/CIAT-DAPA/WFP-profiles/main/migration/_get_climate4regions_districts.R') # Filter climate for districts of interest

get_5_dry_spell <- function(prec){
  
  n_days <- rle(prec)$lengths[rle(prec)$values]
  ret <-  sum(round(n_days[n_days >= 5]/5), na.rm = T)
  return(ret)
}


 
home_dir <- "/home/acmendez/"

data_dir <-  "/cluster01/Workspace/ONECGIAR/Data/"
## Defining country parameters
# Country
iso3     <- 'ZWE'     # 'TZA'
season_type = "type_1"



list.files(paste0(home_dir, "era5_extracted/"), full.names = T ) %>% 
  lapply(., function(i){
    x <- fst::read_fst(i)
  }) %>% 
  bind_rows() %>%
  fst::write_fst(., "/home/acmendez/ZWE_climate.fst")

infile  <- paste0(home_dir, "ZWE_climate.fst") # infile <- flt_clm(iso = iso, country = country)
soilfl  <- paste0(home_dir, "soil_fst/input_soil_", iso3, ".fst")
outfile <- paste0(home_dir, "climatic_index/", iso3, "_indices.fst")
spi_out <- paste0(home_dir, "climatic_index/", iso3, "_spis.fst")


climate = infile
soil    = soilfl
ncores  = 15
outfile = outfile



Soil <- soil %>%
  tidyft::parse_fst(path = .) %>%
  tidyft::select_fst(id,x,y,scp,ssat) %>%
  base::as.data.frame()
head(Soil)


if(nrow(Soil[is.na(Soil$scp),]) > 0){
  NAs <- Soil[is.na(Soil$scp),]
  for(i in 1:nrow(NAs)){
    dst <- geosphere::distm(x = Soil[,c('x','y')], y = NAs[i,c('x','y')]) %>% as.numeric()
    Soil$scp[which(Soil$id == NAs$id[i])] <- Soil$scp[which(as.integer(rank(dst)) == 2)]
    Soil$ssat[which(Soil$id == NAs$id[i])] <- Soil$ssat[which(as.integer(rank(dst)) == 2)]
  }; rm(i, NAs)
}

clim_data <- climate %>%
  tidyft::parse_fst(path = .) %>%
  base::as.data.frame()


clim_data$year <- NULL
clim_data <- clim_data %>%
  dplyr::select(id,x,y,date,prec,tmax,tmean,tmin,srad,rh) %>%
  dplyr::mutate(id1 = id) %>%
  tidyr::nest(Climate = c('id','date','prec','tmax','tmean','tmin','srad','rh')) %>% # 'wind'
  dplyr::rename(id = 'id1') %>%
  dplyr::select(id, dplyr::everything(.))


px <- intersect(clim_data$id, Soil$id)
clim_data <- clim_data[clim_data$id %in% px,]



spi <- clim_data %>% tidyr::unnest() %>%
  dplyr::select(-id1) %>%
  dplyr::mutate(year  = lubridate::year(date),
                month = lubridate::month(date)) %>%
  dplyr::group_by(id, year, month) %>%
  dplyr::summarize(ATR = calc_atrMP(PREC = prec)) %>%
  dplyr::arrange(year, month) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(SPI = calc_spi(DP = ATR)$fitted) %>%
  setNames(c('id','month','year','ATR','SPI')) %>%
  dplyr::mutate(SPI = as.numeric(SPI))
fst::write_fst(x = spi, path = spi_out)
rm(spi)


future::plan(future::multicore, workers = 15)

clim_data2 <- clim_data %>%
  dplyr::mutate(climatic_index = furrr::future_map2(.options = furrr_options(seed = TRUE), .x =  Climate, .y=  id, function(.x, .y){
    
    tbl <-  .x
    if(!all(is.na(tbl$tmax))){
      soilcp <- Soil$scp[Soil$id == .y]
      soilst <- Soil$ssat[Soil$id == .y]
      watbal_loc <- watbal_wrapper(tbl, soilcp, soilst)
      watbal_loc$IRR <- watbal_loc$ETMAX - watbal_loc$prec
      tbl <- tbl %>%
        dplyr::mutate(ERATIO  = watbal_loc$ERATIO,
                      TAV     = (watbal_loc$tmin + watbal_loc$tmax)/2,
                      IRR     = watbal_loc$IRR,
                      LOGGING = watbal_loc$LOGGING,
                      GDAY    = ifelse(TAV >= 6 & ERATIO >= 0.35, yes = 1, no = 0))
      
      
      
      if(season_type == "type_1"){
        dts <- unique(lubridate::year(as.Date(tbl$date)))
        seasons <- lapply(1:(length(dts) - 1), function(i){
          ret <- data.frame( date_seqs =  as.character(seq(as.Date(paste0(dts[i],"-06-02")), as.Date(paste0(dts[i+1],"-06-01")), 1)) ,
                             time_spam =  paste0(paste0(dts[i],"-06-02"), "_", paste0(dts[i+1],"-06-01")) )
          return(ret)
        }) %>% bind_rows()
        
        
      }else if(season_type == "type_2"){
        
        seasons <-  data.frame( date_seqs =  tbl$date ,
                                semester = lubridate::semester(as.Date(tbl$date)) ) %>%
          dplyr::mutate(time_spam = paste0(lubridate::year(as.Date(tbl$date)), "_", semester   )) %>%
          dplyr::select(-semester)
        
      }else{
        
        seasons <-  data.frame( date_seqs =  tbl$date ,
                                time_spam =  lubridate::year(as.Date(tbl$date)))
        
      }
      
      x <- lapply(unique(seasons$time_spam), function(k){
        
        df <- tbl %>%
          dplyr::left_join(., seasons, by = c("date" = "date_seqs") ) %>%
          dplyr::filter(time_spam == k) 
        
        df$to_proc <- FALSE
        
        start <- 1:nrow(df)
        end <- start + 120
        stop <- which(end == nrow(df))
        start <- start[1:stop]
        end <- end[1:stop]
        sums <- c()
        for(j in 1:stop){
          sums[j] <- sum(df$prec[seq(start[j], end[j])], na.rm = T)
        }
        
        df$to_proc[seq(start[which.max(sums)], end[which.max(sums)] )] <- TRUE
        
        df <- df %>%
          dplyr::filter(to_proc)
        
        ret <- tibble(                        time_spam = k,
                                              ATR  = calc_atrMP(PREC = df$prec), 
                                              AMT  = calc_amtMP(TMEAN = df$tmean),
                                              NDD  = calc_nddCMP(PREC = df$prec),
                                              P5D  = calc_p5dCMP(PREC = df$prec),
                                              P95  = calc_p95CMP(PREC = df$prec),
                                              NT_X = calc_htsCMP(tmax = df$tmax, t_thresh = 35),
                                              NDWS = calc_wsdays(df$ERATIO, season_ini = 1, season_end = length(df$ERATIO), e_thresh=0.5),
                                              NWLD = calc_NWLDMP(LOGG = df$LOGGING),
                                              NWLD50 = calc_NWLD50MP(LOGG = df$LOGGING, sat = soilst),
                                              NWLD90 = calc_NWLD90MP(LOGG = df$LOGGING, sat = soilst), 
                                              IRR  = sum(df$IRR, na.rm = T),
                                              SHI  = calc_SHIMP(tmax = df$tmax, RH = df$rh),
                                              calc_HSIMP(tmax = df$tmax, RH = df$rh),
                                              calc_THIMP(tmax = df$tmax, RH = df$rh),
                                              CSDI= calc_csdiMP(TMIN = df$tmin),
                                              dry_5_days_spell = get_5_dry_spell(df$prec)
        ) 
        
        return(ret)
        
      }) %>%
        dplyr::bind_rows()
      
    }else{
      
      ret <- tibble(                    time_spam = NA,
                                        ATR  = NA, 
                                        AMT  = NA,
                                        NDD  = NA,
                                        P5D  = NA,
                                        P95  = NA,
                                        NT_X = NA,
                                        NDWS = NA,
                                        NWLD = NA,
                                        NWLD50 =NA,
                                        NWLD90 =NA, 
                                        IRR  = NA,
                                        SHI  = NA,
                                        HSI  = NA, 
                                        HSI_0 = NA, 
                                        HSI_1 = NA, 
                                        HSI_2 = NA, 
                                        HSI_3 = NA,
                                        THI = NA, 
                                        THI_0 = NA,
                                        THI_1 = NA, 
                                        THI_2 = NA, 
                                        THI_3 = NA,
                                        dry_5_days_spell = NA) 
      
    }
    
    
    
  }) 
  
  ) 

future::plan(future::sequential)


saveRDS(clim_data2, "/home/acmendez/climate_idex_halfway.rds")








