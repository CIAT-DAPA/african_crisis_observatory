#' ############################################################################################
#' Author: Benson Kenduiywo
#' ############################################################################################

rm(list=ls(all=TRUE))
g <- gc(reset = T); 
#install.packages("remotes")
#remotes::install_github("rspatial/luna")
list.of.packages <- c("luna","terra","geodata")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)

path <-  '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO2/FCM/Data/raw/avhrr/'
#path <- '//CATALOGUE.CGIARAD.ORG/WFP_ClimateRiskPr/1.Data/'
start <- "1981-01-01" 
end <- "2024-03-04"


#' List all AVHRR files available for download
#' This could be an internal function to download AVHRR files
#'	
.listAVHRR <- function(path, update=FALSE, baseurl) {
  cat("Creating index of available AVHRR files on", as.character(Sys.Date()), "\n")
  # Two-day delay in ingestion
  filename <- paste0("avhrr_files_", Sys.Date(),".rds")
  filename <- file.path(path, filename)
  
  if (!file.exists(filename) | update){
    startyear <- 1981
    endyear <- format(as.Date(Sys.Date(), format="%d/%m/%Y"),"%Y")
    years <- seq(startyear, endyear)
    
    ff <- list()
    for (i in 1:length(years)){
      url <- file.path(baseurl, years[i])
      wp <- xml2::read_html(url)
      dvns <- rvest::html_attr(rvest::html_nodes(wp, "a"), "href")
      #VIIRS-Land_v001-preliminary_NPP13C1_S-NPP_20190101_c20220418131738.nc	
      
      
      ds <- grep("AVHRR-Land|VIIRS-Land_*.*.nc", dvns, value = TRUE)#grep("^AVHRR-Land_*.*.nc", dvns, value = TRUE)
      ff[[i]] <- ds
    }
    ff <- unlist(ff)
    dates <- sapply(strsplit(ff,"_"), "[[", 5)
    dates <- as.Date(dates, format = "%Y%m%d")
    ff <- data.frame(filename = ff, date = dates, stringsAsFactors = FALSE, row.names = NULL)
    saveRDS(ff, filename)
  } else {
    ff <- readRDS(filename)
  }
  return(ff)
}


getAVHRR <- function(start_date, end_date, path, overwrite = FALSE, update = FALSE, ...) {
  
  if(missing(start_date)) stop("provide a start_date")
  if(missing(end_date)) stop("provide an end_date")
  
  baseurl <- "https://www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access"
  # url to access 8 different ways of downloading the data
  # baseurl <- "https://www.ncei.noaa.gov/thredds/catalog/cdr/ndvi/files"
  #path <- .getPath(path)
  
  # list of AVHRR files
  pp <- .listAVHRR(path = path, baseurl = baseurl, update = FALSE)
  
  # TODO: alternate search through CMR
  # https://cmr.earthdata.nasa.gov/search/concepts/C1277746140-NOAA_NCEI
  
  # subset the files by dates
  pp <- pp[pp$date >= start_date & pp$date <= end_date, ]
  
  if(nrow(pp) == 0) {stop("No AVHRR file available for the date range provided")}
  
  # to store output file names
  
  for (i in 1:nrow(pp)){
    ff <- pp[i,]
    fname <- ff$filename
    #year <- .yearFromDate(ff$date)
    year <- format(as.Date(ff$date, format="%d/%m/%Y"),"%Y")
    furl <- file.path(baseurl, year, fname)
    filename <- file.path(path, fname)
    
    # is ok, if file exists or overwrite is TRUE
    #ok <- (file.exists(filename) | overwrite)
    
    # what if the download is bad; less than 50 mb
    # there must be a better way
    if(file.exists(filename)){
      fsz <- round(file.size(filename)/(1024^2))
      if (fsz < 50) ok <- FALSE
    }
    
    if (!file.exists(filename)){
      cat("Downloading AVHRR tile for", as.character(ff$date), "\n")
      ff <- try(utils::download.file(furl, filename, mode = "wb", quiet = TRUE)) 
    } 
    
    if (inherits(ff, "try-error")) next
  }
}

getAVHRR(start_date=start, end_date= end, path = path)
