# download CHIRPS tiles

chirpsTable <- function(){
  
  region <- c(rep("africa", 7), rep("global", 7))
  intervals <- c("daily", "pentad", "dekad", "monthly", "2-monthly", "3-monthly")
  interval <- c("6-hourly", intervals, intervals, "annual")
  format <- c(NA, rep("tifs", 13))
  # add start and end dates (for something like chirps use last date of previous month considering CHC update schedule)
  
  dd <- data.frame(region, interval, format)
  dd$res <- ifelse(dd$interval == "6-hourly", "0.1", NA)
  dd$res <- ifelse(dd$interval == "daily", "0.05,0.25", dd$res)
  
  # https://data.chc.ucsb.edu/products/CHIRPS-2.0/moving_01pentad/tifs/Anomaly_01PentAccum_Current.tif
  # movingpentad <- data.frame(var = "moving", interval = c("01","02","06","12","18"), suffix = "pentad", format = "tifs")
  return(dd)
}

rain_chirps <- function(region = "global", startdate, enddate, res, path, asRaster = FALSE,...) {
  # check directory
  stopifnot(dir.exists(path))
  
  # check arguments
  dd <- chirpsTable()
  if(!region %in% c(unique(dd$region))){ stop("region name not avaiable; see chirpsTable()") }
  
  # check for resolution
  res <- as.character(res)
  reso <- dd$res[dd$region == region & dd$interval == interval]
  if(!is.na(reso)){
    reso <- unlist(strsplit(reso, ","))
    if(!res %in% reso){
      stop("region/interval/resolution combination not avaiable; see chirpsTable()")
    }
    resolution <- gsub("p0.", "p", paste0("p", res))
  }
  
  # check format
  format <- dd$format[dd$region == region & dd$interval == interval]
  
  # directory structure
  folder <- file.path(path, "chirps", paste0(region,"_",interval), format, resolution)
  
  # file urls
  baseurl <- "https://data.chc.ucsb.edu/products/CHIRPS-2.0"
  u <- file.path(baseurl, paste0(region, "_", interval), format, resolution)
  
  # setup file names
  if(interval == "daily"){
    by <- "day"
    seqdate <- seq.Date(as.Date(startdate), as.Date(enddate), by = by)
    years <- format(seqdate, format="%Y")
    dates <- gsub("-","\\.",seqdate)
    if(format == "tifs"){fnames <- file.path(years, paste0("chirps-v2.0.", dates, ".tif.gz"))}
  } else if (interval == ""){
    
  }
  # download urls
  ff <- lapply(fnames, getGZfile, u, folder)
  
  if(asRaster){
    rr <- terra::rast(unlist(ff))
  }
}
getGZfile <- function(fname, u, folder){
  fgz <- file.path(folder, fname)
  dir.create(dirname(fgz), FALSE, TRUE)
  
  if (!file.exists(fgz)) {
    utils::download.file(file.path(u, fname),fgz, mode="wb")
    if (!file.exists(fgz)) {stop("download failed")}
    fz <- try(R.utils::gunzip(fgz, remove = FALSE))
    if (class(fz) == "try-error") {stop("download failed")}
  }
  tfile <- gsub(".gz$", "", fgz)
  return(tfile)
}

setwd('/cluster01/Workspace/ACO/1.Data/chirps')
rain_chirps(region = "global", startdate = '1981-01-01', enddate = '2020-12-31', res = 0.05, path = ".")

# library(RCurl)
# library(xml2)
# library(rvest)
# 
# doc <- rvest::read_html(u)
# hr <- html_attr(html_nodes(doc, "a"), "href")
# # all trailing / are folders, except products 
# hr <- grep("/$", hr, value = T)
# hr <- grep("^/", hr, value = T, invert = TRUE)