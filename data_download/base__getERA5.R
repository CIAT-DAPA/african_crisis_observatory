# CDS dataset description
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/sis-agrometeorological-indicators?tab=overview

# Has a number of dependencies, working on alternatives
# devtools::install_github("bluegreen-labs/ecmwfr")
# install.packages("ecmwfr")
library("ecmwfr")

# credentials
UID = "129570"
key = "a2dd82bf-3cd2-4216-9f8f-610a2409e349"

# save key for CDS
wf_set_key(user = UID,
           key = key,
           service = "cds")

getERA5 <- function(i, qq, year, month, datadir){
  q <- qq[i,]
  format <- "zip"
  ofile <- paste0(paste(q$variable, q$statistics, year, month, sep = "-"), ".",format)
  
  if(!file.exists(file.path(datadir,ofile))){
    ndays <- lubridate::days_in_month(as.Date(paste0(year, "-" ,month, "-01")))
    ndays <- 1:ndays
    ndays <- sapply(ndays, function(x) ifelse(length(x) == 1, sprintf("%02d", x), x))
    ndays <- dput(as.character(ndays))
    
    cat("Downloading", q[!is.na(q)], "for", year, month, "\n"); flush.console();
    
    request <- list("dataset_short_name" = "sis-agrometeorological-indicators",
                    "variable" = q$variable,
                    "statistic" = q$statistics,
                    "year" = year,
                    "month" = month,
                    "day" = ndays,
                    "area" = "90/-180/-90/179.9", # download global #c(ymax,xmin,ymin,xmax)? 
                    "time" = q$time,
                    "format" = format,
                    "target" = ofile)
    
    request <- Filter(Negate(anyNA), request)
    
    file <- wf_request(user     = UID,   # user ID (for authentification)
                       request  = request,  # the request
                       transfer = TRUE,     # download the file
                       path     = datadir)  
  } else {
    cat("Already exists", q[!is.na(q)], "for", year, month, "\n"); flush.console();
  }
  return(NULL)
}

########################################################################################################
# change data directory
datadir <- "//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/ERA5"
# datadir <- ""
dir.create(datadir, FALSE, TRUE)

# combinations to download
qq <- data.frame(variable = c("solar_radiation_flux",rep("2m_temperature",3),
                              "10m_wind_speed", "2m_relative_humidity"),
                 statistics = c(NA, "24_hour_maximum", "24_hour_mean", "24_hour_minimum",
                                "24_hour_mean", NA),
                 time = c(NA,NA,NA,NA,NA, "12_00"))

# temporal range
syear <- 2024
eyear <- 2024
years <- as.character(2024:format(Sys.time(), "%Y"))
months <- c(paste0("0", 1:9), 10:12)

# all download
for (i in 1:nrow(qq)){
  for (year in years){
    for (month in months){
      tryCatch(getERA5(i, qq, year, month, datadir), error = function(e)NULL)
    }
  }
}

# for (year in as.character(1999:2001)){
#   for (month in months){
#     getERA5(i=5, qq, year, month, datadir)
#   }
# }


# unzip
zz <- list.files(datadir, ".zip$", full.names = TRUE)

vars <- c("solar_radiation_flux","10m_wind_speed","2m_temperature-24_hour_maximum",
          "2m_temperature-24_hour_mean","2m_temperature-24_hour_minimum","2m_relative_humidity")

extractNC <- function(var, zz, datadir, ncores = 1){
  z <- grep(var, zz, value = TRUE)
  fdir <- file.path(dirname(datadir), var)
  dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  parallel::mclapply(z, function(x){unzip(x, exdir = fdir)}, mc.cores = ncores, mc.preschedule = FALSE)
  return(NULL)
} 

for (var in vars){
  extractNC(var, zz, datadir, ncores = 1)
}

