# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, ecmwfr, terra))

cds.usr <- "63618"
cds.key <- "0168398e-3f9f-4a6a-9430-01d176e61e90"
ecmwfr::wf_set_key(user = cds.usr, key = cds.key, service = "cds")

request <- list(
  format = "zip",
  variable = "2m_temperature",
  statistic = "24_hour_maximum",
  year = "2021",
  month = "12",
  day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31"),
  dataset_short_name = "sis-agrometeorological-indicators",
  target = "D:/tmax.zip"
)

file <- ecmwfr::wf_request(user     = cds.usr,
                           request  = request,
                           transfer = TRUE,
                           path     = "~",
                           verbose  = TRUE)

# Unzip the downloaded file and then read it

root <- 'C:/Users/haachicanoy/Downloads/tmax_ERA5' # I've manually changed the location of the files
r <- terra::rast(x = list.files(path = root, full.names = T))
plot(r[[1]] - 273.15) # To show Celsius degrees
