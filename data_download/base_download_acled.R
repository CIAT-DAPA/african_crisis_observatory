### Script to donwload ACLED data
#' @author: Ewaut and Benson Kenduiywo
rm(list=ls(all=TRUE))
out <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/_global/conflict/"

require(httr)

# Documentation:    https://acleddata.com/acleddatanew/wp-content/uploads/2020/10/ACLED_API-User-Guide_2020.pdf
# isoCodes:         https://www.acleddata.com/download/3987/ 


## Credentials ############################

KEY = "scnE%21a2zulyDKmwO3-Gj"
MAIL = 'climatesecurity@cgiar.org'


## Functions ##############################

.downloadPages = function(url)
{
  baseUrl = url
  
  # Init
  page = 1
  nRow = 5000
  dtOut = data.frame()
  
  maxRow = 5000
  while(nRow == maxRow)
  {
    cat(sprintf("Downloading results %s to %s\n", (page-1)*maxRow + 1, page * maxRow))
    url = sprintf("%s&page=%s", baseUrl, page)
    res = httr::GET(url)
    dt = data.frame(httr::content(res))
    dtOut = rbind(dtOut, dt)
    nRow = nrow(dt)
    page = page + 1
  }
  
  return(dtOut)
  
}

# Get all regions, codes and number of events
regions = function()
{
  
  res = httr::GET(
    sprintf(
      "https://api.acleddata.com/region/read.csv?terms=accept&key=%s&email=%s",
      KEY,
      MAIL
    )
  )
  
  data.frame(httr::content(res))
  
}

# Get all countries, codes and number of events
countries = function()
{
  
  res = httr::GET(
    sprintf(
      "https://api.acleddata.com/country/read.csv?terms=accept&key=%s&email=%s",
      KEY,
      MAIL
    )
  )
  
  data.frame(httr::content(res))
  
}

#' @param region filter for region code as listed in the documentation
#' @param iso filter for iso code as listed in the iso code document
#' @param updatedSince only download data that was updated since the selected date, e.g. for more efficient querying in the future
#' @param field Only returns the data columns in the field argument
events = function(region = NULL, iso = NULL, updatedSince = NULL, field = NULL)
{
  
  # Build URL
  url = sprintf(
    "https://api.acleddata.com/acled/read.csv?terms=accept&key=%s&email=%s",
    KEY,
    MAIL
  )
  
  if(!is.null(region)) url = sprintf("%s&region=%s", url, paste0(region, collapse = "|"))
  if(!is.null(iso)) url = sprintf("%s&iso=%s", url, paste0(iso, collapse = "|"))
  if(!is.null(updatedSince)) url = sprintf("%s&timestamp=%s", url, updatedSince)
  if(!is.null(field)) url = sprintf("%s&fields=%s", url, paste0(field, collapse = "|"))
  
  # Download URL data by going through the pages
  return(.downloadPages(url))
  
}



## Examples ##############################

## List regions and codes
dtRegion = regions() 
head(dtRegion)
cat('Region names:', sort(unique(dtRegion$region_name)), sep='\n')
reg.names <- sort(unique(dtRegion$region_name))

#' Download several countries in a region
#' Note this may take some time  donwload
donwloadRegion <- function(regionName){
  regionCodes <- dtRegion[grep(regionName, dtRegion$region_name), 'region']
  #' Donwload  All events for region
  return(events(region = regionCodes)) # Will take a long time
}

regionName <- 'East-Asia-Pacific'
eastAsia  <- donwloadRegion('East Asia')
centralAsia <- donwloadRegion("Caucasus and Central Asia")
southAsia <- donwloadRegion("South Asia")
southEastAsia <- donwloadRegion("Southeast Asia")
oceania <- donwloadRegion("Oceania")
asiaPacific <- do.call(rbind, list(eastAsia, centralAsia, southAsia, southEastAsia, oceania))
names(asiaPacific) <- toupper(names(asiaPacific))
library("writexl")
write_xlsx(asiaPacific, paste0(out, regionName, '.xlsx'))
#data.table::fwrite faster for csv

regionName <- 'Africa'
WA  <- donwloadRegion('Western Africa')
SA <- donwloadRegion("Southern Africa")
MA <- donwloadRegion("Middle Africa")
Na <- donwloadRegion("Northern Africa")
EA <- donwloadRegion("Eastern Africa")
africa <- do.call(rbind, list(WA, SA, MA, Na, EA))
names(africa) <- toupper(names(africa))
#write_xlsx(africa, paste0(out, regionName, '_', Sys.Date(),'.xlsx'))
write_xlsx(africa, paste0(out, regionName, '.xlsx'))

regionName <- 'MENA'
ME  <- donwloadRegion('Middle East')

mena <- do.call(rbind, list(Na, ME))
names(mena) <- toupper(names(mena))
write_xlsx(mena, paste0(out, regionName, '.xlsx'))

regionName <- 'Americas'
cAmer  <- donwloadRegion('Central America')
sAmer  <- donwloadRegion('South America')
nAmer  <- donwloadRegion('North America')

americas <- do.call(rbind, list(cAmer, sAmer, nAmer))
names(americas) <- toupper(names(americas))
write_xlsx(americas, paste0(out, regionName, '.xlsx'))

regionName <- 'Europe'
EU  <- donwloadRegion('Europe')

names(EU) <- toupper(names(EU))
write_xlsx(EU, paste0(out, regionName, '.xlsx'))





#' List all countries
#dtCountries = countries()
#head(dtCountries)

#' List All events for the world
#data1 = events() # Will take a long time
#str(data1)



#' OR
#data2 = events(region = 1:5) # Will take a long time
#str(data2)

## All events for Africa but only return the data that was updated since yesterday 
#data3 = events(region = 1:5, updatedSince = Sys.Date()-1)
#str(data3)

## All events for Africa but only return the data that was updated since yesterday, and only return 3 columns of selected fields
#data4 = events(region = 1:5, field = c("iso", "fatalities", "event_date"), updatedSince = Sys.Date()-1)
#str(data4)

## All events that occurred in Kenya
#isoKenya = dtCountries[dtCountries$country=="Kenya","iso"]
#data5 = events(iso = isoKenya)
#str(data5)

## All events that occurred in multiple countries
#isoSelectedCountries = dtCountries[dtCountries$country%in%c("Kenya","Tanzania"),"iso"]
#data6 = events(iso = isoSelectedCountries)
#str(data6)
