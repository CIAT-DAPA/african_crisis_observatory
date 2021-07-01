lowest_gadm <- function(iso = 'KEN', out = NULL){
  suppressMessages(library(raster))
  levels <- 5:1
  for(i in 1:length(levels)){
    tryCatch(expr = {
      shp <- raster::getData(name = 'GADM', country = iso, level = levels[i])
      break
    },
    error = function(e){
      cat(paste0("Getting GADM level ",levels[i]," failed... Trying a higher level\n"))
      return("\n")
    })
  }
  if(!is.null(out)){
    raster::shapefile(shp, out)
  }
  return(shp)
}
