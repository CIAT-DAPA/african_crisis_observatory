##kmeans
set.seed(2000)
xx=kmncluster <- kmeans(na.omit(aa), centers = 3, iter.max = 500, nstart = 5, algorithm="Lloyd")
##========================
aa=as.data.frame(to_cluster)
head(aa)
aa=subset(aa, select=-c(geometry, ov))
library(randomForest)
rf2 <- randomForest(x = aa, mtry = 2, ntree = 300, proximity = TRUE)
prox <- rf2$proximity
pam.rf <- pam(prox, 3)
pam.km <- stats::kmeans(prox, 3)
bb=aa
bb$clust <- pam.rf$clustering
table(bb$clust)

x11()
boxplot(EVENTS~clust, data=bb)
x11()
boxplot(knl~clust, data=bb)
table(bb$clust)

##Coupling kmeans with RF
library(randomForest)
ff <- as.data.frame(pca_w$ind$coord)
ff$clust <- as.factor(km$cluster)
library(caret)
train.index <- createDataPartition(ff$clust, p = .7, list = FALSE)
train <- ff[ train.index,]
test  <- ff[-train.index,]

#aa$clust <- as.factor(km$cluster)
rf3.model <- randomForest(clust~., data=train, mtry = 2, ntree = 300, proximity = TRUE)

rf3 <- predict(rf3.model, ff)
table(rf3)

##Validation
rf4 <- predict(rf3.model, test)
conf <- table(rf4,test$clust)
#Overall accuracy
oa <- sum(diag(conf))/sum(conf)
print(paste('Overal accuracy ', round(oa*100,2)))



#================================================
shp <- rgeos::gUnaryUnion(shp)
xx <- rgeos::gBuffer(shp, width = 0.01)


results[[i]] <- sp.df@polygons[[1]]@Polygons[[i]]@coords

n <- length(shp@polygons[[1]]@Polygons)
shp@polygons[[1]]@Polygons[[1]]@coords <- GrowPolygon(shp@polygons[[1]]@Polygons[[1]]@coords ,buffer= 0.01)

yy <- GrowPolygon(shp@polygons[[1]]@Polygons[[1]]@coords ,buffer= 0.005)
plot(shp@polygons[[1]]@Polygons[[1]]@coords, type = "l", 
     xlim = c(119.43, 119.49),
     ylim = c(4.5, 4.7))
lines(yy$x, yy$y, col = "red")


shp <- rgeos::gUnaryUnion(shp)
shp2 <- shp
for(i in 1:length(shp2@polygons[[1]]@Polygons)){
  res <- GrowPolygon(shp2@polygons[[1]]@Polygons[[i]]@coords, buffer = 0.005)
  shp2@polygons[[1]]@Polygons[[i]]@coords <- cbind(res$x, res$y)
}

shapefile(shp2, paste0(root,'Buffer_PHL.shp'))

x11();plot(shp);lines(shp2, col = "red")
PolygonCentre <- function(x, y = NULL) {
  xy <- xy.coords(x, y)
  area <- PolygonArea(xy, positive = FALSE)
  
  x <- xy$x
  x1 <- c(x[-1], x[1])
  y <- xy$y
  y1 <- c(y[-1], y[1])
  prod <- (x * y1) - (x1 * y)
  
  cbind(x = sum((x + x1) * prod),
        y = sum((y + y1) * prod)) / (6 * area)
}

PolygonArea <- function(x, y = NULL, positive = TRUE) {
  xy <- xy.coords(x, y)
  x <- c(xy$x, xy$x[1])
  y <- c(xy$y, xy$y[1])
  area <- (sum((xy$x * y[-1])) - sum((xy$y * x[-1]))) / 2
  
  # Return:
  if (isTRUE(positive)) {
    abs(area)
  } else {
    area
  }
}

GrowPolygon <- function(x, y = NULL, buffer = 0) {
  xy <- xy.coords(x, y)
  cent <- PolygonCentre(x, y)
  
  x0 <- xy$x - cent[, "x"]
  y0 <- xy$y - cent[, "y"]
  hyp <- sqrt((x0 * x0) + (y0 * y0))
  stretch <- (hyp + buffer) / hyp
  x1 <- x0 * stretch
  y1 <- y0 * stretch
  
  xy$x <- x1 + cent[, "x"]
  xy$y <- y1 + cent[, "y"]
  
  xy
}
