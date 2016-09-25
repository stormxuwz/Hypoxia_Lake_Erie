library(sp)
require(rgdal)
require(raster)
require(dismo)

findBathy <- function(spData,bathyRasterFile){
	# find the bathymatry of spData
	coordinates(spData)=~longitude+latitude
	bathymetry_raster = raster(bathyRasterFile)
	bathymetry <- extract(bathymetry_raster,spData)
	return(bathymetry)
}


createGrid <- function(loggerInfo, by.x = 0.01, by.y = 0.01){
	longitudeRange <- range(loggerInfo$longitude)
	latitudeRange <- range(loggerInfo$latitude)
	
	grid <- expand.grid(longitude=seq(longitudeRange[1],longitudeRange[2],by = by.x),latitude=seq(latitudeRange[1],latitudeRange[2],by = by.y)) %>% lonlat2UTM()
	convexHullModel <- convHull(loggerInfo[,c("longitude","latitude")])
	
	attr(grid,"pixSize") = 22
	
	totalArea = diff(range(grid$x))*diff(range(grid$y))

	grid$convexIndex <- predict(convexHullModel,grid)
	
	attr(grid,"totalArea") <- totalArea*sum(grid$convexIndex)/nrow(grid)  # km^2
	print(attr(grid,"totalArea"))

	return(grid)
}


# The script contains the several helper functions
lonlat2UTM <- function(xy){
	# xy is a data frame contains longitude and latitude
	xy0 <- xy 
	coordinates(xy)=~longitude+latitude
	proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
	res <- spTransform(xy, CRS("+proj=utm +zone=17T ellps=WGS84"))
	
	UTMxy <- as.data.frame(coordinates(res))/1000
	names(UTMxy) <- c("x","y")
	
	return(cbind(xy0,UTMxy))
}

scale2Unit <- function(s){
	# scale to 0 to 1
	return((s-min(s,na.rm=T))/(max(s,na.rm=T)-min(s,na.rm=T)))
}     

ft2meter <- function(ft){
	# feet to meter
	return(ft/3.2808399)
}