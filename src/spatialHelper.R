require(rgdal)
require(raster)
require(dismo)

findBathy <- function(spData,bathyRasterFile){
	coordinates(spData)=~Long+Lat
	bathymetry_raster = raster(bathyRasterFile)
	bathymetry <- extract(bathymetry_raster,spData)
	return(bathymetry)
}


createGrid <- function(loggerInfo){
	longitudeRange <- range(loggerInfo$longitude)
	latitudeRange <- range(loggerInfo$latitude)
	grid <- expand.grid(longitude=seq(longitudeRange[1],longitudeRange[2],by=0.01),latitude=seq(latitudeRange[1],latitudeRange[2],by=0.01))
	convexHullModel<-convHull(loggerInfo[,c("longitude","latitude")])
	convexIndex <- predict(convexHullModel,grid)
	grid <- subset(grid,convexIndex==1)
	return(grid)
}