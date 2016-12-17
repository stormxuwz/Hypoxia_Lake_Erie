# The script contains the several helper functions

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

transformGeo <- function(geo, isGrid = TRUE){
	# function to do projection
	if(!isGrid){
		geo <- arrange(geo,loggerID)
	}
	
	coordinates(geo)=~longitude+latitude
	raster::projection(geo)=CRS("+init=epsg:4326")
	return(geo)
}

createGrid <- function(loggerInfo, by.x = 0.01, by.y = 0.01){
	# create grid based on sampling locations
	longitudeRange <- range(loggerInfo$longitude)
	latitudeRange <- range(loggerInfo$latitude)
	
	grid <- expand.grid(longitude=seq(longitudeRange[1],longitudeRange[2],by = by.x),latitude=seq(latitudeRange[1],latitudeRange[2],by = by.y)) %>% 
	lonlat2UTM()

	convexHullModel <- convHull(loggerInfo[,c("longitude","latitude")])
	

	totalArea = diff(range(grid$x))*diff(range(grid$y)) # may need to modify slightly

	grid$convexIndex <- predict(convexHullModel,grid)
	
	attr(grid,"totalArea") <- totalArea*sum(grid$convexIndex)/nrow(grid)  # km^2
	print(attr(grid,"totalArea"))

	return(grid)
}


lonlat2UTM <- function(xy){
	# xy is a data frame contains longitude and latitude, add UTM x and y
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


st_variogram <- function(zooDF,logger_geo, detrendExpr = "~1",...){
	# zooDF is a zoo object with columns as each logger's data
	# logger_geo has columns of loggerID, lon,lat, bathy, x and y
	
	require(spacetime)
	require(gstat)
	
	args <- list(...)
	logger_geo <- transformGeo(logger_geo)
	logger_time <- as.POSIXct(index(zooDF),origin = "1970-1-1")
	
	n <- ncol(zooDF) # the number of sensors
	T <- nrow(zooDF) # the number of sampling periods
	
	trend <- 0
	
	res <- (zooDF - trend) %>% as.data.frame() 
	
	res$samplingTime <- logger_time
	res <- melt(res, id.vars = c("samplingTime")) %>%
		arrange(samplingTime,variable)
	
	timeDF <- STFDF(sp=logger_geo,time=logger_time,data=data.frame(value = res$value))
	# vST <- variogramST(value~longitude+latitude+bathymetry+I(bathymetry^2),timeDF,tlags=0:4,boundaries=c(0,25,50,75))
	vST <- variogramST(as.formula(paste("value", detrendExpr)), timeDF, tlags = 0:10, boundaries = seq(0,100,10))
	prodSumModel <- vgmST("productSum",space = vgm(1, "Exp", 150, 0.5),time = vgm(1, "Exp", 5, 0.5),k = 50) 
	vST_fit <- fit.StVariogram(vST, prodSumModel, fit.method=6)
	
	return(list(vgmModel = vST, fit_vgmModel = vST_fit, timeDF = timeDF))
}



