require(fields)
require(gstat)
require(sp)
require(dplyr)
source("src/spatialHelper.R")
source("src/plot.R")
library(geoR)
source("src/database.R")

spatial_interpolation <- function(df,grid,method = "IDW"){
	# df and grid is a dataframe that contains longitude and latitude and value as columns
	pred <- rep(NA,nrow(grid))

	convexIndex <- grid$convexIndex
	grid <- subset(grid,convexIndex==1)

	if(method == "IDW"){
		coordinates(df) = ~x + y
		coordinates(grid) = ~x + y
		pred[convexIndex == 1] <- krige(value ~  1 , df, grid)$var1.pred
	}
	else if(method == "loglik"){
		# Using log likelihood to fit covariance
		df <- df[,c("x","y","value","bathymetry")] %>% as.geodata()
		ml <- likfit(df, ini = c(5,40), fix.nugget = T, lik.method = "ML",cov.model = "exponential",trend = "1st")
		pred[convexIndex == 1] <- krige.conv(df, locations = grid[,c("x","y")], krige = krige.control(obj.m = ml))$predict
	}else{
		
	}
	
	return(pred)
}



stKriging <- function(df,logger_geo, grid,detrendMethod = NULL,...){
	# df is a zoo object with columns as each logger's data
	# logger_geo has columns of loggerID, lon,lat, bathy, x and y
	args <- list(...)
	coordinates(logger_geo)=~longitude+latitude
	raster::projection(logger_geo)=CRS("+init=epsg:4326")
	
	logger_geo <- arrange(logger_geo,loggerID)
	logger_time <- as.POSIXct(index(df))

	n <- ncol(df) # the number of sensors
	T <- nrow(df) # the number of sampling periods

	if(is.null(detrendMethod)){
		trend <- matrix(0,T,n)

	}else{

	}

	res <- (df - trend) %>% as.data.frame() 
	
	res$samplingTime <- logger_time
	res <- melt(res, id.vars = c("samplingTime")) %>%
			arrange(samplingTime,variable)
			

	timeDF <- STFDF(sp=logger_geo,time=logger_time,data=data.frame(value = res$value))

	vST <- variogramST(value~longitude+latitude+time+bathymetry+I(bathymetry^2),timeDF,tlags=0:4,boundaries=c(0,25,50,75))

}


calulateHypoxiaExtent <- function(data,locationInfo,method = "IDW"){
	# assume the spatial interpolation grid doesn't change along time
	# data is a zoo data frame
	data <- na.omit(data)  # only remain the time where all data are available
	times <- index(data)
	grid <- createGrid(locationInfo)  # will also return an area
	interpolationRes <- matrix(0,nrow = nrow(grid),ncol = nrow(data))
	loggerNames  <- as.numeric(colnames(data))
	print("# of interpolation:")
	print(nrow(data))
	for(i in 1:nrow(data)){
		subData <- data.frame(logger = loggerNames,DO = as.numeric(data[i,]))
		subData <- merge(subData, locationInfo,by.x = "logger",by.y = "loggerID") %>% rename(value = DO)

		grid$pred <- spatial_interpolation(subData,grid)
		interpolationRes[,i] <- ifelse(grid$pred<0,0,grid$pred)
	}

	interpolationRes <- interpolationRes[grid$convexIndex == 1,] # remove the locations that are NA

	totalPx <- nrow(interpolationRes)

	hypoxia_2 <- colSums(interpolationRes<2)/totalPx
	hypoxia_0 <- colSums(interpolationRes<0.01)/totalPx
	hypoxia_4 <- colSums(interpolationRes<4)/totalPx

	hypoxiaExtent <- zoo(data.frame(below_0.01 = hypoxia_0, below_2 = hypoxia_2, below_4 = hypoxia_4),order.by = times)

	# attr(hypoxiaExtent,"pixSize") <- attr(grid,"pixSize")
	attr(hypoxiaExtent,"totalArea") <- attr(grid,"totalArea")

	return(hypoxiaExtent)
}



interpolation_controller <- function(data,locationInfo,method = "IDW"){
	allTimes <- unique(data$Time)
	interpolationResults <- list()
	
	for(i in 1:length(allTimes)){
	#for(i in 1:3){
		hourlyTime = allTimes[i]
		subData <- subset(data,Time == hourlyTime) %>% na.omit()
		subData <- merge(subData, locationInfo,by.x = "logger",by.y = "loggerID") %>% rename(value = DO)
		
		
		grid <- createGrid(subData) # form the grid
		# grid$bathy <- findBathy(grid,"../input/erie_lld/erie_lld.asc")
		
		grid$pred <- spatial_interpolation(subData,grid)
		grid$pred <- ifelse(grid$pred<0,0,grid$pred)
		attributes(grid)$time <- paste(as.character(hourlyTime),"GMT")
		interpolationResults[[i]] <- grid
		print(summary(grid))
	}
	return(interpolationResults)
}


spatialTemporalKriging <- function(year){
	require(SpatioTemporal)
	loggerInfo <- retriveGeoData(year,"B")
	loggerInfo <- lonlat2UTM(loggerInfo)
	DOdata <- retriveLoggerData(loggerInfo$loggerID,year,"DO","daily","AVG",transform = TRUE)
	
	DOdata <- na.omit(data)
	time <- index(DOdata)
	
	DOdata <- as.matrix(DOdata)
	rownames(DOdata) <- as.character(time)
	
	names(loggerInfo)[1] <- "ID"
	DO_class <- createSTdata(obs = DOdata, covars = loggerInfo)
	
	# D <- createDataMatrix(DO_class)
	# SVD.cv <- SVDsmooth(D,1:4)
	
	DO_class <- updateTrend(DO_class,n.basis = 4)
	# plot(DO_class,"obs",ID ="10523447")
	
	beta.lm <- estimateBetaFields(DO_class)
	
	# set up covariance function
	cov.beta <- list(covf = "exp", nugget = FALSE)
	cov.nu <- list(covf = "exp", random.effect = FALSE)
	
	# LUR <- list(~bathymetry,~bathymetry,~bathymetry,~bathymetry)
	LUR <- NULL
	locations <- list(coordis = c("x","y"),long.lat = c("longitude","latitude"))
	DO.model <- createSTmodel(DO_class,LUR = LUR,cov.beta = cov.beta, cov.nu = cov.nu,locations = locations)
	
	# parameter estimation
	dim <- loglikeSTdim(DO.model)  # dim$nparam.cov = 13, 
	x.init <- cbind(c(rep(2,dim$nparam.cov-1),0),
									c(rep(c(1,-3),dim$m+1),-3,0))
	
	rownames(x.init) <- loglikeSTnames(DO.model,all = FALSE)
	est.DO.model <- estimate(DO.model, x.init, type = "p",hessian.all = TRUE)
		
}


# Test
testFunc <- function(){
	year <- 2014
	loggerInfo <- retriveGeoData(year,"B")
	data <- retriveLoggerData(loggerInfo$loggerID,year,"DO","daily","AVG",transform = TRUE) %>% na.omit()
	interpolation_controller(data,loggerInfo) %>% 
	plot_spatial(locationInfo = loggerInfo,outputFolder = "../output/2015spPlots/")
}

