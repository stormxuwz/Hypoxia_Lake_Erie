require(fields)
require(gstat)
require(sp)
require(dplyr)
source("src/spatialHelper.R")
source("src/plot.R")

spatial_interpolation <- function(df,grid,method = "IDW"){
	# df and grid is a dataframe that contains longitude and latitude and value as columns
	pred <- rep(NA,nrow(grid))

	convexIndex <- grid$convexIndex
	grid <- subset(grid,convexIndex==1)

	coordinates(df) = ~longitude + latitude
	coordinates(grid) = ~longitude + latitude
	
	#print(df)
	#print(grid)

	if(method == "IDW"){
		pred[convexIndex == 1] <- krige(value ~  1 , df, grid)$var1.pred
	}
	else{
		# pred = NA
	}
	
	return(pred)
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
	year <- 2015
	loggerInfo <- retriveGeoData(year,"B")
	data <- retriveLoggerData(loggerInfo$loggerID,year,"DO","daily","AVG",transform = FALSE)
	interpolation_controller(data,loggerInfo) %>% 
	plot_spatial(locationInfo = loggerInfo,outputFolder = "../output/2015spPlots/")
}

