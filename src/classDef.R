# Define the necessary DO classes
source("src/database.R")
library(dplyr)
# define the generic function

getLakeDO <- function(year, depthLocation = "B", timeAggType = "hourly", ...){
	loggerInfo <- retriveGeoData(year,depthLocation) %>% arrange(loggerID)
	loggerInfo$loggerID <- as.character(loggerInfo$loggerID)

	timeRange <- list(...)$timeRange
	data <- 
		retriveLoggerData(loggerInfo$loggerID,year,"DO",timeAggType,"AVG",timeRange = timeRange,transform = TRUE)  # %>% 
		# na.omit()  # remove the missing data
	data <- data[,loggerInfo$loggerID]
	data <- data * (as.matrix(data) > 0) # put it to zero

	# add a filter to remove the wrong data
	if(year == 2015){
		mask <- index(data) > as.POSIXct("2015-07-23",tz = "GMT")
		data <- data[mask,]
	}

	lakeDO <- as.lakeDO(data, loggerInfo)

	return(lakeDO)
}

as.lakeDO <- function(samplingData, loggerInfo){
	X <- list(samplingData = samplingData, loggerInfo = loggerInfo)
	stopifnot(all(names(X$samplingData) == X$loggerInfo$loggerID)) # make sure the order is correct
	class(X) <- "lakeDO"
	return(X)
}

na.omit.lakeDO <- function(lakeDOObj){
	lakeDOObj$samplingData <- na.omit(lakeDOObj$samplingData)
	stopifnot(length(unique(diff(index(lakeDOObj$samplingData))))==1)
	return(lakeDOObj)
}

subset.lakeDO <- function(X, loggerID_not){
	# remove certain id
	print("removing logger:")
	print(loggerID_not)
	loggerID_not <- as.character(loggerID_not)
	remainLogger <- !names(X$samplingData) %in% loggerID_not
	X$samplingData <- X$samplingData[,remainLogger]
	X$loggerInfo <- X$loggerInfo[remainLogger,]
	stopifnot(all(names(X$samplingData) == X$loggerInfo$loggerID))
	return(X)
}

idwModel <- function(x, metaFolder, nmax){
	# x is the lakeDO class
	# separate X into different time steps
	
	# data <- list([[time1]]:df[x,y, value])
	if(!is.null(metaFolder)) createFolder(metaFolder)
	model <- list()
	loggerInfo <- x$loggerInfo
	for(i in 1:nrow(x$samplingData)){
		loggerInfo$value <- as.numeric(x$samplingData[i,])
		model[[i]] <- na.omit(loggerInfo)
	}

	config <- list(simNum = 1, trend = "value ~ 1", nmax = nmax)
	class(model) <- "krigModelIdw"
	attr(model,"metaFolder") <- metaFolder
	attr(model,"config") <- config
	return(model)
}


basisModel <- function(x, trend, fitMethod, r, metaFolder){
	# x is a object of "lakeDO" class
	# separate X into different prediction parts
	# x$samplingData <- na.omit(x$samplingData)
	x$samplingData <- na.omit(x$samplingData)

	createFolder(metaFolder)
	if(!fitMethod %in% c("Reml", "Baye")){
		stop("fitMethod is not implemented")
	}
	# create the residuals
	
	
	SVDmodel <- SVD_basis(x$samplingData, r)
	basisCoef <- SVDmodel$coef
	basis <- SVDmodel$basis
	loggerInfo <- x$loggerInfo

	# create residuals
	residuals <- x$samplingData - SVDmodel$fit	
	residuals <- as.lakeDO(residuals, loggerInfo)

	# create models
	model <- list()
	for(i in 1:ncol(basis)){
		loggerInfo$value <- basisCoef[i,]
		model[[i]] <- loggerInfo
	}

	# add attributes
	config <- list(
		trend = trend,
	 	simNum = 100, 
	 	cov.model = "exponential")

	class(model) <- paste0("krigModel", fitMethod)
	attr(model,"metaFolder") <- metaFolder
	attr(model, "config") <- config
	attr(model, "basis") <- basis
	attr(model, "residuals") <- residuals 
	return(list(model = model, residuals = residuals))
}

