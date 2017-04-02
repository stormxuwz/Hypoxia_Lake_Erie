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

# Indicator kriging to determine the probability of whether DO is above thermocline or below thermocline
require(gstat)
require(sp)
require(raster)
source("src/database.R")

checkAboveThreshold <- function(series,threshold){
	return(ifelse(series>threshold,1,0))
}

indicatorKriging <- function(data,sp,grid){
	# data is a data frame with time, value and logger as columns
	coordinates(sp)=~Longitude+Latitude
	data$value <- checkAboveThreshold(data$value)
}

calculateEmpericalVariogram <- function(data,sp){
	var <- names(data)[2]
	data <- dcast(data,Time~logger,value.var=var)
	time <- data[,1]
	data <- data[,-1]
	loggerNames <- names(data)
	
	geoData <- merge(data.frame(loggerNames = loggerNames),sp, by.x = "loggerNames",by.y = "loggerID",all.x = TRUE)
	
	coordinates(geoData) = ~longitude+latitude
	projection(geoData)=CRS("+init=epsg:4326")
	
	combinations <- combn(length(loggerNames),2)
	
	allDist <- spDists(geoData)
	
	pairDiffList <- list()
	for(i in 1:ncol(combinations)){
		pair <- combinations[,i]
		pairNames <- loggerNames[pair]
		pairDiff <- data[,pair[1]]-data[,pair[2]]
		pairDiffList[[i]] <- list(name = pairNames, idx = pair, diff = pairDiff,dist = allDist[pair[1],pair[2]])
	}
	
	# combine all
	vgmAll = data.frame()
	for(content in pairDiffList){
		vgmAll <- rbind(vgmAll,data.frame(diff = content$diff,dist = content$dist)) 
	}
	return(vgmAll)
}

calculateEmpericalVariogram2 <- function(data,sp,i){
	var <- names(data)[2]
	data <- dcast(data,Time~logger,value.var=var)
	time <- data[,1]
	data <- data[,-1]
	loggerNames <- names(data)
	
	geoData <- merge(data.frame(loggerNames = loggerNames),sp, by.x = "loggerNames",by.y = "loggerID",all.x = TRUE)
	
	coordinates(geoData) = ~longitude+latitude
	projection(geoData)=CRS("+init=epsg:4326")
	
	geoData$value <- as.numeric(data[i,])
	v <- variogram(value~longitude+latitude+bathymetry,data = geoData,cloud = TRUE)
	plot(v)
}


thermoclineProbability <- function(){
	var <- "Temp"
	timeRange <- c("2014-07-01","2014-10-05")
	data <- retriveLoggerData(allGeoData$loggerID,year,var,groupRange,dataType,timeRange=timeRange,transform = FALSE)
	sp  <- retriveGeoData(2014,"B")
	data$Temp <- checkAboveThreshold(data$Temp,14)
	
	vgmAll <- calculateEmpericalVariogram(data,sp)
	
	data2 <- dcast(data,Time~logger,value.var="Temp")
	#index <- seq(1,dim(vgmAll)[1],by = 2208)
	#qplot(dist,diff,data = vgmAll[index,])
}


# spatialTemporalKriging <- function(year){
# 	require(SpatioTemporal)
# 	loggerInfo <- retriveGeoData(year,"B")
# 	loggerInfo <- lonlat2UTM(loggerInfo)
# 	DOdata <- retriveLoggerData(loggerInfo$loggerID,year,"DO","daily","AVG",transform = TRUE)
	
# 	DOdata <- na.omit(data)
# 	time <- index(DOdata)
	
# 	DOdata <- as.matrix(DOdata)
# 	rownames(DOdata) <- as.character(time)
	
# 	names(loggerInfo)[1] <- "ID"
# 	DO_class <- createSTdata(obs = DOdata, covars = loggerInfo)
	
# 	# D <- createDataMatrix(DO_class)
# 	# SVD.cv <- SVDsmooth(D,1:4)
	
# 	DO_class <- updateTrend(DO_class,n.basis = 4)
# 	# plot(DO_class,"obs",ID ="10523447")
	
# 	beta.lm <- estimateBetaFields(DO_class)
	
# 	# set up covariance function
# 	cov.beta <- list(covf = "exp", nugget = FALSE)
# 	cov.nu <- list(covf = "exp", random.effect = FALSE)
	
# 	# LUR <- list(~bathymetry,~bathymetry,~bathymetry,~bathymetry)
# 	LUR <- NULL
# 	locations <- list(coordis = c("x","y"),long.lat = c("longitude","latitude"))
# 	DO.model <- createSTmodel(DO_class,LUR = LUR,cov.beta = cov.beta, cov.nu = cov.nu,locations = locations)
	
# 	# parameter estimation
# 	dim <- loglikeSTdim(DO.model)  # dim$nparam.cov = 13, 
# 	x.init <- cbind(c(rep(2,dim$nparam.cov-1),0),
# 									c(rep(c(1,-3),dim$m+1),-3,0))
	
# 	rownames(x.init) <- loglikeSTnames(DO.model,all = FALSE)
# 	est.DO.model <- estimate(DO.model, x.init, type = "p",hessian.all = TRUE)
		
# }


