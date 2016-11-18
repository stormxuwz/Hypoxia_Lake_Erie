# rm(list = ls())
require(fields)
require(gstat)
require(sp)
require(dplyr)
require(geoR)
require(spacetime)
source("src/spatialHelper.R")
source("src/plot.R")
source("src/database.R")
source("src/basisDecomposition.R")

library(geoR)

basis_interpolation_main <- function(data, locationInfo, basis, simNum, intMethod, r){
	
	if(r>nrow(locationInfo)){
		throw("Too large r, please select r<= # of location")
	}

	grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)
	trendPrediction <- basis_interpolation(data,locationInfo,grid,basis = basis,simNum = simNum, intMethod = intMethod, r = r)

	# save the intermediate results
	saveRDS(list(basisIntRes = trendPrediction, grid = grid),"../meta/basis_prediction.rds")
	
	trendSim <- trendPrediction$trend 	# trend simulations, shape of (simNum, TimeN, nrow(grid)
	varExpl <- trendPrediction$varExpl  # variance explained

	print("doing spatial temporal kriging on residuals")
	#residual_stKriging <- spatial_temporal_krigingInterpolation(residual, locationInfo, grid = grid)
	residual_idw <- idw_interpolation_main(trendPrediction$res, locationInfo)

	saveRDS(residual_idw,"../meta/residual_prediction.rds")

	for(i in 1:simNum){
		trendSim[i,,] <- trendSim[i,,] + residual_idw[1,,]
	}
	# return(list(vgmModel = vST, fit_vgmModel = vST_fit, timeDF = timeDF))
	saveRDS(trendSim,"../output/final_prediction.rds")
	return(trendSim)
}



random_interpolation <- function(data, locationInfo, nsim = 100){
	# function to interpolate the map based on random values. i.e. no correlation
	grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)
	data <- as.numeric(data)

	for(i in 1:nsim){}

}


summarySimulation <- function(hypoxiaExtentPred,summaryName){

	variance <- sqrt(apply(hypoxiaExtentPred, 2, var))
	hypoxiaSummary <- data.frame(m = colMeans(hypoxiaExtentPred))
	hypoxiaSummary$upper <- hypoxiaSummary$m+2*variance
	hypoxiaSummary$lower <- hypoxiaSummary$m-2*variance

	names(hypoxiaSummary) <- c(paste(summaryName,c("Mean","Upper","Lower"),sep=""))
	return(hypoxiaSummary)
}


summaryHypoxia <- function(predictionMatrix, timeIndex){
	# predictionMatrix is a 3D array that has shape of c(simNum, TimeN, nrow(grid)))
	predictShape <- dim(predictionMatrix)
	simNum <- predictShape[1]
	TimeN <- predictShape[2]
	locationNum <- predictShape[3]

	# initialize the hypoxia extent corresponding to below 0.01, 2 and 4 mg/L
	hypoxiaExtent_0 <- matrix(0,simNum,TimeN)  
	hypoxiaExtent_2 <- hypoxiaExtent_0	
	hypoxiaExtent_4 <- hypoxiaExtent_0

	for(i in 1:simNum){
		hypoxiaExtent_0[i,] <- rowSums(predictionMatrix[i,,]<0.01,na.rm = TRUE)
		hypoxiaExtent_2[i,] <- rowSums(predictionMatrix[i,,]<2,na.rm = TRUE)
		hypoxiaExtent_4[i,] <- rowSums(predictionMatrix[i,,]<4,na.rm = TRUE)
	}

	if(simNum>1){
		# predictionMatrix contains multiple simulations, with uncertainty
		hypoxiaExtent <- summarySimulation(hypoxiaExtent_0,"less0") %>% 
						cbind(summarySimulation(hypoxiaExtent_2,"less2")) %>%
						cbind(summarySimulation(hypoxiaExtent_4,"less4")) %>%
						zoo(order.by = timeIndex)
	}
	else{
		# predictionMatrix is a deterministic interpolation
		hypoxiaExtent <- data.frame(less0 = hypoxiaExtent_0[1,],less2 = hypoxiaExtent_2[1,],less4 = hypoxiaExtent_4[1,]) %>% 
		zoo(order.by = timeIndex)
	}

	return(hypoxiaExtent)
}

idw_interpolation_main <- function(data, locationInfo){
	data <- na.omit(data)  # only remain the time where all data are available
	times <- index(data)
	grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)  # will also return an area
	
	prediction <- array(0, dim = c(1, length(times), nrow(grid)))

	loggerNames  <- as.numeric(colnames(data))
	print("# of interpolation:")
	print(nrow(data))

	for(i in 1:nrow(data)){
		subData <- data.frame(logger = loggerNames,DO = as.numeric(data[i,]))
		subData <- merge(subData, locationInfo,by.x = "logger",by.y = "loggerID") %>% rename(value = DO)

		pred <- spatial_interpolation(subData,grid)[,1]
		prediction[1,i,] <- ifelse(pred<0,0,pred)
	}

	return(prediction)
}

spatial_interpolation <- function(df,grid,method = "IDW", simNum = 1,tmpName = "tmp",...){
	# df and grid is a dataframe that contains longitude and latitude and value as columns
	pred <- array(NA,dim=c(nrow(grid),simNum)) # shapes as the [gridNum, simNum]
	convexIndex <- grid$convexIndex
	grid <- subset(grid,convexIndex==1)
	
	# check the input validation
	if(method=="IDW" & simNum!=1){
		throw("IDW shouldn't should have simNum =1")
	}

	if(simNum<1){
		throw("Sim num should be positive")
	}

	if(method == "IDW"){
		coordinates(df) = ~x + y
		coordinates(grid) = ~x + y
		pred[convexIndex == 1,1] <- krige(value ~  1 , df, grid)$var1.pred
	}
	else if(method == "loglik"){
		# Using log likelihood to fit covariance
		df <- df[,c("x","y","value","bathymetry")] %>% as.geodata(covar.col = 4)
		# trend <- as.formula(~bathymetry + I(bathymetry^2))
		trend <- "cte"
		ml <- likfit(df, ini = c(1000,70), fix.nugget = T, lik.method = "REML",cov.model = "spherical",trend = trend)  # the nugget is defaulted as zero
		
		# plot the variogram and data
		png(paste("../meta/",tmpName,"loglik_df.png",sep = ""))
		print(plot(df))
		dev.off()
		
		png(paste("../meta/",tmpName, "loglik_variogram.png",sep=""))
		print(plot(variog(df,trend = trend,option = "cloud")))
		print(lines(ml))
		dev.off()

		# kriging control
		KC <- krige.control(obj.m = ml, trend.d = trend, trend.l = trend)
		
		if(simNum ==1){
			pred[convexIndex == 1, 1] <- krige.conv(df, locations = grid[,c("x","y")], krige = KC)$predict
		}else{
			# do conditional simulation
			pred[convexIndex == 1, ] <- krige.conv(df, locations = grid[,c("x","y")], krige = KC, output = output.control(n.predictive = simNum))$simulations
		}
		
	}
	else if(method == "baye"){
		# conduct the bayesian kriging
		trend = "cte"
		MC <- model.control(cov.model = "exponential", trend.l = trend, trend.d = trend) 
		
		# specify the priors
		PC <- prior.control(phi.discrete=seq(30,100,10),  # range is discreted
							beta.prior = "flat",  # beta is flat 
							sigmasq.prior = "reciprocal",
							tausq.rel.prior = "fixed",
							tausq.rel = 0)  # sigma^2 is 
		 
		# specify the output control
		OC <- output.control(n.pos = simNum, # the number of samples taking from posterior distribution
							n.pred = simNum,   # sample to taken from the predictive distribution
							 signal = FALSE)
		
		df <- df[,c("x","y","value","bathymetry")] %>% as.geodata(covar.col = 4)
		
		predRes <- krige.bayes(geodata = df, 
								 locations = grid[,c("x","y")], 
								 model = MC,
								 prior = PC,
								 output = OC)
		pred[convexIndex == 1,] <- predRes$predictive$simulations
		
		png(paste("../meta/",tmpName,"baye_df.png",sep = ""))
		print(plot(df))
		dev.off()
		my.summary <- function(x){quantile(x, prob = c(0.05, 0.5, 0.95))}
		
		png(paste("../meta/",tmpName, "baye_variogram.png",sep=""))
		print(plot(variog(df,option = "cloud",trend = trend)))
		print(lines(predRes, summ = my.summary, ty="l", lty=c(2,1,2), col=1))
		dev.off()
		
	}
	else{
		# specify the covariance matrix
	}
	return(pred)
}


transformGeo <- function(geo, isGrid = TRUE){
	if(!isGrid){
		geo <- arrange(geo,loggerID)
	}

	coordinates(geo)=~longitude+latitude
	raster::projection(geo)=CRS("+init=epsg:4326")
	return(geo)
}


spatial_temporal_krigingInterpolation <- function(data, locationInfo,grid = NULL){
	stv <- st_variogram(data,locationInfo)
	saveRDS(stv,file = "../meta/st_variogram.rds")
	
	timeDF <- stv$timeDF
	vgmModel <- stv$vgmModel
	fit_vgmModel <- stv$fit_vgmModel

	logger_time <- as.POSIXct(index(data),origin = "1970-1-1")

	if(is.null(grid)){
		grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)
	}

	grid <- transformGeo(grid)
	newData <- STFDF(sp=grid,time=logger_time)
	pred <- krigeST(value~1, data = timeDF, modelList= fit_vgmModel, newdata = grid)
	return(pred)
}

st_variogram <- function(zooDF,logger_geo, detrendExpr = "~1",...){
	# zooDF is a zoo object with columns as each logger's data
	# logger_geo has columns of loggerID, lon,lat, bathy, x and y
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



basis_interpolation <- function(DOdata, logger_geo, grid, timeIndex = NULL, basis = "fda",intMethod = "likfit", simNum = 1,...){
	# function to decompose the data into basis and then do interpolation
	# DOdata is the dataframe
	# basis can take "fda" or "svd"
	# time index is a vector that specify how long data is used
	args <- list(...)
	
	ID <- as.numeric(colnames(DOdata))
	TimeN <- nrow(DOdata)  # Totally how many timestamps
	
	if(basis == "fda"){
		decompRes <- B_spline(DOdata.args$knots)
	}else{
		decompRes <- SVD_basis(DOdata,args$r)
		nBasis <- dim(decompRes$coef)[1] # columns of the coefficent matrix
		varExpl <- decompRes$varExpl
	}
	
	basis <- decompRes$basis
	saveRDS(basis, file = "../meta/basis.rds")
	# construct complete coefficient matrix
	coef <- decompRes$coef %>% t() %>% as.data.frame()
	coef$ID <- ID
	coef_df <- merge(coef,logger_geo,by.x = "ID",by.y = "loggerID")
	
	# do interpolation on each basis coefficients
	# prediction <- matrix(0,T,nrow(grid)) # a matrix T * nBasis, row is time, columns is locations
	prediction <- array(0, dim = c(simNum, TimeN, nrow(grid)))
	
	for(i in 1:(nBasis+1)){
		print(sprintf("doing basis %d",i))
		spData <- coef_df[,c("x","y",paste("X",i,sep = ""),"bathymetry","longitude","latitude")]
		names(spData)[3] <- "value"
		
		pred <- spatial_interpolation(df=spData,grid=grid,method=intMethod, simNum =simNum,tmpName = paste("basis_",i,sep=""))
		# pred is matrix with [# of grid, numSim]
		for(sim in 1:dim(pred)[2]){
			prediction[sim,,] <- prediction[sim,,]+basis[,i] %*% t(pred[,sim])
		}
	}
	
	# do interpolation on the residuals
	DO_res <- DOdata - decompRes$fit  # a matrix of shape T, nLocation
	
	return(list(trend = prediction,res = DO_res,varExpl = varExpl))
}


interpolation_main <- function(data,locationInfo,method = "IDW",...){
	timeIndex <- index(data)
	args <- list(...)

	if(method == "IDW"){
		prediction <- idw_interpolation_main(data,locationInfo)
	}
	else if(method == "basis"){
		prediction <- basis_interpolation_main(data,locationInfo, 
			basis = args$basis_method, 
			simNum = args$simNum, 
			intMethod = args$intMethod, 
			r = args$r)
	}else{
		throw("not implemented")
	}

	hypoxiaExtent <- summaryHypoxia(prediction, timeIndex)
	return(hypoxiaExtent)
}


# Test
# testFunc <- function(){
year <- 2015
loggerInfo <- retriveGeoData(year,"B")
data <- retriveLoggerData(loggerInfo$loggerID,year,"DO","hourly","AVG",transform = TRUE) %>% na.omit()  # remove the data
	# IDW_hypoxiaExtent <- interpolation_main(data,loggerInfo,"IDW")
basis_hypoxiaExtent <- interpolation_main(data,loggerInfo,"basis",basis_method = "svd", 
																					simNum = 100, intMethod = "baye", r = 10)
# method : loglik, baye, IDW
plot(basis_hypoxiaExtent)
	# interpolation_controller(data,loggerInfo) %>% 
	# plot_spatial(locationInfo = loggerInfo,outputFolder = "../output/2015spPlots/")
# }







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
