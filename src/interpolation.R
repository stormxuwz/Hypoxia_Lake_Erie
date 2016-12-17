# main api for different interpolation method
# rm(list = ls())
require(fields)
require(gstat)
require(sp)
require(dplyr)
require(geoR)

source("src/spatialHelper.R")
source("src/basisDecomposition.R")
source("src/interpolation_base.R")
source("src/interp_summary.R")


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

st_interpolation_main <- function(data, locationInfo,grid = NULL){
	stv <- st_variogram(data,locationInfo)
	saveRDS(stv,file = paste(metaFolder,"st_variogram.rds",sep = ""))
	
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


crossValidation <- function(data, locationInfo, basis, rList = c(5:12)){
	
	crossPred <- function(i,j){
		trainData <- data[,-j]
			testData <- data[,j]
			trainLocation <- locationInfo[-j,]
			targetLocation <- locationInfo[j,c("longitude","latitude","x","y")]
			targetLocation$convexIndex <- 1
			trendPrediction <- basis_interpolation(trainData,trainLocation, targetLocation, basis = basis, simNum = 1,intMethod = "loglik", r =rList[[i]], saveMeta = FALSE)
			pred <- trendPrediction$trend[1,,]
			residual_interp <- idw_interpolation_main(trendPrediction$res, locationInfo)
			# residual_interp  <- 0
			return(data.frame(pred = pred+residual_interp, pred_trend =trendPrediction, label = as.numeric(testData)))
	}

	for(i in 1:length(rList)){
		for(j in 1:nrow(locationInfo)){
			cp <- crossPred(i,j)
			saveRDS(cp,sprintf("%s_cv_r%d_%s.rds",metaFolder,i,locationInfo[j,"loggerID"]))
		}
	}
}




basis_interpolation_main <- function(data, locationInfo, basis, simNum, intMethod, r, residualMethod = "IDW"){
	
	if(r>nrow(locationInfo)){
		throw("Too large r, please select r<= # of location")
	}
	
	grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)
	trendPrediction <- basis_interpolation(data,locationInfo,grid,basis = basis,simNum = simNum, intMethod = intMethod, r = r)
	
	# save the intermediate results
	saveRDS(list(basisIntRes = trendPrediction, grid = grid),paste(metaFolder,"trend_prediction_basis.rds",sep=""))
	
	trendSim <- trendPrediction$trend 	# trend simulations, shape of (simNum, TimeN, nrow(grid)
	varExpl <- trendPrediction$varExpl  # variance explained
	
	print("doing spatial temporal kriging on residuals")
	
	if(residualMethod == "IDW"){
		# do IDW on the residuals
		residual_interp <- idw_interpolation_main(trendPrediction$res, locationInfo)
	}
	else if(residualMethod == "stKrig"){
		# do stkriging on the residuals
		residual_interp <- st_interpolation_main(trendPrediction$res, locationInfo)
	}
	else{ # ignore the residuals
		residuals_interp <- 0
	}
	
	saveRDS(residual_interp,paste(metaFolder,"residual_prediction_basis.rds",sep =""))
	
	for(i in 1:simNum){
		trendSim[i,,] <- trendSim[i,,] + residual_interp[1,,]
	}
	# trendSim[i,,] <- ifelse(trendSim[i,,]>0, trendSim[i,,], 0)

	saveRDS(trendSim,paste(outputFolder,"final_prediction_basis.rds",sep =""))
	return(trendSim)
}	


interpolation_main <- function(data,locationInfo,method = "IDW",...){
	timeIndex <- index(data)
	args <- list(...)

	if(method == "IDW"){
		prediction <- idw_interpolation_main(data,locationInfo)
		saveRDS(prediction,paste(outputFolder,"final_prediction_IDW.rds",sep =""))
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
