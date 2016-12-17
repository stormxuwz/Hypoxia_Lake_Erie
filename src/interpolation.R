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
	
	require(doParallel)
	cl <- makeCluster(6)
	registerDoParallel(cl)
	
	res <- foreach(i=1:nrow(data)) %dopar% { 
		library(dplyr)		
		source("src/spatialHelper.R")
		source("src/basisDecomposition.R")
		source("src/interpolation_base.R")
		subData <- data.frame(logger = loggerNames,DO = as.numeric(data[i,]))
		subData <- merge(subData, locationInfo,by.x = "logger",by.y = "loggerID") %>% rename(value = DO)
		pred <- spatial_interpolation(subData,grid)[,1]
		pred <- ifelse(pred<0,0,pred)
		pred
	}
	stopCluster(cl)
	prediction[1,,] <- t(matrix(unlist(res),nrow = nrow(grid)))
	gc()


	# for(i in 1:nrow(data)){
		# subData <- data.frame(logger = loggerNames,DO = as.numeric(data[i,]))
		# subData <- merge(subData, locationInfo,by.x = "logger",by.y = "loggerID") %>% rename(value = DO)
		
		# pred <- spatial_interpolation(subData,grid)[,1]
		# prediction[1,i,] <- ifelse(pred<0,0,pred)
	# }
	
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
	
	crossPred <- function(r,j){
		trainData <- data[,-j]
			testData <- data[,j]
			trainLocation <- locationInfo[-j,]
			targetLocation <- locationInfo[j,c("longitude","latitude","x","y")]
			targetLocation$convexIndex <- 1
			
			grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)

			res <- basis_interpolation_step1(data, locationInfo, basis, simNum=100, intMethod="baye", r=r, residualMethod = "IDW",FALSE)

			pred <- basis_interpolation_step2(r[[1]],r[[2]],nSim = 100,TRUE,FALSE)
			return(data.frame(pred = pred, pred_res =r[[2]], label = as.numeric(testData)))
	}

	for(r in 1:length(rList)){
		for(j in 1:nrow(locationInfo)){
			cp <- crossPred(r,j)
			saveRDS(cp,sprintf("%s_cv_r%d_%s.rds",metaFolder,r,locationInfo[j,"loggerID"]))
		}
	}
}




basis_interpolation_step1 <- function(data, locationInfo, basis, simNum, intMethod, r, residualMethod = "IDW",saveMeta = TRUE){
	
	if(r>nrow(locationInfo)){
		throw("Too large r, please select r<= # of location")
	}
	
	grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)
	trendBasisCoeff <- basis_interpolation(data,locationInfo,grid,basis = basis,simNum = simNum, intMethod = intMethod, r = r)
	
	# save the intermediate results
	if(saveMeta)
		saveRDS(list(basisIntRes = trendBasisCoeff, grid = grid),paste(metaFolder,"trend_basis.rds",sep=""))
	

	# Interpolation on residuals
	if(residualMethod == "IDW"){
		# do IDW on the residuals
		print("doing IDW on residuals")
		residual_interp <- idw_interpolation_main(trendBasisCoeff$res, locationInfo)
	}
	else if(residualMethod == "stKrig"){
		# do stkriging on the residuals
		print("doing spatial temporal kriging on residuals")
		residual_interp <- st_interpolation_main(trendBasisCoeff$res, locationInfo)
	}
	else{ # ignore the residuals
		residuals_interp <- 0
	}
	if(saveMeta)
		saveRDS(residual_interp,paste(metaFolder,"residual_prediction_basis.rds",sep =""))

	return(list(trendBasisCoeff=trendBasisCoeff,residual_interp=residual_interp))
}	

basis_interpolation_step2 <- function(trendBasisCoeff,residual_interp,nSim,parallel,returnHypoxia){
	# trend_basis <- readRDS(paste0(metaFolder,"trend_basis.rds"))
	# residual_interp <- readRDS(paste0(metaFolder,"residual_prediction_basis.rds"))

	coeffPredList <- trendBasisCoeff$trendCoeff
	basis <- trendBasisCoeff$basis
	
	if(parallel){
		require(doParallel)
		cl <- makeCluster(6)
		registerDoParallel(cl)
		
		res <- foreach(sim = 1:nSim) %dopar% {
			prediction = 0
			for(i in 1:ncol(basis)){
				pred <- coeffPredList[[i]]
				ind <- sample(1:dim(pred)[2],1,replace= TRUE)
				prediction <- prediction+basis[,i] %*% t(pred[,ind])
			}
			prediction <- prediction + residual_interp[1,,]
			
			if(returnHypoxia){
				hypoxiaExtent_0 <- rowSums(prediction<0.01,na.rm =TRUE)
				hypoxiaExtent_2 <- rowSums(prediction<2,na.rm =TRUE)
				hypoxiaExtent_4 <- rowSums(prediction<4,na.rm =TRUE)
				c(hypoxiaExtent_0,hypoxiaExtent_2,hypoxiaExtent_4)	
			}
			else{
				prediction
			}
		}
		stopCluster(cl)
		res <- data.frame(res)
	}
	return(res)
}

basis_interpolation_main <- function(data, locationInfo, basis, simNum, intMethod, r, timeIndex,residualMethod = "IDW"){
	res <- basis_interpolation_step1(data, locationInfo, basis, simNum, intMethod, r, residualMethod = "IDW")
	return(basis_interpolation_step2(res[[1]],res[[2]],nSim = 1000,TRUE,TRUE))
}


interpolation_main <- function(data,locationInfo,method = "IDW",...){
	timeIndex <- index(data)
	args <- list(...)

	if(method == "IDW"){
		prediction <- idw_interpolation_main(data,locationInfo)
		saveRDS(prediction,paste(outputFolder,"final_prediction_IDW.rds",sep =""))
		hypoxiaExtent <- summaryHypoxia(prediction, timeIndex)
	}
	else if(method == "basis"){
		hypoxiaExtent <- basis_interpolation_main(data,locationInfo, 
			basis = args$basis_method, 
			simNum = args$simNum, 
			intMethod = args$intMethod, 
			r = args$r,timeIndex=timeIndex)
	}else{
		throw("not implemented")
	}
	
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
