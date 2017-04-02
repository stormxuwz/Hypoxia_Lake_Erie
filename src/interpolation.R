
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


idw_interpolation_main <- function(data, locationInfo, grid = NULL){
	data <- na.omit(data)  # only remain the time where all data are available
	times <- index(data)

	if(is.null(grid)){
		grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)  
	}
	
	prediction <- list()
	loggerNames  <- as.numeric(colnames(data))
	# print("sampling data number:")
	# print(nrow(data))
	
	require(doParallel)
	
	cl <- makeCluster(6)
	registerDoParallel(cl)
	
	print("Parallel doing IDW interpolation")
	res <- foreach(i=1:nrow(data)) %dopar% { 
		library(dplyr)		
		source("src/spatialHelper.R")
		source("src/basisDecomposition.R")
		source("src/interpolation_base.R")
		subData <- data.frame(logger = loggerNames,DO = as.numeric(data[i,]))
		subData <- merge(subData, locationInfo,by.x = "logger",by.y = "loggerID") %>% rename(value = DO)
		pred <- spatial_interpolation(subData,grid, method="IDW")[,1]
		pred <- ifelse(pred<0,0,pred)
		pred
	}
	stopCluster(cl)
	prediction[[1]] <- t(matrix(unlist(res),nrow = nrow(grid)))
	gc()

	return(prediction)
}

st_interpolation_main <- function(data, locationInfo,grid = NULL,...){
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


basis_interpolation_step1 <- function(data,locationInfo, basisDecomp, simNum, fitMethod, r, residualMethod, grid = NULL, saveMeta = TRUE)
{
	# function to predict the coefficients at every grid point
	# Save the coefficient predictions

	if(r > nrow(locationInfo)){
		throw("Too large r, please select r<= # of location")
	}
	
	if(is.null(grid))
		grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)

	trendBasisCoeff <- basis_interpolation(
		DOdata = data,
		logger_geo = locationInfo,
		grid = grid, 
		basisDecomp = basisDecomp, 
		simNum = simNum, 
		fitMethod = fitMethod, 
		r = r)
	
	# save the intermediate results, the coefficients of the basis of the trend
	if(saveMeta)
		saveRDS(list(basisIntRes = trendBasisCoeff, grid = grid),paste0(metaFolder,"trend.rds"))
	
	# Interpolation on residuals
	if(residualMethod == "IDW"){
		# do IDW on the residuals
		print("doing IDW on residuals")
		residual_interp <- idw_interpolation_main(trendBasisCoeff$res, locationInfo, grid = grid)
	}
	else if(residualMethod == "stKrig"){
		# do stkriging on the residuals
		print("doing spatial temporal kriging on residuals")
		residual_interp <- st_interpolation_main(trendBasisCoeff$res, locationInfo, grid = grid)
	}
	else{ # ignore the residuals
		residual_interp <- list()
		residual_interp[[1]] <- matrix(0, nrow = length(index(data)), ncol = nrow(grid))
	}

	if(saveMeta) # save the residuals 
		saveRDS(residual_interp,paste0(metaFolder,"residual_prediction.rds"))

	return(list(trendBasisCoeff=trendBasisCoeff,residual_interp=residual_interp))
}	

basis_interpolation_step2 <- function(trendBasisCoeff,residual_interp,nSim,parallel,returnHypoxia, saveMeta = TRUE){
	# reconstruct from the coefficient interpolation
	# trend_basis <- readRDS(paste0(metaFolder,"trend_basis.rds"))
	# residual_interp <- readRDS(paste0(metaFolder,"residual_prediction_basis.rds"))

	coeffPredList <- trendBasisCoeff$trendCoeff  # list of [#grid, nSim]
	basis <- trendBasisCoeff$basis
	
	inds <- sample(1:dim(coeffPredList[[1]])[2],nSim*ncol(basis),replace= TRUE)
	inds <- matrix(inds,nrow = nSim, ncol = ncol(basis))

	if(saveMeta)
		saveRDS(inds,paste0(metaFolder,"inds.rds"))

	if(parallel){
		print(sprintf("Parallel calculating %d simulations",nSim))
		require(doParallel)
		cl <- makeCluster(6)
		registerDoParallel(cl)
		
		res <- foreach(sim = 1:nSim) %dopar% {  # res is a list of the results of each iteration
			prediction = 0 
			for(i in 1:ncol(basis)){
				# get the time series of each grid
				prediction <- prediction + 
					basis[,i] %*% t(coeffPredList[[i]][,inds[sim,i]]) 
			}
			# add residual interpolation part
			# for deterministic interpolation on the residuals, choose [1]
			prediction <- prediction + residual_interp[[1]]  
			prediction <- prediction*(prediction>0)

			if(returnHypoxia){
				hypoxiaExtent_0 <- rowSums(prediction<0.01,na.rm =TRUE)
				hypoxiaExtent_2 <- rowSums(prediction<2,na.rm =TRUE)
				hypoxiaExtent_4 <- rowSums(prediction<4,na.rm =TRUE)
				data.frame(less0 = hypoxiaExtent_0,
					less2 = hypoxiaExtent_2,
					less4 = hypoxiaExtent_4)	
			}
			else{
				prediction
			}
		}
		stopCluster(cl)
	}
	return(res)
}

basis_interpolation_main <- function(data, locationInfo, basisDecomp, simNum, fitMethod, r, residualMethod, grid = NULL){
	# simNum: the number of simulations for spatial basis coefficients 
	# basis: how to decompose the data
	# fitMethod: interpolation method
	# r: 

	res <- basis_interpolation_step1(data, locationInfo, basisDecomp, simNum, fitMethod, r, residualMethod = residualMethod, grid = grid)
	return(basis_interpolation_step2(res[[1]],res[[2]],
		nSim = 1000, 
		parallel = TRUE,
		returnHypoxia = TRUE))
}


interpolation_main <- function(data,locationInfo,method, grid, ...){
	timeIndex <- index(data)
	args <- list(...)
	
	if(method == "IDW"){
		hypoxiaExtent <- idw_interpolation_main(data,locationInfo,grid = grid)
	}
	else if(method == "basis"){
		hypoxiaExtent <- basis_interpolation_main(
			data = data,
			locationInfo = locationInfo, 
			basisDecomp = args$basisDecomp, 
			simNum = args$simNum, 
			fitMethod = args$fitMethod, 
			r = args$r,
			residualMethod = args$residualMethod,
			grid = grid)
	}else{
		throw("not implemented")
	}
	
	return(hypoxiaExtent)
}

