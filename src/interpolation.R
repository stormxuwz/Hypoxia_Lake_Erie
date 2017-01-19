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



idw_interpolation_main <- function(data, locationInfo,...){
	data <- na.omit(data)  # only remain the time where all data are available
	times <- index(data)

	grid <- list(...)$grid

	if(is.null(grid)){
		grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)  # will also return an area
	}
	
	
	prediction <- array(0, dim = c(1, length(times), nrow(grid)))
	
	loggerNames  <- as.numeric(colnames(data))
	print("sampling data number:")
	print(nrow(data))
	
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


crossValidation <- function(data, locationInfo, basis, rList = c(5:12),method = "IDW"){
	
	crossPred <- function(r,j){
		# r is the basis, j is the left out logger index
		trainData <- data[,-j]
		testData <- data[,j]
		trainLocation <- locationInfo[-j,]
		targetLocation <- locationInfo[j,c("longitude","latitude","x","y")]
		targetLocation$convexIndex <- 1
		
		res <- basis_interpolation_step1(data, trainLocation, basis, simNum=100, intMethod="baye", r=r, residualMethod = "IDW",saveMeta = FALSE, grid = targetLocation)
		pred <- basis_interpolation_step2(res[[1]],res[[2]],nSim = 100,TRUE,FALSE)
		
		nsim <- length(pred) # total simulations

		pred  <- data.frame(lapply(pred,function(x) return(x[,1]))) # extract the first column which should only have the one column
		colnames(pred) <- paste0("sim_",as.character(1:nsim))

		return(list(trend = pred, pred_res = res[[2]][,,1],label = as.numeric(testData)))
	}

	for(r in rList){
		for(j in 1:nrow(locationInfo)){
			cp <- crossPred(r,j)
			fileName = sprintf("%scv_r%d_%s.rds",metaFolder,r,locationInfo[j,"loggerID"])
			saveRDS(cp,fileName)
			crossValidation_summary(fileName,withResidual = TRUE)
		}
	}
}


crossValidation_summary <- function(savedFile,withResidual= TRUE){
	d <- readRDS(savedFile)
	info <- strsplit(strsplit(basename(savedFile),"[.]")[[1]][1],"_")[[1]]
	nsim = length(d$trend)
	
	if(withResidual){
		for(i in 1:nsim){
			d$trend[,i] = d$trend[,i]+d$pred_res
		}
		m = rowMeans(d$trend)
		std = apply(d$trend, 1, var)
	}

	plotDF <- data.frame(m = m, upper = m+std, lower = m-std,true = d$label,time = 1:length(m))

	yRange <- range(m,m+std,m-std,d$label)

	p = ggplot(data = plotDF) + geom_point(aes(time,m),color = "red")+geom_point(aes(time,true),color = "black")+geom_point(aes(time,upper),color = "blue")+geom_point(aes(time,lower),color = "blue")+ylim(yRange)+ylab("DO")

	png(sprintf("%splot_%s.png",metaFolder,paste(info,collapse = "_")))
	print(p)
	dev.off()
}



basis_interpolation_step1 <- function(data,locationInfo, basis, simNum, intMethod, r, residualMethod = "IDW",saveMeta = TRUE,...)
{
	# function to predict the coefficients at every grid point
	# Save the coefficient predictions

	if(r>nrow(locationInfo)){
		throw("Too large r, please select r<= # of location")
	}
	
	if(is.null(list(...)$grid))
		grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)
	else{
		grid <- list(...)$grid
	}

	trendBasisCoeff <- basis_interpolation(data,locationInfo,grid,basis = basis,simNum = simNum, intMethod = intMethod, r = r)
	
	# save the intermediate results, the coefficients of the basis of the trend
	if(saveMeta)
		saveRDS(list(basisIntRes = trendBasisCoeff, grid = grid),paste(metaFolder,"trend_basis.rds",sep=""))
	
	# Interpolation on residuals
	if(residualMethod == "IDW"){
		# do IDW on the residuals
		print("doing IDW on residuals")
		residual_interp <- idw_interpolation_main(trendBasisCoeff$res, locationInfo,...)
	}
	else if(residualMethod == "stKrig"){
		# do stkriging on the residuals
		print("doing spatial temporal kriging on residuals")
		residual_interp <- st_interpolation_main(trendBasisCoeff$res, locationInfo,...)
	}
	else{ # ignore the residuals
		residual_interp <- array(0, dim = c(1, length(index(data)), nrow(grid)))
	}

	if(saveMeta) # save the residuals 
		saveRDS(residual_interp,paste(metaFolder,"residual_prediction_basis.rds",sep =""))

	return(list(trendBasisCoeff=trendBasisCoeff,residual_interp=residual_interp))
}	

basis_interpolation_step2 <- function(trendBasisCoeff,residual_interp,nSim,parallel,returnHypoxia){
	#trend_basis <- readRDS(paste0(metaFolder,"trend_basis.rds"))
	#residual_interp <- readRDS(paste0(metaFolder,"residual_prediction_basis.rds"))

	coeffPredList <- trendBasisCoeff$trendCoeff  # list of [#grid,# of sim of coefficients]
	basis <- trendBasisCoeff$basis
	
	inds <- sample(1:dim(coeffPredList[[1]])[2],nSim*ncol(basis),replace= TRUE)

	inds <- matrix(inds,nrow = nSim, ncol = ncol(basis))


	if(parallel){
		print("Parallel calculating simulations")

		require(doParallel)
		cl <- makeCluster(6)
		registerDoParallel(cl)
		
		res <- foreach(sim = 1:nSim) %dopar% {  # res is a list of the results of each iteration
			prediction = 0 

			for(i in 1:ncol(basis)){
				pred <- coeffPredList[[i]]
				ind <- inds[sim,i] # take a sample of the coefficients interpolation grid
				prediction <- prediction+basis[,i] %*% t(pred[,ind]) # get the time series of each grid
			}
			
			prediction <- prediction + residual_interp[1,,]  # for deterministic interpolation on the residuals, choose [1,,]
			
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
		if(returnHypoxia)
			res <- data.frame(res)
	}
	return(res)
}

basis_interpolation_main <- function(data, locationInfo, basis, simNum, intMethod, r, timeIndex,residualMethod = "IDW"){
	# simNum: the number of simulations of spatial basis coefficients 
	# basis: how to decompose the data
	# intMethod: interpolation method
	# r: 

	res <- basis_interpolation_step1(data, locationInfo, basis, simNum, intMethod, r, residualMethod = "None")
	return(basis_interpolation_step2(res[[1]],res[[2]],nSim = 1000,parallel = TRUE,returnHypoxia = FALSE))
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
