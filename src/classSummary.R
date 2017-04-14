summary.idwModel <- function(model){
	# return the hypoxia extent in the idw model
	DOvalue <- reConstruct(model)
	return(hypoxiaCount(DOvalue))
}


summary.basisModel <- function(trendModel, residualPrediction, parallel, totalSim, indMatrix = NULL){
	# return the hypoxia extent in the basis model
	availableSimNum <- ncol(trendModel$predictions[[1]]$simulations)
	
	grid <- trendModel$grid
	nTime <- nrow(residualPrediction)

	if(is.null(indMatrix)){
		print("Create new sampling index matrix")
		indMatrix <- base::sample(1:availableSimNum, totalSim*basisNum,replace= TRUE) %>%
				matrix(nrow = basisNum)
	}else{
		stopifnot(totalSim <= ncol(indMatrix))
		indMatrix <- indMatrix[, 1:totalSim]
	}

	basisNum <- length(trendModel$predictions)
	if(parallel){
		require(doParallel)
		cl <- makeCluster(6)
		registerDoParallel(cl)		

		res <- foreach(simIdx = 1:totalSim) %dopar% {
			source("src/classReConstruct.R")
			source("src/classSummary.R")
			DOvalue <- reConstruct(trendModel, residualPrediction, simIdx, indMatrix)$predValue
			hypoxiaCount(DOvalue)
		}

		stopCluster(cl)
		gc()
	}
	# res is a list of data frame, each data frame containing count of 0,2,4
	hypoxiaExtent_0 <- matrix(0, nTime, totalSim)
	hypoxiaExtent_2 <- hypoxiaExtent_0
	hypoxiaExtent_4 <- hypoxiaExtent_0

	for(i in 1:totalSim){
		hypoxiaExtent_0[,i] <- res[[i]]$less0
		hypoxiaExtent_2[,i] <- res[[i]]$less2
		hypoxiaExtent_4[,i] <- res[[i]]$less4
	}

	hypoxiaSummary <- lapply(list(hypoxiaExtent_0,hypoxiaExtent_2,hypoxiaExtent_4),summaryStatistics)
	names(hypoxiaSummary) <- c("less0","less2","less4")
	return(hypoxiaSummary)
}

summaryStatistics <- function(x){
	return(data.frame(
		median = apply(x,1, median),
		upper = apply(x,1, quantile, probs = c(0.95), na.rm = TRUE),
		lower = apply(x,1, quantile, probs = c(0.05), na.rm = TRUE)))
}


hypoxiaCount <- function(prediction){
	prediction <- prediction * (prediction>0)

	hypoxiaExtent_0 <- rowSums(prediction<0.01,na.rm =TRUE)
	hypoxiaExtent_2 <- rowSums(prediction<2,na.rm =TRUE)
	hypoxiaExtent_4 <- rowSums(prediction<4,na.rm =TRUE)
	res <- data.frame(less0 = hypoxiaExtent_0,
			less2 = hypoxiaExtent_2,
			less4 = hypoxiaExtent_4)
	return(res)
}


cvUncertainty <- function(prediction, trueValue){
	nsim <- length(prediction$predValue)
	nTime <- length(prediction$predValue[[1]])

	allPrediction <- matrix(NA, nrow = nTime, ncol = nsim)

	for(i in 1:nsim){
		allPrediction[,i] <- prediction$predValue[[i]] * (prediction$predValue[[i]] > 0)
	}

	
	tmp <- summaryStatistics(allPrediction)
	tmp$true <- as.numeric(trueValue)
	tmp$withinBound <- (tmp$true < tmp$upper) & (tmp$true > tmp$lower)
	attr(tmp, "boundRatio") <- sum(tmp$withinBound)/nrow(tmp)
	return(tmp)
}


