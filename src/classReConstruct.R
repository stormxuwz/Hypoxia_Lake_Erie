reConstruct <- function(x,...){
	UseMethod("reConstruct")
}


reConstruct.idwModel <- function(model){
	# res is a list containing predictions (list of pred) and grid

	prediction <- matrix(NA, nrow = length(model$predictions), ncol = nrow(model$grid))

	for(i in 1:length(model$predictions)){
		prediction[i,] <- model$predictions[[i]]$pred[,1]
	}
	return(prediction)
}


reConstruct.basisModel <- function(
	trendModel,
	residualPrediction, 
	simulationNum = -1,
	indMatrix = NULL,
	parallel = FALSE){
	# trendModel: list containing 
		# predictions (list of pred, simulations)
		# grid
		# basis
	# residualPrediction is a T * n_grid matrix
	# simulationNum 

	# return a list of 1 matrix:T * n_grid and 2 variance n_grid

	availableSimNum <- ncol(trendModel$predictions[[1]]$simulations)
	basisNum <- length(trendModel$predictions)
	basis <- trendModel$basis

	prediction = 0
	variance = 0
	totalSim <- abs(simulationNum)

	if(is.null(indMatrix)){
		print("Create new sampling index matrix")
		indMatrix <- sample(1:availableSimNum, totalSim*basisNum,replace= TRUE) %>%
				matrix(nrow = basisNum)
	}
	
	stopifnot(totalSim <= ncol(indMatrix)) # check totalSim is feasible

	if(simulationNum == 0){
		# reconstruct the BLUE estimates
		for(i in 1:basisNum){
			coeff <- trendModel$predictions[[i]]$pred[,1]
			prediction <- prediction + basis[,i] %*% t(coeff)
			variance <- variance + trendModel$predictions[[i]]$pred[,2]
		}
		prediction <- prediction + residualPrediction

	}else if(simulationNum < 0){
		# reConstruct all simulations
		if(parallel){
			require(doParallel)
			cl <- makeCluster(6)
			registerDoParallel(cl)		

			prediction <- foreach(simIdx = 1:totalSim) %dopar% {
					source("src/classReConstruct.R")
					source("src/classSummary.R")
					reConstruct(trendModel, residualPrediction, simIdx, indMatrix)$predValue
					# T * n_grid and 2 variance n_grid
			}
			variance <- reConstruct(trendModel, residualPrediction, 1, indMatrix)$predVariance
			stopCluster(cl)
			gc()		
		}
	}else{
		# generate from each basis predictions
		ind <- indMatrix[, simulationNum]
		for(i in 1:basisNum){
			coeff <- trendModel$predictions[[i]]$simulations[, ind[i]]
			prediction <- prediction + basis[,i] %*% t(coeff)
			variance <- variance + trendModel$predictions[[i]]$pred[,2]
		}
		prediction <- prediction + residualPrediction
	}

	return(list(predValue = prediction, predVariance = variance))
}
