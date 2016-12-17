
summarySimulation <- function(hypoxiaExtentPred,summaryName){
	
	std <- sqrt(apply(hypoxiaExtentPred, 2, var))
	hypoxiaSummary <- data.frame(m = colMeans(hypoxiaExtentPred))
	hypoxiaSummary$upper <- hypoxiaSummary$m+2*std
	hypoxiaSummary$lower <- hypoxiaSummary$m-2*std
	
	names(hypoxiaSummary) <- c(paste(summaryName,c("Mean","Upper","Lower"),sep=""))
	return(hypoxiaSummary)
}

summaryHypoxia <- function(predictionMatrix, timeIndex, raw = FALSE){
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
	if(raw){
		return(list(hypoxiaExtent_0,hypoxiaExtent_2,hypoxiaExtent_4))
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

