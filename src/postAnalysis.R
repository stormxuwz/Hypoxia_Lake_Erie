# calculate the hypoxia days

hypoxiaAnalysis <- function(DO_predictions, DOthreshold){
	# Input:
	# 	DO_predictions: matrix of shape T * N
	#	DOthreshold: the threshold of DO to define hypoxia
	# return :
	# 	A list containing
	# 	hypoxiaTime: the hypoxia time of each point
	# 	longestTime: the longtest continuous hypoxia time
	hypoxiaMatrix <- DO_predictions < DOthreshold

	hypoxiaTime <- colSums(hypoxiaMatrix)

	longestTime <- rep(0, ncol(DO_predictions))

	for(i in 1:ncol(DO_predictions)){
		tmp <- rle(hypoxiaMatrix[,i])
		longestTime[i] <- max(tmp$lengths[tmp$values == "TRUE"])
	}

	return(list(hypoxiaTime = hypoxiaTime, longestTime = longestTime))
}


identifyKeyLoggersToHypoxia <- function(lakeDOObj, hypoxiaExtentStream, DOthreshold){
	# quick way to estimate hypoxia extent
	# to identitfy which loggers plays a key role in determining the area of hypoxia
	# so that worth to put sensor on 
	loggerNames <- names(lakeDOObj$samplingData)

	# using whether loggers are less than DO threshold
	hypoxiaMatrix <- as.data.frame(lakeDOObj$samplingData) < DOthreshold
	# or using the pure values (which needs a non-linear model to do)
	hypoxiaMatrix$Y <- hypoxiaExtentStream

	model <- lm(Y~., data = hypoxiaMatrix)

	fit <- predict(model, hypoxiaMatrix)

	return(list(hypoxiaFit = fit, featureImp <- model))
}

getCoeffMatrix <- function(krigModelObj){
	basisNum <- length(krigModelObj)
	numLogger <- nrow(krigModelObj[[1]])
	basisCoeffMatrix <- matrix(0, nrow=numLogger, ncol=basisNum)

	for(i in 1:basisNum){
		basisCoeffMatrix[,i] <- krigModelObj[[i]]$value
	}
	return(basisCoeffMatrix)
}

cluster <- function(lakeDOObj, basis_r, numCluster){
	require(ggplot2)
	nmfDecomp <- lakeDOObj$samplingData %>% na.omit() %>% NMF_basis(basis_r)
	coeff <- nmfDecomp$coef %>% t() %>% scale()
	clusterRes <- kmeans(coeff, numCluster)

	# lakeDOObj$loggerInfo$cluster <- as.factor(clusterRes$cluster)
	# baseMap <- readRDS("erieGoogleMap_2014.rds")
	# baseMap + geom_point(aes(longitude, latitude, color = cluster), data = lakeDOObj$loggerInfo)

	return(list(basis = nmfDecomp$basis, cluster = clusterRes$cluster))
}



