# calculate the hypoxia days

hypoxiaAnalysis <- function(DO_predictions, DOthreshold){
	# Return the longest hypoxia of a point
	# Input:
	# 	DO_predictions: matrix of shape T * N
	#	DOthreshold: the threshold of DO to define hypoxia
	# return :
	# 	A list containing
	# 	hypoxiaTime: the hypoxia time of each point
	# 	longestTime: the longest continuous hypoxia time
	hypoxiaMatrix <- DO_predictions < DOthreshold

	hypoxiaTime <- colSums(hypoxiaMatrix)

	longestTime <- rep(0, ncol(DO_predictions))

	for(i in 1:ncol(DO_predictions)){
		tmp <- rle(hypoxiaMatrix[,i])
		if(!TRUE %in% tmp$values){
			longestTime[i] <- 0
		}else{
			longestTime[i] <- max(tmp$lengths[tmp$values == TRUE])
		}
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

	return(list(hypoxiaFit = fit, featureImp = model))
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


clusterByNMF <- function(lakeDOObj, basis_r){
	# cluster sensor raw data by decomposing by NMF and cluster based on NMP coefficients
	# Input:
	# 	lakeDOObj: lakeDO object
	# 	basis_r: the rank of the lower matrix
	# 	numCluster: how many cluster to make

	require(ggplot2)
	nmfDecomp <- lakeDOObj$samplingData %>% na.omit() %>% NMF_basis(basis_r)
	coeff <- nmfDecomp$coef %>% t() %>% scale()
	print(dim(coeff))

	clusterRes <- list()

	for(nc in c(2,3,4,5,6)){
		if(nc >= basis_r){
			next
		}
		clusterRes[[paste0("cluster_", nc) ] ] <- kmeans(coeff, nc)$cluster
	}
	

	# lakeDOObj$loggerInfo$cluster <- as.factor(clusterRes$cluster)
	# baseMap <- readRDS("erieGoogleMap_2014.rds")
	# baseMap + geom_point(aes(longitude, latitude, color = cluster), data = lakeDOObj$loggerInfo)

	return(list(basis = nmfDecomp$basis, cluster = clusterRes))
}


resultSummary <- function(){
	# summary rmse and withinBoundRatio of each logger of each method in CV
	fullRes <- data.frame()
	for(year in c(2014, 2015)){
		for(aggType in c("daily","hourly")){
			erieDO <- getLakeDO(year, "B", aggType)
			erieDO$samplingData <- na.omit(erieDO$samplingData)
			for (cv_loggerID in erieDO$loggerInfo$loggerID){
				for(method  in c("idw","Reml","Baye")){
					for(r in c(5,10,15)){
						if(method == "idw" & r>5){
							next
						}else{
							metaList <- list(year = year, method = method, aggType = aggType, cv_loggerID = cv_loggerID, r = r)
							res <- c(metaList,plot_cv(year, method, aggType, cv_loggerID, r)) %>% data.frame()
							fullRes <- rbind(fullRes, res)
						}
					}
				}
			}
		}
	}

	saveRDS(fullRes,sprintf("%s/results/fullRes.rds",outputBaseName))
	return(fullRes)
}

plotCVOnMap <- function(year_, aggType_, method_, r_){
	# function to plot CV rmse and withinBoundRatio on the map
	library(ggmap)

	if(method_ == "idw") r_ = 5
	erieDO <- getLakeDO(year_, "B", aggType_)
	res <- readRDS(sprintf("%s/results/fullRes.rds",outputBaseName)) %>% 
	subset(year == year_ & aggType == aggType_ & method == method_ & r == r_) %>% 
	merge(erieDO$loggerInfo, by.x = "cv_loggerID", by.y = "loggerID", all.x = TRUE)

	# lonRange <- range(res$longitude)
	# latRange <- range(res$latitude)
	# bbox <- make_bbox(lonRange,latRange,f = 0.2)
	# myMap <- get_map(location=bbox, source="google",crop=FALSE) %>% ggmap()
	# saveRDS(myMap,"erieGoogleMap.rds")

	myMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",year_))
	if(method_ == "idw"){
		p <- myMap + geom_point(aes(longitude, latitude, color = rmse), data = res, size = 5)
	}else{
		p <- myMap + geom_point(aes(longitude, latitude, size = withinBoundRatio, color = rmse), data = res)
	}
	p <- p + scale_color_gradientn(colors = terrain.colors(10))
	return(p)
	#p_rmse <-  myMap + geom_point(aes(longitude, latitude, color = rmse), data = res)
	#p_inBound <- myMap + geom_point(aes(longitude, latitude, color = withinBoundRatio), data = res)
}


RMSE <- function(x,y){
	return(round(mean((x-y)^2),digits = 3))
}

plot_cv <- function(year, method, aggType, cv_loggerID, r = NULL){
  # function to plot time seris of CV resutls and return rmse, withinBoundRatio
	if(method == "idw"){
		metaFolder <- sprintf("%s%d_%s_idw/cv/%s/",outputBaseName, year, aggType, cv_loggerID)
		figName <- sprintf("%s/results/%d/%s_%s_%s.pdf", outputBaseName, year, cv_loggerID, method, aggType)

	}else{
		metaFolder <- sprintf("%s%d_%s_%s_%d/cv/%s/",outputBaseName, year, aggType, method, r, cv_loggerID)
		figName <- sprintf("%s/results/%d/%s_%s_%s_r%d.pdf", outputBaseName, year, cv_loggerID,method, aggType, r)
	}

	df <- readRDS(sprintf("%s/simulations_stat.rds",metaFolder))
	t <- index(df)
	
	df <- as.data.frame(df)
	df$time <- t

	withinBoundRatio <- sum(df$withinBound)/nrow(df)
	rmse <- RMSE(df$median,df$true)

	pdf(figName)
	p <- ggplot(data = df) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,median),color = "black") + 
		geom_line(aes(time,true),color = "red",alpha = 0.8) + 
		ylim(c(-1,20))+ylab("DO (mg/L)")+ xlab("Time") + 
		ggtitle(sprintf("Logger: %s, RMSE:%.3f, inBoundRatio: %.3f",cv_loggerID, rmse, withinBoundRatio))

	print(p)
	dev.off()
	return(list(rmse = rmse, withinBoundRatio = withinBoundRatio))
}	



