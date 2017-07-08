
getHypoxiaExtent <- function(year, aggType, method, r){
	# function to plot hypoxia time series
	IDWFolderName <- sprintf("../output/%d_%s_idw/",year, aggType)
	
	folderName <- sprintf("../output/%d_%s_%s_%d/",year, aggType, method, r)
	
	model <- readRDS(paste0(folderName,"basisModelRes.rds"))
	trendModel <- readRDS(paste0(folderName, "trendModel.rds"))
	grid <- trendModel$grid
	n <- sum(grid$convexIndex)
	timeIndex <- index(model$residuals$samplingData)
	HEList <- readRDS(paste0(folderName,"extent.rds"))
	
	idwHE <- readRDS(paste0(IDWFolderName,"extent.rds"))
	
	units <- list(less0 = "0.01 mg/L", less2 = "2 mg/L", less4 = "4 mg/L")
	
	n <- sum(grid$convexIndex)
	
	for(threshold in c("less0","less2","less4")){
		HE <- HEList[[threshold]]
		HE$idw <- idwHE[[threshold]]
		HE$time <- timeIndex
		
		p <- ggplot(HE) + geom_ribbon(aes(time,ymin = lower/n, ymax = upper/n), fill = "grey70") +
			geom_line(aes(time,median/n),color = "black") + geom_line(aes(time, idw/n), color = "blue", alpha = 0.7, size = 0.5)+
			ylim(c(0,1))+ylab("Hypoxia Extent")+theme_bw()+
			ggtitle(paste0("When DO <", units[[threshold]])) +
			theme(axis.title.x = element_text(size=9),axis.title.y = element_text(size=9), plot.title = element_text(size=9))
		
		pdf(sprintf("../output/results/%d_%s_%s_%d_extent_%s.pdf",year, aggType, method, r, threshold), 
				width = 5, height = 3)
		print(p)
		dev.off()
	}
}


getDecompositionResults <- function(year, aggType, r){
	# function to plot basis info
	folderName <- sprintf("../output/%d_%s_%s_%d/",year, aggType, "Baye", r)
	
	model <- readRDS(paste0(folderName,"basisModelRes.rds"))
	timeIndex <- index(model$residuals$samplingData)
	residuals <- model$residuals$samplingData
	
	loggerNames <- names(residuals)
	names(residuals) <- c(1:ncol(residuals))
	
	basis <-attr(model$model,"basis")[,1:r] %>% data.frame()
	names(basis) <- paste0("basis_",c(1:r))
	
	
	basis <- zoo(basis, order.by = timeIndex)
	
	# plot the basis functions
	pdf(sprintf("../output/results/%d_%s_r_%d_basis.pdf",year,aggType, r),	width = 8, height = 3)
	print(plot(basis[,1:3], xlab = "Time",nc = 3 , yax.flip = FALSE,cex.axis = 1.5, cex.lab = 2, main = "",oma = c(5, 0, 2, 0),
		mar = c(0, 5.1, 0, 2.1)))
	dev.off()
		
	# plot coefficients on the map
	myMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",year)) + labs(x = "Longitude", y = "Latitude") 
	
	plot.list1 <- lapply(seq(1,3,1), function(x){
		subData <- model$model[[x]]
		myMap + geom_point(aes(longitude,latitude, color = value),data = subData, size = 4) + 
			scale_color_gradientn(colours = topo.colors(10), name = "Basis\nCoefficients")+ 
			ggtitle(paste0("Basis_", x)) + theme(axis.title.x = element_text(size=11),axis.title.y = element_text(size=11))
	})

	# plot.list2 <- lapply(seq(2,r,2), function(x){
	# 	subData <- model$model[[x]]
	# 	myMap + geom_point(aes(longitude,latitude, color = value),data = subData, size = 4) + 
	# 		scale_color_gradientn(colours = topo.colors(10), name = "Basis\nCoefficients")+ 
	# 		ggtitle(paste0("Basis_", x))
	# })
	args.list <- c(c(plot.list1),list(nrow=1))
	# args.list <- c(c(plot.list1,plot.list2),list(nrow=2))
	pdf(sprintf("../output/results/%d_%s_r_%d_basisCoeff.pdf",year,aggType, r),	
			width = 3*4, height = 3)
	print(do.call(grid.arrange, args.list))
	dev.off()
	
	
	# plot residuals
	pdf(sprintf("../output/results/%d_%s_r_%d_residuals.pdf",year,aggType, r),
			width = 6.5, height = 7)
	print(plot(residuals, xlab = "Time", main = "",nc = 3, oma = c(5, 0, 2, 0), mar = c(0, 4.6, 0, 1.1)))
	dev.off()
 	
 	newRes <- data.frame(residuals)
	names(newRes) <- c(1:ncol(newRes))

	pdf(sprintf("../output/results/%d_%s_r_%d_resBoxplot.pdf",year,aggType, r), width = 5, height = 3)
	par(oma = c(0, 0, 0, 0), mar = c(2, 4, 0.5, 0.5))
	print(boxplot(newRes, ylab = "Residuals (mg/L)"))
	dev.off()
	
	
	# plot the spatio-temporal variogram
	require(ggplot2)
	residuals <- model$residuals$samplingData
	
	#stv <- st_variogram(model$residuals$samplingData, model$residuals$loggerInfo)
	#saveRDS(stv, sprintf("../output/results/%d_%s_r_%d_stVgm.rds",year,aggType, r))
	stv <- readRDS(sprintf("../output/results/%d_%s_r_%d_stVgm.rds",year,aggType, r))
	
	pdf(sprintf("../output/results/%d_%s_r_%d_stVgm.pdf",year,aggType, r),
			width = 5, height = 3)
	par(oma = c(0, 0, 0, 0), mar = c(2, 4, 0.5, 0.5))
	print(plot(stv$vgmModel, map = F))
	dev.off()
}


getSensorWithHypoxiaExtent <- function(year, aggType, method, r){
	# plot the sensor importances in determining hypoxia extent
	erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()
	loggerInfo <- erieDO$loggerInfo
	
	if(method == "idw"){
		if(r>5){
			return()
		}
		folderName <- sprintf("../output/%d_%s_idw/",year, aggType)
		trendModel <- readRDS(paste0(folderName, "trendPredictions.rds"))
		HEList <- readRDS(paste0(folderName,"extent.rds"))
		
	}else{
		folderName <- sprintf("../output/%d_%s_%s_%d/",year, aggType, method, r)
		HEList <- readRDS(paste0(folderName,"extent.rds"))
		trendModel <- readRDS(paste0(folderName, "trendModel.rds"))
	}
	
	grid <- trendModel$grid
	n <- sum(grid$convexIndex)
	threshold <- list(less4 = 4, less2 = 2, less0 = 0.01)
	
	myMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",year))
	
	modelList <- list()
	
	for(thres in c("less0","less2","less4")){
		if(method == "idw"){
			hypoxiaExtentStream <- HEList[[thres]]/n
		}else{
			hypoxiaExtentStream <- HEList[[thres]]$median/n
		}
		
		lmModel <- identifyKeyLoggersToHypoxia(erieDO, hypoxiaExtentStream, threshold[[thres]])
		
		loggerInfo$value <- lmModel$coefficients
		
		# plot the coefficients 
		p <- myMap + geom_point(aes(longitude,latitude, color = value),data = loggerInfo, size = 4) + 
			scale_color_gradientn(colours = topo.colors(10), name = "Coefficients")
		
		png(sprintf("../output/results/%d_%s_%s_%d_%s_sensorLMCoeff.png",year,aggType, method,r, thres),
				width = 600, height = 400, pointsize = 20)
		print(p)
		dev.off()
		
		modelList[[thres]] <- lmModel
		
	}
	
	saveRDS(modelList, sprintf("../output/results/%d_%s_%s_%d_sensorLM_HE_fit.rds",year,aggType, method, r))
	
}



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
	loggerNames <- paste0("logger_", loggerNames)
	
	# using whether loggers are less than DO threshold
	hypoxiaMatrix <- (lakeDOObj$samplingData < DOthreshold) %>% 
				data.frame() %>% 
				data.matrix() %>% 
				data.frame()
	
	names(hypoxiaMatrix) <- loggerNames
	# or using the pure values (which needs a non-linear model to do)
	hypoxiaMatrix$Y <- hypoxiaExtentStream

	formu <- paste0("Y~",paste(loggerNames, collapse = "+"),"-1") %>% formula()
	
	model <- lm(formu, data = hypoxiaMatrix)

	return(model)
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


	return(mergedDF)
	# pdf(sprintf("../output/results/%d_%s_r_%d_basisCoeff.pdf",year,aggType, r),	
	# pdf("~/Desktop/tmp.pdf")
			# width = 4* basis_r, height = 6)
	# print(do.call(grid.arrange, args.list))
	# dev.off()
	
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
					for(r in rList){
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
	fullRes <- readRDS(sprintf("%s/results/fullRes.rds",outputBaseName))
	fullRes$method = as.character(fullRes$method)
	# plot the CV boxplot 
	fullRes[fullRes$method == "idw","method"] <- "IDW"
	fullRes[fullRes$method == "Reml","method"] <- "MLE"

	fullRes$method <- factor(fullRes$method, levels = c("IDW","MLE","Baye"))
	
	for(year_ in c(2014,2015)){
		for(aggType_ in c("daily","hourly")){

			p_rmse <- dplyr::filter(fullRes, aggType == aggType_, year == year_) %>% 
				ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = rmse, fill = method),size = I(0.5), position = position_dodge(width = 0.8),outlier.size = 0.5) + 
				xlab("Model Parameter") + ylab("RMSE") + theme_bw() + 
				scale_fill_manual(name="Method",
                         breaks=c("IDW", "MLE","Baye"),
                         values = c("#999999","#E69F00","#56B4E9")) + 
				theme(axis.title.x = element_text(size=8),axis.title.y = element_text(size=8))
				
			pdf(sprintf("%s/results/%d_%s_CV_rmse.pdf",outputBaseName, year_, aggType_),width = 3.25, height = 2)
			print(p_rmse)
			dev.off()
			
			p_withBound <- dplyr::filter(fullRes, aggType == aggType_, year == year_, method != "IDW") %>% 
				ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = withinBoundRatio, fill = method), size = I(0.5), position = position_dodge(width = 0.8),outlier.size = 0.5) + 
				xlab("Model Parameter") + ylab("CI Coverage") + theme_bw() + 
				scale_fill_manual(name="Method",
                         breaks=c("MLE","Baye"),
                         values = c("#E69F00","#56B4E9") ) + 
				theme(axis.title.x = element_text(size=9),axis.title.y = element_text(size=9))
			
			pdf(sprintf("%s/results/%d_%s_CV_withinBound.pdf", outputBaseName, year_, aggType_),width = 3.25, height = 2)
			print(p_withBound)
			dev.off()	
			
		}
	}
		#dplyr::filter(fullRes, aggType == "daily", year == 2015) %>% 
		#	ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = rmse, fill = method), position = position_dodge(width = 0.8)) + 
		#	xlab("Basis Number (with out bias basis)")
		
		#dplyr::filter(fullRes, aggType == "hourly", year == 2015) %>% 
		#	ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = rmse, fill = method), position = position_dodge(width = 0.8))
		
		
		#dplyr::filter(fullRes, aggType == "daily", year == 2014) %>% 
		#	ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = rmse, fill = method), position = position_dodge(width = 0.8))
		
		#dplyr::filter(fullRes, aggType == "hourly", year == 2014) %>% 
		#	ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = rmse, fill = method), position = position_dodge(width = 0.8)) +
		#	xlab("Basis Number (with out bias basis)")
		
		
		
		#dplyr::filter(fullRes, aggType == "daily", year == 2015, method != "idw") %>% 
		#	ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = withinBoundRatio, fill = method), position = position_dodge(width = 0.8))
		
		#dplyr::filter(fullRes, aggType == "hourly", year == 2015, method != "idw") %>% 
		#	ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = withinBoundRatio, fill = method), position = position_dodge(width = 0.8))
		
		
		#dplyr::filter(fullRes, aggType == "daily", year == 2014, method != "idw") %>% 
		#	ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = withinBoundRatio, fill = method), position = position_dodge(width = 0.8))
		
		#dplyr::filter(fullRes, aggType == "hourly", year == 2014, method != "idw") %>% 
		#	ggplot(data = .)+geom_boxplot(aes(x = factor(r), y = withinBoundRatio, fill = method), position = position_dodge(width = 0.8)) +
		#	xlab("Basis Number (with out bias basis)")
			

	return(fullRes)
}

plotCVOnMap <- function(year_, aggType_, method_, r_){
	# function to plot CV rmse and withinBoundRatio on the map
	library(ggmap)

	if(method_ == "idw") r_ = 5
	erieDO <- getLakeDO(year_, "B", aggType_) %>% na.omit()
	res <- readRDS(sprintf("%s/results/fullRes.rds",outputBaseName)) %>% 
	subset(year == year_ & aggType == aggType_ & method == method_ & r == r_) %>% 
	merge(erieDO$loggerInfo, by.x = "cv_loggerID", by.y = "loggerID", all.x = TRUE)

	# lonRange <- range(res$longitude)
	# latRange <- range(res$latitude)
	# bbox <- make_bbox(lonRange,latRange,f = 0.2)
	# myMap <- get_map(location=bbox, source="google",crop=FALSE) %>% ggmap()
	# saveRDS(myMap,"erieGoogleMap.rds")

	myMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",year_)) + labs(x = "Longitude", y = "Latitude")
	if(method_ == "idw"){
		p <- myMap + geom_point(aes(longitude, latitude, color = rmse), data = res, size = 5)
	}else{
		# p <- myMap + geom_point(aes(longitude, latitude, size = withinBoundRatio, color = rmse), data = res)
		p <- myMap + geom_point(aes(longitude, latitude, size = rmse, color = withinBoundRatio), data = res)
	}
	p <- p + scale_color_gradientn(colours=rainbow(6,end=4/6), name = "CI Coverage") + scale_size_continuous(name = "RMSE")
	return(p)
	#p_rmse <-  myMap + geom_point(aes(longitude, latitude, color = rmse), data = res)
	#p_inBound <- myMap + geom_point(aes(longitude, latitude, color = withinBoundRatio), data = res)
}


RMSE <- function(x,y){
	return(round(sqrt(mean((x-y)^2)),digits = 3))
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

	
	p <- ggplot(data = df) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,median),color = "black",size = 0.5) + 
		geom_line(aes(time,true),color = "red", alpha = 0.7, size = 0.5) + 
		ylim(c(-1,20))+ylab("DO (mg/L)")+ xlab("") + theme_bw() +
		ggtitle(sprintf("Logger: %s, RMSE:%.3f, CI Coverage: %.3f",cv_loggerID, rmse, withinBoundRatio)) + 
		theme(axis.title.x = element_text(size=9),axis.title.y = element_text(size=9), plot.title = element_text(size=9))

	pdf(figName, width = 3.5, height = 3)
	print(p)
	dev.off()
	return(list(rmse = rmse, withinBoundRatio = withinBoundRatio))
}	



