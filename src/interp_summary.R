summaryHypoxia <- function(predictionSimulations,timeIndex){
	# predictionSimulations is a list, each element is a [time, grid] matrix, 
	# representing the interpolation data

	simNum <- length(predictionSimulations)
	TimeN <- nrow(predictionSimulations[[1]])
	locationNum <- ncol(predictionSimulations[[1]])
	
	# initialize the hypoxia extent corresponding to below 0.01, 2 and 4 mg/L
	hypoxiaExtent_0 <- matrix(0,TimeN,simNum)  
	hypoxiaExtent_2 <- hypoxiaExtent_0	
	hypoxiaExtent_4 <- hypoxiaExtent_0
	
	hypoxiaExtent <- list()

	for(i in 1:simNum){
		hypoxiaExtent_0[,i] <- rowSums(predictionSimulations[[i]]<0.01,na.rm = TRUE)
		hypoxiaExtent_2[,i] <- rowSums(predictionSimulations[[i]]<2,na.rm = TRUE)
		hypoxiaExtent_4[,i] <- rowSums(predictionSimulations[[i]]<4,na.rm = TRUE)
	}
	
	hypoxiaExtent[["less0"]] <- summarySimulation(hypoxiaExtent_0, timeIndex)
	hypoxiaExtent[["less2"]] <- summarySimulation(hypoxiaExtent_2, timeIndex) 
	hypoxiaExtent[["less4"]] <- summarySimulation(hypoxiaExtent_4, timeIndex) 
	
	return(hypoxiaExtent)
}


summaryHypoxia_2 <- function(predictionSimulations, timeIndex){
	simNum <- length(predictionSimulations)
	TimeN <- nrow(predictionSimulations[[1]])

	hypoxiaExtent_0 <- matrix(0,TimeN,simNum)  
	hypoxiaExtent_2 <- hypoxiaExtent_0	
	hypoxiaExtent_4 <- hypoxiaExtent_0
	hypoxiaExtent <- list()
	for(i in 1:simNum){
		hypoxiaExtent_0[,i] <- predictionSimulations[[i]]$less0
		hypoxiaExtent_2[,i] <- predictionSimulations[[i]]$less2
		hypoxiaExtent_4[,i] <- predictionSimulations[[i]]$less4
	}

	hypoxiaExtent[["less0"]] <- summarySimulation(hypoxiaExtent_0, timeIndex)
	hypoxiaExtent[["less2"]] <- summarySimulation(hypoxiaExtent_2, timeIndex) 
	hypoxiaExtent[["less4"]] <- summarySimulation(hypoxiaExtent_4, timeIndex) 

	return(hypoxiaExtent)
} 


summarySimulation <- function(hypoxiaSimulations,timeIndex){
	# hypoxiaSimulations is a table of shape [Time, nSimulations]

	if(ncol(hypoxiaSimulations)>1){
		hypoxiaSummary <- data.frame(
			median = apply(hypoxiaSimulations,1, median),
			upper = apply(hypoxiaSimulations,1, quantile, probs = c(0.95), na.rm = TRUE),
			lower = apply(hypoxiaSimulations,1, quantile, probs = c(0.05), na.rm = TRUE)
		)
	}else{
		hypoxiaSummary <- data.frame(
			median = hypoxiaSimulations[,1],
			upper = hypoxiaSimulations[,1],
			lower = hypoxiaSimulations[,1]
		)
	}

	# hypoxiaSummary <- zoo(hypoxiaSummary,order.by = timeIndex)
	hypoxiaSummary$time <- timeIndex
	return(hypoxiaSummary)
}


getOneSimulation <- function(year,timeAggType,method,simInd,timeInd,...){
	require(reshape2)
	require(dplyr)

	loggerInfo <- retriveGeoData(year,"B") %>% arrange(loggerID)
	# the data used to interpolation
	#timeRange = c("2014-07-01","2014-09-01")
	timeRange = NULL # no specific time range
	data <- retriveLoggerData(loggerInfo$loggerID,year,"DO",timeAggType,"AVG",timeRange = timeRange,transform = TRUE) %>% na.omit()  # remove the data
	time <- index(data)
	loggerNames  <- as.numeric(colnames(data))
	subData <- data.frame(logger = loggerNames,DO = as.numeric(data[timeInd,]))
	subData <- merge(subData, loggerInfo,by.x = "logger",by.y = "loggerID") %>% rename(value = DO)

	if(method == "basis"){
		metaFolder <<- sprintf("../meta_%d_%s_%s_%s_%d/", 
			year,timeAggType,method,list(...)$fitMethod,list(...)$r)
		outputFolder <<- sprintf("../output_%d_%s_%s_%s_%d/", 
			year,timeAggType,method,list(...)$fitMethod, list(...)$r)
	}else{
		metaFolder <<- sprintf("../meta_%d_%s_%s/", year,timeAggType,method)
		outputFolder <<- sprintf("../output_%d_%s/", year,timeAggType,method)
	}

	tmp <- readRDS(paste0(metaFolder,"trend.rds"))
	residual_interp <- readRDS(paste0(metaFolder,"residual_prediction.rds"))

	trendBasisCoeff <- tmp$basisIntRes
	grid <- tmp$grid

	coeffPredList <- trendBasisCoeff$trendCoeff  # list of [#grid, nSim]
	basis <- trendBasisCoeff$basis

	inds <- readRDS(paste0(metaFolder,"inds.rds"))

	print("read results complete, calculting DO")

	prediction = 0 
	for(i in 1:ncol(basis)){
		# get the time series of each grid
		prediction <- prediction + 
			basis[,i] %*% t(coeffPredList[[i]][,inds[simInd,i]]) 
	}

	# add residual interpolation part
	# for deterministic interpolation on the residuals, choose [1]
	prediction <- prediction + residual_interp[[1]]  
	prediction <- prediction*(prediction>0)
	area <- attr(grid,"totalArea")
	grid$value <- prediction[timeInd,]

	# bbox <- make_bbox(range(grid$longitude),range(grid$latitude),f = 0.1)
	# myMap <- get_map(location=bbox, source="google",crop=FALSE)
	# p <- ggmap(myMap)
	p <- readRDS("map.rds")
	p <- p+geom_tile(aes(longitude,latitude,fill = value),data = subset(grid,convexIndex == 1 & bathymetry<(-12)))
	p <- p+scale_fill_gradientn(colours = terrain.colors(10),limit = c(0,15))
	p <- p+geom_point(aes(longitude,latitude,fill = value),data = subData,shape = I(21), colour = I("black"),size = I(3),stroke = I(1))
	p <- p+ggtitle(paste("Interpolation Area:",round(area,2),"km^2"))
	print(p)
	print(sum(subset(grid,convexIndex == 1 & bathymetry<(-12))$value<2))
	print(sum(subData$value<2))
	return(list(grid = grid, logger = subData, time = time[timeInd]))
}	




summary_plot <- function(year,timeAggType,method,...){
	loggerInfo <- retriveGeoData(year,"B") %>% arrange(loggerID)

	if(method == "basis"){
		metaFolder <<- sprintf("../meta_%d_%s_%s_%s_%d/", 
			year,timeAggType,method,list(...)$fitMethod,list(...)$r)
		outputFolder <<- sprintf("../output_%d_%s_%s_%s_%d/", 
			year,timeAggType,method,list(...)$fitMethod, list(...)$r)
	}else{
		metaFolder <<- sprintf("../meta_%d_%s_%s/", year,timeAggType,method)
		outputFolder <<- sprintf("../output_%d_%s/", year,timeAggType,method)
	}

	predictionRes <- readRDS(paste0(outputFolder,"hypoxiaExtent.rds"))
	
	predictionSimulations <- predictionRes$res
	timeIndex <- predictionRes$timeIndex
	grid <- predictionRes$grid

	hypoxiaExtent <- summaryHypoxia_2(predictionSimulations,timeIndex)

	area <- attr(grid,"totalArea")
	n <- sum(grid$convexIndex)

	bbox <- make_bbox(range(grid$longitude),range(grid$latitude),f = 0.1)
	myMap <- get_map(location=bbox, source="google",crop=FALSE)
	p <- ggmap(myMap)
	p <- p+geom_tile(aes(longitude,latitude),data = subset(grid,convexIndex == 1), alpha = I(0.5))
	p <- p+geom_point(aes(longitude,latitude),data = loggerInfo)+
	ggtitle(paste("Interpolation Area:",round(area,2),"km^2"))
	
	png(paste(outputFolder,"interpolationArea.png",sep = ""))
	print(p)
	dev.off()

	hypoxiaExtent$less0[,1:3] <- hypoxiaExtent$less0[,1:3]/n
	hypoxiaExtent$less2[,1:3] <- hypoxiaExtent$less2[,1:3]/n
	hypoxiaExtent$less4[,1:3] <- hypoxiaExtent$less4[,1:3]/n

	write.csv(x = hypoxiaExtent$less0,file = paste0(outputFolder,"hypoxia_less0.csv"))
	write.csv(x = hypoxiaExtent$less2,file = paste0(outputFolder,"hypoxia_less2.csv"))
	write.csv(x = hypoxiaExtent$less4,file = paste0(outputFolder,"hypoxia_less4.csv"))


	p_less0 = ggplot(data = hypoxiaExtent$less0) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,median),color = "black") + 
		ylim(c(0,1))+ylab("DO")+
		ggtitle(paste0("Hypoxia extent DO<0.01 mg/L"))

	p_less2 = ggplot(data = hypoxiaExtent$less2) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,median),color = "black") + 
		ylim(c(0,1))+ylab("DO")+
		ggtitle(paste0("Hypoxia extent DO<2 mg/L"))

	p_less4 = ggplot(data = hypoxiaExtent$less4) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,median),color = "black") + 
		ylim(c(0,1))+ylab("DO")+
		ggtitle(paste0("Hypoxia extent DO<4 mg/L"))

	pdf(sprintf("%splot_%s.pdf",outputFolder,"HE_less0.pdf"),width = 6, height = 4)
	print(p_less0)
	dev.off()

	pdf(sprintf("%splot_%s.pdf",outputFolder,"HE_less2.pdf"),width = 6, height = 4)
	print(p_less2)
	dev.off()

	pdf(sprintf("%splot_%s.pdf",outputFolder,"HE_less4.pdf"),width = 6, height = 4)
	print(p_less4)
	dev.off()
}



