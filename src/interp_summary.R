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

	hypoxiaSummary <- zoo(hypoxiaSummary,order.by = timeIndex)
	return(hypoxiaSummary)
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

	hypoxiaExtent <- summaryHypoxia(predictionSimulations,timeIndex)

	area <- attr(grid,"totalArea")
	n <- sum(grid$convexIndex)

	bbox <- make_bbox(range(grid$longitude),range(grid$latitude),f = 0.1)
	myMap <- get_map(location=bbox, source="google",crop=FALSE)
	p <- ggmap(myMap)
	p <- p+geom_tile(aes(longitude,latitude),data = subset(grid,convexIndex == 1), alpha = I(0.5))
	p <- p+geom_point(aes(longitude,latitude),data = loggerInfo)+ggtitle(paste("Interpolation Area:",round(area,2),"km^2"))
	
	png(paste(outputFolder,"interpolationArea.png",sep = ""))
	print(p)
	dev.off()

	write.csv(x = hypoxiaExtent$less0/n,file = paste0(outputFolder,"hypoxia_less0.csv"))
	write.csv(x = hypoxiaExtent$less2/n,file = paste0(outputFolder,"hypoxia_less2.csv"))
	write.csv(x = hypoxiaExtent$less4/n,file = paste0(outputFolder,"hypoxia_less4.csv"))
}
