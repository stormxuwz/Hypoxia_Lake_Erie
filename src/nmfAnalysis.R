
NMFReconstruct <- function(year, aggType, basis_r, loggerID, new = FALSE,...){
	erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()
	timeIdx <- index(erieDO$samplingData)

	method <- list(...)$method
	if(is.null(method)){
    	method <- "brunet"
  	}

	methodAlias  <- method
	if(method == "snmf/r"){
	  methodAlias <- "snmf-r"
	}else if(method == "snmf/l"){
	  methodAlias <- "snmf-l"
	}
	
	resFileName <- sprintf("%s/results/cluster/nmfBasis_%d_%s_%d_%s.rds",outputBaseName, year,aggType,basis_r, methodAlias)
	if(file.exists(resFileName) & !new){
		nmfDecomp <- readRDS(resFileName) %>% NMF_scale()
	}else{
		print("no previous results exist, new calculation")
		nmfDecomp <- NMF_basis(erieDO$samplingData,basis_r,...) %>% NMF_scale()
		saveRDS(nmfDecomp,resFileName)
	}


	reconstructed <- nmfDecomp$basis %>% data.frame()

	for(i in 1:ncol(nmfDecomp$basis)){
		reconstructed[,i] <- reconstructed[,i] * nmfDecomp$coef[i, loggerID]
	}

	reconstructed$time <- as.POSIXct(timeIdx)
	reconstructed_basis <- melt(reconstructed, id.vars = "time")
	reconstructed$y <- as.numeric(erieDO$samplingData[,loggerID])

	ggplot() + geom_bar(mapping = aes(x = time, y =value, fill = variable),stat = "identity", data = reconstructed_basis) +
	geom_line(aes(time, y), data = reconstructed)


	# dygraph(reconstructed) %>%  dyOptions(stackedGraph = TRUE)
	
	# reconstructed <- zoo(reconstructed, order.by = timeIdx)

	# res <- data.frame(nmf = reconstructed[,as.character(loggerID)], y = as.numeric(erieDO$samplingData[,loggerID]), time = timeIdx)
	# print(qplot(timeIdx, nmf, data = res) + geom_line(aes(timeIdx, y), data = res, color = "red"))
}

NMF_analysis <- function(year, aggType, basis_r = NULL,explore = FALSE, new = FALSE, ...){
	erieDO <- getLakeDO(year, "B", aggType) %>% filterEPA(changeColumnName = FALSE) %>% na.omit()
	timeIdx <- index(erieDO$samplingData)
	baseMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",year)) + labs(x = "Longitude", y = "Latitude")
		 
	d <- erieDO$samplingData
	
	require(NMF)
	
	method <- list(...)$method
 	if(is.null(method)){
    	method <- "brunet"
  	}
	
	methodAlias  <- method
	if(method == "snmf/r"){
		methodAlias <- "snmf-r"
	}else if(method == "snmf/l"){
		methodAlias <- "snmf-l"
	}

	if(explore){
	  print("doing exploration")
		exploreRank <- nmfEstimateRank(as.matrix(d), c(2:6), method = method)
		saveRDS(exploreRank,sprintf("%s/results/cluster/nmfEstimateRank_%d_%s_%s.rds",outputBaseName, year,aggType, methodAlias))
	}
	# exploreRank <- readRDS(sprintf("%s/results/cluster/nmfEstimateRank_%d_%s.rds",outputBaseName, year,aggType))
	
	if(is.null(basis_r)){
	  if(year == 2015){
	    basis_r = 5
	  }else if(year == 2014){
	    basis_r = 5
	  }else if(year == 2016){
	    basis_r = 6
	  }
	}

	resFileName <- sprintf("%s/results/cluster/nmfBasis_%d_%s_%d_%s.rds",outputBaseName, year,aggType,basis_r, methodAlias)
	if(file.exists(resFileName) & !new){
		nmfDecomp <- readRDS(resFileName) %>% NMF_scale()
	}else{
		print("no previous results exist, new calculation")
		nmfDecomp <- NMF_basis(erieDO$samplingData,basis_r,...) %>% NMF_scale()
		saveRDS(nmfDecomp,resFileName)
	}
	
	tmp <- data.frame(nmfDecomp$coef)
	names(tmp) <- colnames(nmfDecomp$coef)
	tmp$basis  <- 1:nrow(tmp)
	tmp <- melt(tmp,id.vars = "basis") 
	mergedDF <- merge(tmp,erieDO$loggerInfo,by.x = "variable", by.y = "loggerID")
	
	
	# plot 
	plot.list1 <- lapply(1:basis_r, function(r){
		subData <- subset(mergedDF,basis == r)
		baseMap + geom_point(aes(longitude,latitude, color = value),data = subData, size = 4) + 
			scale_color_gradient(low = "white", high = "blue",name = "Basis\nCoefficients",limit = range(0, max(mergedDF$value)+0.1))+ 
			ggtitle(paste0("Basis_", r)) + theme(axis.title.x = element_text(size=11),axis.title.y = element_text(size=11))
	})
	plot.list2 <- lapply(1:basis_r, function(r){
		newDf <- data.frame(DO = nmfDecomp$basis[,r],time = timeIdx)
		ggplot(data = newDf) + geom_line(aes(time, DO)) + theme_bw()
	})
	args.list <- c(c(plot.list1,plot.list2),list(nrow=2))
	
	require(gridExtra)
	pdf(sprintf("%s/results/cluster/cluster_%d_%s_%d_%s.pdf",outputBaseName,year,aggType, basis_r, methodAlias),
			width = 4*basis_r, height = 6)
	print(do.call(grid.arrange, args.list))
	dev.off()
}
