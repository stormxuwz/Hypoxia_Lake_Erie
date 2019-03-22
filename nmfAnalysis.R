rm(list= ls())
setwd("/Users/wenzhaoxu/Developer/Hypoxia/Hypoxia_Lake_Erie")
source("src/database.R")
source("src/classDef.R")
source("src/helper.R")
source("src/ecoAnalysis.R")
source("src/basisDecomposition.R")


dbConfig <- list(dbname = "DO", username="root", password="XuWenzhaO", host="127.0.0.1")
varUnit <- list(DO="DO(mg/L)",Temp="Temperature(C)")

nearShoreSite <- list(
nearShoreSite_2014 = c("ASH003","CLE002","CLE001","GEN001","ESL001","CND001","CND002","CND003"),
nearShoreSite_2015 = c("CLE003","CLE001","ASH003","PRG001","PST006","WTY004","WTY002","LRG001"),
nearShoreSite_2016 = c("CBG_83","CBG_43","CBG_55","CBG_94","CBG_11","CBG_52","CBG_73")
)

NMFReconstruct <- function(year, aggType, basis_r, loggerID, new = FALSE,...){
	# erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()

	erieDO <- getLakeDO(year, "B", aggType) %>% 
	filterSites(siteList = nearShoreSite[[paste("nearShoreSite", year, sep = "_")]],changeColumnName = FALSE) %>% 
	na.omit()

	erieDO$samplingData <- erieDO$samplingData[-c(1:3),]
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
		# print("no previous results exist, new calculation")
		# nmfDecomp <- NMF_basis(erieDO$samplingData,basis_r,...) %>% NMF_scale()
		# saveRDS(nmfDecomp,resFileName)
	}

	reconstructed <- nmfDecomp$basis %>% data.frame()

	for(i in 1:ncol(nmfDecomp$basis)){
		reconstructed[,i] <- reconstructed[,i] * nmfDecomp$coef[i, loggerID]
	}

	reconstructed$time <- as.POSIXct(timeIdx)
	reconstructed_basis <- melt(reconstructed, id.vars = "time")
	reconstructed$y <- as.numeric(erieDO$samplingData[,loggerID])

	print(ggplot() + geom_col(mapping = aes(x = time, y =value, fill = variable),stat = "identity", data = reconstructed_basis) +
	geom_line(aes(time, y), data = reconstructed)+ggtitle(loggerID))


	# dygraph(reconstructed) %>%  dyOptions(stackedGraph = TRUE)
	# reconstructed <- zoo(reconstructed, order.by = timeIdx)
	# res <- data.frame(nmf = reconstructed[,as.character(loggerID)], y = as.numeric(erieDO$samplingData[,loggerID]), time = timeIdx)
	# print(qplot(timeIdx, nmf, data = res) + geom_line(aes(timeIdx, y), data = res, color = "red"))
}

NMF_analysis <- function(year, aggType, erieDO, basis_r = NULL,explore = FALSE, new = FALSE, ...){
	
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
	  	tmpd <- as.matrix(d)
	  	tmpd <- ifelse(tmpd<0.01, 0.01, tmpd)
		exploreRank <- nmfEstimateRank(tmpd, c(2:6), method = method)
		saveRDS(exploreRank,sprintf("%s/results/cluster/nmfEstimateRank_%d_%s_%s.rds",outputBaseName, year,aggType, methodAlias))
		print(plot(exploreRank))
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



methodList <- c("brunet","snmf/l","snmf/r")

# EPA site analysis
outputBaseName <- "/Users/wenzhaoxu/Developer/Hypoxia/output_nonEPA"
createFolder(outputBaseName)
createFolder(paste0(outputBaseName,"/results/cluster"))
agg <- "hourly"	

for(year in c(2014:2016)){
	d <- getLakeDO(year, "B", agg) %>% 
	filterSites(siteList = nearShoreSite[[paste("nearShoreSite", year, sep = "_")]],changeColumnName = FALSE) %>% 
	na.omit()
	d$samplingData <- d$samplingData[-c(1:3),]
	for(r in c(2:5)){
		for(method in methodList){
			NMF_analysis(year, agg,d, basis_r = r, method = method, new = TRUE,explore = FALSE)
		}
	}
}

# EPA site analysis
outputBaseName <- "/Users/wenzhaoxu/Developer/Hypoxia/output_EPA"
createFolder(outputBaseName)
createFolder(paste0(outputBaseName,"/results/cluster"))

for(year in c(2014:2016)){
	d <- getLakeDO(year, "B", agg) %>% 
	filterSites(siteList = EPASite,changeColumnName = FALSE) %>% 
	na.omit()
	d$samplingData <- d$samplingData[-c(1:3),]
	for(r in c(2:5)){
		for(method in methodList){
			NMF_analysis(year, agg,d, basis_r = r, method = method, new = TRUE,explore = FALSE)
		}
	}
}


outputBaseName <- "/Users/wenzhaoxu/Developer/Hypoxia/output_all"
createFolder(outputBaseName)
createFolder(paste0(outputBaseName,"/results/cluster"))

for(year in c(2014:2016)){
	d <- getLakeDO(year, "B", agg) %>% na.omit()
	d$samplingData <- d$samplingData[-c(1:2),]
	for(r in c(2:10)){
		for(method in methodList){
			NMF_analysis(year, agg,d, basis_r = r, method = method, new = TRUE,explore = FALSE)
		}
	}
}

# data2014 <- getLakeDO(2014, "B", "daily") %>% filterSites(siteList = nearShoreSite_2014,changeColumnName = FALSE) %>% na.omit()
# data2014$samplingData <- data2014$samplingData[-c(1:3),]
# NMF_analysis(2014, "daily",data2014, basis_r = 2, new = TRUE,explore = TRUE)
# NMF_analysis(2014, "daily",data2014, basis_r = 2, method = "snmf/r", new = TRUE, explore = TRUE)
# NMF_analysis(2014, "daily",data2014, basis_r = 2, method = "snmf/l", new = TRUE, explore = TRUE)

# data2015 <- getLakeDO(2015, "B", "daily") %>% filterSites(siteList = nearShoreSite_2015, changeColumnName = FALSE) %>% na.omit()
# data2015$samplingData <- data2015$samplingData[-c(1:3), ]
# NMF_analysis(2015, "daily",data2015, basis_r = 2, new = TRUE, explore = TRUE)
# NMF_analysis(2015, "daily",data2015, basis_r = 2, method = "snmf/r", new = TRUE, explore = TRUE)
# NMF_analysis(2015, "daily",data2015, basis_r = 2, method = "snmf/l", new = TRUE, explore = TRUE)

# # data2016 <- getLakeDO(2016, "B", "daily") %>% filterSites(siteList = nearShoreSite_2016, changeColumnName = FALSE) %>% na.omit()
# # data2016$samplingData <- data2016$samplingData[-c(1:3),]
# NMF_analysis(2016, "daily",data2016, basis_r = 3, new=TRUE, explore = TRUE)
# NMF_analysis(2016, "daily",data2016, basis_r = 3, method = "snmf/r", new = TRUE, explore = TRUE)
# NMF_analysis(2016, "daily",data2016, basis_r = 3, method = "snmf/l", new = TRUE, explore = TRUE)




