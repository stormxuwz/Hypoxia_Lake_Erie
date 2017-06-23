rm(list= ls())
setwd("/Users/wenzhaoxu/Developer/Hypoxia/Hypoxia_Lake_Erie")
source("src/database.R")
source("src/classPrediction.R")
source("src/classDef.R")
source("src/classReConstruct.R")
source("src/classSummary.R")
source("src/helper.R")
source("src/postAnalysis.R")
source("src/basisDecomposition.R")

dbConfig <- list(dbname = "DO", username="root", password="XuWenzhaO", host="127.0.0.1")
varUnit <- list(DO="DO(mg/L)",Temp="Temperature(C)")
outputBaseName <- "/Users/wenzhaoxu/Developer/Hypoxia/output/"
mapDx <- 0.025
mapDy <- 0.025

rList <- c(5,10,15)

main <- function(year = 2014, aggType = "daily"){
	#year <- 2014
	#aggType <- "daily"
	#r <- 5
	# trend <- ~coords[,"x"]+coords[,"y"] + bathymetry + I(bathymetry^2)
	trend <- ~coords[,"x"]+ coords[,"y"] + bathymetry + I(bathymetry^2)

	erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()
	grid <- createGrid(erieDO$loggerInfo, mapDx, mapDy) # grid size: 0.025, 0.025 long,lat
	
	# using IDW
	print("interpolationg using IDW")
	hypoxia_idw <- predict(obj = erieDO, 
				grid = grid, 
				method = "idw", 
				predictType = "extent",
				metaFolder = sprintf("%s%d_%s_idw/",outputBaseName, year, aggType), 
				nmax = 5)

	saveRDS(hypoxia_idw, sprintf("%s%d_%s_idw/extent.rds",outputBaseName, year, aggType))

	for(r in rList){
		print(paste("r:",r))
		print("interpolating using Reml")
		hypoxia_reml <- predict(obj = erieDO,
					grid = grid, 
					method = "Reml", 
					predictType = "extent", 
					trend = trend, 
					r = r, 
					totalSim = 1000,
					nmax = 5,
					metaFolder = sprintf("%s%d_%s_Reml_%d/",outputBaseName, year, aggType, r))
		saveRDS(hypoxia_reml, sprintf("%s%d_%s_Reml_%d/extent.rds",outputBaseName, year, aggType, r))
	
		print("interpolating using Baye")
		hypoxia_baye <- predict(obj = erieDO, 
					grid = grid, 
					method ="Baye", 
					predictType = "extent", 
					trend = trend, 
					r = r, 
					totalSim = 1000,
					nmax = 5,
					metaFolder = sprintf("%s%d_%s_Baye_%d/",outputBaseName, year, aggType, r))

		saveRDS(hypoxia_baye, sprintf("%s%d_%s_Baye_%d/extent.rds",outputBaseName, year, aggType, r))
	}
}

main_CV <- function(year = 2014, aggType = "daily"){
	trend <- ~coords[,"x"]+coords[,"y"]+bathymetry+I(bathymetry^2)
	erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()

	for(i in 1:nrow(erieDO$loggerInfo)) {
	#for(i in 1:2) {
		cv_loggerID <- erieDO$loggerInfo$loggerID[i]
		subErieDO <- subset(erieDO, cv_loggerID)
		grid <- erieDO$loggerInfo[i,]
		grid$convexIndex <- 1

		trueDO <- erieDO$samplingData[,cv_loggerID]
		timeIdx <- index(trueDO)

		metaFolder <- sprintf("%s%d_%s_idw/cv/%s/",outputBaseName, year, aggType, cv_loggerID)
		
		hypoxia_idw <- predict(obj = subErieDO, 
					grid = grid, 
					method = "idw", 
					predictType = "simulations", 
					metaFolder = metaFolder,
					nmax = 5)
		
		statsSummary <- data.frame(median = hypoxia_idw, lower = hypoxia_idw, 
			upper = hypoxia_idw, true = as.numeric(trueDO)) %>% 
			zoo(order.by = timeIdx)

		saveRDS(statsSummary, sprintf("%s/simulations_stat.rds",metaFolder))

		for(r in rList){
			metaFolder <- sprintf("%s%d_%s_Reml_%d/cv/%s/",outputBaseName, year, aggType, r, cv_loggerID)
			
			hypoxia_reml_cv <- predict(obj =subErieDO, 
						grid = grid, 
						method = "Reml", 
						predictType = "simulations", 
						trend = trend, 
						r = r, 
						totalSim = 1000,
						metaFolder = metaFolder,
						nmax = 5)

			# hypoxia_reml_cv <- readRDS(sprintf("%s/simulations.rds",metaFolder))
			statsSummary <- cvUncertainty(hypoxia_reml_cv,trueDO) %>% zoo(order.by = timeIdx)
			saveRDS(statsSummary, sprintf("%s/simulations_stat.rds",metaFolder))
			saveRDS(hypoxia_reml_cv, sprintf("%s/simulations.rds",metaFolder))

			metaFolder <- sprintf("%s%d_%s_Baye_%d/cv/%s/",outputBaseName, year, aggType, r, cv_loggerID)
			hypoxia_baye_cv <- predict(obj =subErieDO, 
						grid = grid, 
						method = "Baye", 
						predictType = "simulations", 
						trend = trend, 
						r = r, 
						totalSim = 1000,
						nmax = 5,
						metaFolder = metaFolder)

			# hypoxia_baye_cv <- readRDS(sprintf("%s/simulations.rds",metaFolder))
			statsSummary <- cvUncertainty(hypoxia_baye_cv,trueDO) %>% zoo(order.by = timeIdx)
			saveRDS(statsSummary, sprintf("%s/simulations_stat.rds",metaFolder))
			saveRDS(hypoxia_baye_cv, sprintf("%s/simulations.rds",metaFolder))
		}
	}
}

main_removingLogger <- function(year = 2014, aggTye = "daily", ...){
	trend <- ~coords[,"x"]+coords[,"y"]+bathymetry+I(bathymetry^2)
	erieDO <- getLakeDO(year, "B", aggType)
	
	grid <- createGrid(erieDO$loggerInfo, 0.025, 0.025)

	removedLoggerID <- c()
	r <- 5

	erieDO <- subset(erieDO, removedLoggerID)

	hypoxia_idw <- predict(erieDO, grid, "idw")
	hypoxia_reml_cv <- predict(erieDO, grid, "Reml", predictType = "expected", trend = trend, r = r, totalSim = 10)
	hypoxia_baye_cv <- predict(erieDO, grid, "Baye", predictType = "extent", trend = trend, r = r, totalSim = 1000)
}

main_analysis <- function(year,aggType){
	# cluster analysis
	# year = 2014
	# aggType = "hourly"

	hypoxiaTimeSummary <- function(predictions, filePrefix){	
		tmp <- hypoxiaAnalysis(predictions,2)
		grid$totolHypoxiaTime <- tmp$hypoxiaTime
		grid$longestHypoxiaTime <- tmp$longestTime
		
		p <- baseMap + 
			geom_tile(aes(longitude,latitude,fill = totolHypoxiaTime),data = subset(grid, convexIndex == 1))+
			scale_fill_gradient(name = "Total\nhypoxic\ntime\n(hour)",low = "cyan",high = "red")

		pdf(sprintf("%s/results/%s_totalHypoxiaCount.pdf",outputBaseName, filePrefix), width = 4, height =3)
		par(oma = c(0, 0, 0, 0))
		print(p)
		dev.off()

		pdf(sprintf("%s/results/%s_longestHypoxiaCount.pdf",outputBaseName, filePrefix),width = 4, height =3)
		par(oma = c(0, 0, 0, 0))
		p <- baseMap + 
			geom_tile(aes(longitude,latitude,fill = longestHypoxiaTime),data = subset(grid, convexIndex == 1))+
			scale_fill_gradient(name = "Longest\nhypoxic\ntime\n(hour)",low = "cyan",high = "red")
		print(p)
		dev.off()
	}


	erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()
	timeIdx <- index(erieDO$samplingData)
	baseMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",year)) + labs(x = "Longitude", y = "Latitude")
	
	grid <- createGrid(erieDO$loggerInfo, mapDx, mapDy)
	
	d <- erieDO$samplingData %>% na.omit() 

	####################
	# perform NMF clustering analysis
	#####################
	if(year == 2015){
		basis_r = 5
	}else{
		basis_r = 3
	}
	require(NMF)
	exploreRank <- nmfEstimateRank(as.matrix(d), c(2:10))
	saveRDS(exploreRank,sprintf("%s/results/cluster/nmfEstimateRank_%d_%s.rds",outputBaseName, year,aggType))
	nmfDecomp <- NMF_basis(d,basis_r)
	
	tmp <- data.frame(nmfDecomp$coef)
	names(tmp) <- colnames(nmfDecomp$coef)
	tmp$basis  <- 1:nrow(tmp)
	tmp <- melt(tmp,id.vars = "basis") 
	
	mergedDF <- merge(tmp,erieDO$loggerInfo,by.x = "variable", by.y = "loggerID")

	plot.list1 <- lapply(1:basis_r, function(r){
		subData <- subset(mergedDF,basis == r)
		baseMap + geom_point(aes(longitude,latitude, color = value),data = subData, size = 4) + 
			scale_color_gradient(low = "white", high = "blue",name = "Basis\nCoefficients",limit = c(0,1))+ 
			ggtitle(paste0("Basis_", r)) + theme(axis.title.x = element_text(size=11),axis.title.y = element_text(size=11))
	})

	plot.list2 <- lapply(1:basis_r, function(r){
		newDf <- data.frame(DO = nmfDecomp$basis[,r],time = timeIdx)
		ggplot(data = newDf) + geom_line(aes(time, DO)) + theme_bw()
	})

	args.list <- c(c(plot.list1,plot.list2),list(nrow=2))
	
	require(gridExtra)
	pdf(sprintf("%s/results/cluster/cluster_%d_%s_%d2.pdf",outputBaseName,year,aggType, basis_r),
		width = 4*basis_r, height = 6)
	print(do.call(grid.arrange, args.list))
	dev.off()
	

	#######################
	# plot hypoxia curve
	#######################
	for(method in c("Reml","Baye")){
		for(r in rList){
			getHypoxiaExtent(year, aggType, method, r)
		}
	}
	
	###################
	# analyze the decomposition results
	#################
	for(r in rList){
		getDecompositionResults(year, aggType, r)
	}
	

	# # get sensors linear regression
	# for(method in c("Reml","Baye","idw")){
	# 	for(r in rList){
	# 		getSensorWithHypoxiaExtent(year, aggType, method, r)
	# 	}
	# }
	
	##########################
	# analyze the hypoxia time
	##########################
	for(method in c("Reml","Baye")){
		for(r in rList){
			fileFolder <- sprintf("%s/%d_%s_%s_%d/", outputBaseName, year, aggType, method, r)

			residualPrediction <- readRDS(paste0(fileFolder, "residualPredictions.rds"))
			predictions <- readRDS(paste0(fileFolder,"trendModel.rds")) %>% 
				reConstruct(residualPrediction = residualPrediction, simulationNum = 0)

			hypoxiaTimeSummary(predictions$predValue, filePrefix = sprintf("%d_%s_%s_%d",year, aggType, method, r) )

		}
	}

	predictions <- readRDS(sprintf("%s/%d_%s_idw/trendPredictions.rds", outputBaseName, year, aggType)) %>%
				reConstruct()
	hypoxiaTimeSummary(predictions, filePrefix = sprintf("%d_%s_%s",year, aggType, "idw") )

	#######################
	# Summarize CV results
	#####################
	for(method in c("Reml","Baye")){
		for(r in rList){
			filePrefix <- sprintf("%d_%s_%s_%d",year, aggType, method, r)

			p <- plotCVOnMap(year, aggType, method, r)
			pdf(sprintf("%s/results/%s_CV_summary.pdf",outputBaseName, filePrefix), width = 6, height = 4)
			print(p)
			dev.off()
		}
	}
	
	filePrefix <- sprintf("%d_%s_%s",year, aggType, "idw")
	p <- plotCVOnMap(year, aggType, "idw", r)
	pdf(sprintf("%s/results/%s_CV_summary.pdf",outputBaseName, filePrefix), width = 6, height = 4)
	print(p)
	dev.off()
}


# function to calculate the CV resutls


# for(year in c(2015)){
# 	for(aggType in c("daily","hourly")){
# 		print(system.time(main(year = year, aggType = aggType)))
# 		print(system.time(main_CV(year = year, aggType = aggType)))
# 	}
# }

resultSummary()


for(year in c(2015)){
	for(aggType in c("hourly","daily")){
		print(system.time(main_analysis(year = year, aggType = aggType)))
	}
}





# for(year in c(2014, 2015)){
# 	for(aggType in c("hourly","daily")){
		
# 		# get hypoxia extent
# 		# for(method in c("Reml","Baye")){
# 		# 	for(r in rList){
# 		# 		getHypoxiaExtent(year, aggType, method, r)
# 		# 	}
# 		# }
		
# 		# analyze the decomposition results
# 		for(r in rList){
# 			getDecompositionResults(year, aggType, r)
# 		}
		
# 		# get sensors linear regression
# 		# for(method in c("Reml","Baye","idw")){
# 		# 	for(r in rList){
# 		# 		getSensorWithHypoxiaExtent(year, aggType, method, r)
# 		# 	}
# 		# }
# 	}
# }







# system.time(main_CV(year = 2014, aggType = "daily"))
# system.time(main_CV(year = 2014, aggType = "hourly"))
# system.time(main_CV(year = 2015, aggType = "daily"))
# system.time(main_CV(year = 2015, aggType = "hourly"))
# resultSummary()


# function to calculate the hypoxia extent
#print(system.time(main(year = 2014, aggType = "daily")))
#system.time(main(year = 2014, aggType = "hourly"))
#system.time(main(year = 2015, aggType = "daily"))
#system.time(main(year = 2015, aggType = "hourly"))


# function to do post analysis
# main_analysis()

# API testing
# basisModel_REML <- basisModel(erieDO, trend, "Reml", 10,  "tmp")
# trendModel <- basisModel_REML$model
# residualModel <- basisModel_REML$residuals %>% idwModel()

# residualPredictions <- predict(residualModel, grid, TRUE) %>% reConstruct()
# trendPrediction <- predict(trendModel, grid) %>% reConstruct(residualPredictions,-1)

# tmp <- predict(trendModel, grid)
# trendPrediction <- reConstruct(tmp, residualPredictions, 0)
# hypoxiaExtent <- summary(tmp, residualPredictions, parallel = TRUE, totalSim = 10, randomSeed = 0)




