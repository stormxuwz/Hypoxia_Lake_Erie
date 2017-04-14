rm(list= ls())
setwd("/Users/wenzhaoxu/Developer/Hypoxia/Hypoxia_Lake_Erie")
source("src/database.R")
source("src/classPrediction.R")
source("src/classDef.R")
source("src/classReConstruct.R")
source("src/classSummary.R")
source("src/helper.R")

dbConfig <- list(dbname = "DO", username="root", password="XuWenzhaO", host="127.0.0.1")
varUnit <- list(DO="DO(mg/L)",Temp="Temperature(C)")

outputBaseName <- "/Users/wenzhaoxu/Developer/Hypoxia/output/"

main <- function(year = 2014, aggType = "daily", ...){
	#year <- 2014
	#aggType <- "daily"
	#r <- 5

	trend <- ~coords[,"x"]+coords[,"y"]+bathymetry+I(bathymetry^2)
	erieDO <- getLakeDO(year, "B", aggType)
	grid <- createGrid(erieDO$loggerInfo, 0.025, 0.025) 

	hypoxia_idw <- predict(obj = erieDO, 
				grid = grid, 
				method = "idw", 
				predictType = "extent",
				metaFolder = sprintf("%s%d_%s_idw/",outputBaseName, year, aggType), 
				nmax = 5)

	
	saveRDS(hypoxia_idw, sprintf("%s%d_%s_idw/extent.rds",outputBaseName, year, aggType))

	for(r in c(5,10,15)){
		hypoxia_reml <- predict(obj = erieDO,
					grid = grid, 
					method = "Reml", 
					predictType = "extent", 
					trend = trend, 
					r = r, 
					totalSim = 10,
					nmax = 5,
					metaFolder = sprintf("%s%d_%s_Reml_%d/",outputBaseName, year, aggType, r))
		saveRDS(hypoxia_reml, sprintf("%s%d_%s_Reml_%d/extent.rds",outputBaseName, year, aggType, r))

		hypoxia_baye <- predict(obj = erieDO, 
					grid = grid, 
					method ="Baye", 
					predictType = "extent", 
					trend = trend, 
					r = r, 
					totalSim = 10,
					nmax = 5,
					metaFolder = sprintf("%s%d_%s_Baye_%d/",outputBaseName, year, aggType, r))

		saveRDS(hypoxia_baye, sprintf("%s%d_%s_Baye_%d/extent.rds",outputBaseName, year, aggType, r))
	}
}

main_CV <- function(year = 2014, aggType = "daily", ...){
	trend <- ~coords[,"x"]+coords[,"y"]+bathymetry+I(bathymetry^2)
	erieDO <- getLakeDO(year, "B", aggType)
	erieDO$samplingData <- na.omit(erieDO$samplingData)

	for(i in 1:nrow(erieDO$loggerInfo)) {
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

		for(r in c(5,10,15)){
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

resultSummary <- function(){
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

	# analyze 2014
	#res_2014 <- subset(fullRes, year == 2014)

	# analyze 2015
#	res_2015 <- subset(fullRes, year == 2015)

	return(fullRes)
}

tmp <- function(year_, aggType_, method_, r_){
	library(ggmap)

	if(method_ == "idw") r_ = 5
	erieDO <- getLakeDO(year_, "B", aggType_)
	res <- readRDS(sprintf("%s/results/fullRes.rds",outputBaseName)) %>% 
	subset(year == year_ & aggType == aggType_ & method == method_ & r == r_) %>% 
	merge(erieDO$loggerInfo, by.x = "cv_loggerID", by.y = "loggerID", all.x = TRUE)

	lonRange <- range(res$longitude)
	latRange <- range(res$latitude)
	bbox <- make_bbox(lonRange,latRange,f = 0.2)
	myMap <- get_map(location=bbox, source="google",crop=FALSE) %>% ggmap()
	saveRDS(myMap,"erieGoogleMap.rds")

	myMap <- readRDS("erieGoogleMap.rds")
	if(method_ == "idw"){
		p <- myMap + geom_point(aes(longitude, latitude, color = rmse), data = res, size = 5)
	}else{
		p <- myMap + geom_point(aes(longitude, latitude, size = withinBoundRatio, color = rmse), data = res)
	}
	p <- p + scale_color_gradientn(colors = terrain.colors(10))
	print(p)
	#p_rmse <-  myMap + geom_point(aes(longitude, latitude, color = rmse), data = res)
	#p_inBound <- myMap + geom_point(aes(longitude, latitude, color = withinBoundRatio), data = res)
}


RMSE <- function(x,y){
	return(round(mean((x-y)^2),digits = 3))
}

plot_cv <- function(year, method, aggType, cv_loggerID, r = NULL){
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



#main_CV(year = 2014, aggType = "daily")
#main_CV(year = 2014, aggType = "hourly")

#main_CV(year = 2015, aggType = "daily")
#main_CV(year = 2015, aggType = "hourly")


resultSummary()

# basisModel_REML <- basisModel(erieDO, trend, "Reml", 10,  "tmp")
# trendModel <- basisModel_REML$model
# residualModel <- basisModel_REML$residuals %>% idwModel()

# residualPredictions <- predict(residualModel, grid, TRUE) %>% reConstruct()
# trendPrediction <- predict(trendModel, grid) %>% reConstruct(residualPredictions,-1)

# tmp <- predict(trendModel, grid)
# trendPrediction <- reConstruct(tmp, residualPredictions, 0)
# hypoxiaExtent <- summary(tmp, residualPredictions, parallel = TRUE, totalSim = 10, randomSeed = 0)




