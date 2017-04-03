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
	year <- 2014
	aggType <- "daily"
	r <- 5

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

main_CV <- function(year = 2014, aggTye = "daily", ...){
	trend <- ~coords[,"x"]+coords[,"y"]+bathymetry+I(bathymetry^2)
	erieDO <- getLakeDO(year, "B", aggType)

	for(i in 1:nrow(erieDO$loggerInfo)) {
		cv_loggerID <- erieDO$loggerInfo$loggerID[i]
		subErieDO <- subset(erieDO, cv_loggerID)
		grid <- erieDO$loggerInfo[i,]
		
		hypoxia_idw <- predict(subErieDO, grid, "idw")
		
		for(r in c(5,10,15)){
			metaFolder <- sprintf("%s%d_%s_Reml_%d/cv/%s/",outputBaseName, year, aggType, r, cv_loggerID))
			hypoxia_reml_cv <- predict(obj =subErieDO, 
						grid = grid, 
						method = "Reml", 
						predictType = "simulations", 
						trend = trend, 
						r = r, 
						totalSim = 1000,
						metaFolder = metaFolder)

			saveRDS(hypoxia_reml_cv, sprintf("%s/simulations.rds",metaFolder))

			metaFolder <- sprintf("%s%d_%s_Baye_%d/cv/%s/",outputBaseName, year, aggType, r, cv_loggerID))
			hypoxia_baye_cv <- predict(obj =subErieDO, 
						grid = grid, 
						method = "Baye", 
						predictType = "simulations", 
						trend = trend, 
						r = r, 
						totalSim = 1000,
						metaFolder = metaFolder)

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










# basisModel_REML <- basisModel(erieDO, trend, "Reml", 10,  "tmp")
# trendModel <- basisModel_REML$model
# residualModel <- basisModel_REML$residuals %>% idwModel()

# residualPredictions <- predict(residualModel, grid, TRUE) %>% reConstruct()
# trendPrediction <- predict(trendModel, grid) %>% reConstruct(residualPredictions,-1)

# tmp <- predict(trendModel, grid)
# trendPrediction <- reConstruct(tmp, residualPredictions, 0)
# hypoxiaExtent <- summary(tmp, residualPredictions, parallel = TRUE, totalSim = 10, randomSeed = 0)




