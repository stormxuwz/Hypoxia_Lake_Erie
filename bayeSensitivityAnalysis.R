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


bayeSensitivity <- function(year = 2014, aggType = "daily", cv = FALSE){
	trend <- ~coords[,"x"]+ coords[,"y"] + bathymetry + I(bathymetry^2)
	erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()
	grid <- createGrid(erieDO$loggerInfo, mapDx, mapDy) # grid size: 0.025, 0.025 long,lat
	
	for(r in rList){
		print(paste("r:",r))
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

	if(cv){
		for(i in 1:nrow(erieDO$loggerInfo)) {
			cv_loggerID <- erieDO$loggerInfo$loggerID[i]
			subErieDO <- subset(erieDO, cv_loggerID)  #remove cv_loggerID from erieDO
			grid <- erieDO$loggerInfo[i,]
			grid$convexIndex <- 1

			trueDO <- erieDO$samplingData[,cv_loggerID]
			timeIdx <- index(trueDO)		

		for(r in rList){
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
	

}


library(geoR)

myPriorList <- list(

	# effects of phi.prior
	# change phi.prior to reciprocal
	prior1 = prior.control(
			phi.prior = "reciprocal",
			phi.discrete=seq(20,70,5),
			beta.prior = "flat",  # beta is flat
			sigmasq.prior = "reciprocal",
			tausq.rel.prior = "fixed",
			tausq.rel = 0),

	# change phi.prior to uniform
	prior2 = prior.control(
			phi.prior = "uniform",
			# phi.discrete=seq(20,70,5),
			beta.prior = "flat",  # beta is flat
			sigmasq.prior = "reciprocal",
			tausq.rel.prior = "fixed",
			tausq.rel = 0),

	# change phi.prior to squared.reciprocal
	prior3 = prior.control(
		phi.prior = "squared.reciprocal",
		phi.discrete=seq(20,70,5),
		beta.prior = "flat", 
		sigmasq.prior = "reciprocal",
		tausq.rel.prior = "fixed",
		tausq.rel = 0),


	# effects of tarsq.rel
	# change tausq.rel.prior to 0.1
	prior4 = prior.control(
		phi.discrete=seq(20,70,5),
		beta.prior = "flat",
		sigmasq.prior = "reciprocal",
		tausq.rel.prior = "fixed",
		tausq.rel = 0.1),

	# change tausq.rel.prior to 0.2
	prior5 = prior.control(
		phi.discrete=seq(20,70,5),
		beta.prior = "flat", 
		sigmasq.prior = "reciprocal",
		tausq.rel.prior = "fixed",
		tausq.rel = 0.2),

	# change tausq.rel.prior to uniform
	prior6 = prior.control(
		phi.discrete=seq(20,70,5),
		beta.prior = "flat",
		sigmasq.prior = "reciprocal",
		tausq.rel.prior = "uniform",
		tausq.rel.discrete = seq(0, 1, 0.1)),

	# change tausq.rel.prior to reciprocal
	prior7 = prior.control(
		phi.discrete=seq(20,70,5),
		beta.prior = "flat",
		sigmasq.prior = "reciprocal",
		tausq.rel.prior = "reciprocal",
		tausq.rel.discrete = seq(0, 1, 0.1),
		tausq.rel = 0),


	# change cov.model to spherical
	prior8 = "spherical.cov",

	# change sigma to sc.inv.chisq with df  = 1
	prior9 = "df_1_sc.inv.chisq.sigmasq",
	
	# change sigma to sc.inv.chisq with df  = 3
	prior10 = "df_3_sc.inv.chisq.sigmasq"

	# change sigma
	# prior11 = "fix.sigmasq"	
)

mapDx <- mapDy <- 0.025
# rList <- c(2,5,10,15)
rList <- c(10)
# mapDx <- mapDy <- 0.1

# for(i in 1:10){
# 	print(paste("parameter set", i))
# 	outputBaseName <- paste0("/Users/wenzhaoxu/Developer/Hypoxia/output_sensitivity_", i,"/")
# 	# myPrior <- myPriorList[[paste0("prior",i)]]
# 	# print(system.time(bayeSensitivity(year = 2014, aggType = "hourly", cv = TRUE)))
# 	resultSummary(aggList = c("hourly"), yearList = c(2014), methodList = c("Baye"))
# }

# read CV results
fullRes <- data.frame()
for(i in 1:10){
	print(paste("parameter set", i))
	outputBaseName <- paste0("/Users/wenzhaoxu/Developer/Hypoxia/output_sensitivity_", i,"/")
	
	tmpRes <- readRDS(sprintf("%s/results/fullRes.rds", outputBaseName))
	tmpRes$parameter <- i

	p <- plotCVOnMap(year_ = 2014, aggType_ = "hourly", method_ = "Baye", r_ = 10)
	filePrefix <- sprintf("%d_%s_%s_%d",2014, "hourly", "Baye", 10)
	pdf(sprintf("%s/results/%s_CV_summary.pdf",outputBaseName, filePrefix), width = 6, height = 4)
	print(p)
	dev.off()

	fullRes <- rbind(fullRes, tmpRes)
}

# read original

baseRes <- readRDS("/Users/wenzhaoxu/Developer/Hypoxia/output/results/fullRes.rds") %>% 
	dplyr::filter(year == 2014, method == "Baye", aggType == "hourly", r == 10) %>% 
	dplyr::mutate(parameter = 0)

diffRes <- merge(fullRes, baseRes, by.x = "cv_loggerID", by.y = "cv_loggerID", all.x = T) %>% 
	mutate(rmse_diff = rmse.x - rmse.y, CI_coverage_diff = withinBoundRatio.x - withinBoundRatio.y) %>% # altervative - base, negative rmse or positive CI means altervative better
	rename(parameter = parameter.x)
diffRes$parameter <- as.factor(diffRes$parameter)

pdf("/Users/wenzhaoxu/Developer/Hypoxia/bayeSensitivity_dRMSE.pdf", width = 5.8, height = 2.5)
ggplot(data = diffRes) + geom_boxplot(aes(x = as.factor(parameter), y = rmse_diff),size = I(0.5), position = position_dodge(width = 0.8),outlier.size = 0.5) + 
xlab("Hyper parameter set") + ylab(TeX("$\\Delta$ RMSE")) + theme_bw()
dev.off()

pdf("/Users/wenzhaoxu/Developer/Hypoxia/bayeSensitivity_dCICover.pdf", width = 5.8, height = 2.5)
ggplot(data = diffRes) + geom_boxplot(aes(x = as.factor(parameter), y = CI_coverage_diff),size = I(0.5), position = position_dodge(width = 0.8),outlier.size = 0.5) + 
xlab("Hyper parameter Set") + ylab(TeX("$\\Delta$ CI coverage")) + theme_bw()
dev.off()



ggplot(data = diffRes) + geom_boxplot(aes(x = as.factor(parameter), y = rmse_diff),size = I(0.5), position = position_dodge(width = 0.8),outlier.size = 0.5) + 
xlab("Hyper parameter set") + ylab(TeX("$\\Delta$ RMSE")) + theme_bw() + geom_text(aes(parameter, rmse_diff, label = cv_loggerID))

ggplot(data = diffRes) + geom_boxplot(aes(x = as.factor(parameter), y = CI_coverage_diff),size = I(0.5), position = position_dodge(width = 0.8),outlier.size = 0.5) + 
xlab("Hyper parameter set") + ylab(TeX("$\\Delta$ RMSE")) + theme_bw() + geom_text(aes(parameter, CI_coverage_diff, label = cv_loggerID))

fullRes$parameter <- as.factor(fullRes$parameter)
# saveRDS(fullRes,"/Users/wenzhaoxu/Developer/Hypoxia/bayeSensitivityRes.rds")
# fullRes <- readRDS("/Users/wenzhaoxu/Developer/Hypoxia/bayeSensitivityRes.rds")

# fullRes <- rbind(fullRes, baseRes)
# ggplot(data = fullRes) + geom_boxplot(aes(x = as.factor(parameter), y = rmse, fill = parameter),size = I(0.5), position = position_dodge(width = 0.8),outlier.size = 0.5) + 
# xlab("Bayesian Prior Parameter Set") + ylab("RMSE") + theme_bw()

# ggplot(data = fullRes) + geom_boxplot(aes(x = as.factor(parameter), y = withinBoundRatio, fill = parameter),size = I(0.5), position = position_dodge(width = 0.8),outlier.size = 0.5) + 
# xlab("Bayesian Prior Parameter Set") + ylab("RMSE") + theme_bw()

fullRes <- readRDS(sprintf("%s/results/fullRes.rds",outputBaseName))
EPA <- dplyr::filter(fullRes,site %in% EPASite, aggType == "hourly", r == 10, method == "Baye")
ggplot(data = EPA) + geom_bar(aes(x = site, y = rmse, fill = factor(year)), stat = "identity", position=position_dodge()) + 
theme_bw()

ggplot(data = EPA) + geom_bar(aes(x = site, y = withinBoundRatio, fill = factor(year)), stat = "identity", position=position_dodge()) + 
theme_bw()