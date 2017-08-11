require(fields)
require(gstat)
require(sp)
require(dplyr)
require(geoR)

predict.krigModelIdw <- function(model, grid, parallel){
	# model is a list of data frames with x,y and value
	# grid is a data frame with x,y and other covariates
	metaFolder <- attr(model, "metaFolder")
	convexIndex <- grid$convexIndex
	
	config <- attr(model, "config")
	trend <- as.formula(config$trend) # "value ~ 1"
	simNum <- config$simNum  # how many simulations for this prediction


	grid0 <- grid
	predictions <- list()
	pred <- matrix(NA, nrow = nrow(grid), ncol = simNum)

	grid <- subset(grid, convexIndex == 1)
	coordinates(grid) = ~x + y
	
	if(!parallel){
		for(i in 1:length(model)){
			df <- model[[i]]
			coordinates(df) = ~x + y
			pred[convexIndex == 1,1] <- idw(trend , df, grid, nmax = config$nmax)$var1.pred
			predictions[[i]] <- list(pred = pred)
		}
	}else{
		require(doParallel)
		cl <- makeCluster(2)
		registerDoParallel(cl)

		predictions <- foreach(i=1:length(model)) %dopar% {
			require(fields)
			require(gstat)
			require(sp)
			
			df <- model[[i]]
			coordinates(df) = ~x + y
			pred[convexIndex == 1,1] <- idw(trend , df, grid, nmax = config$nmax)$var1.pred
			list(pred = pred)
		}
		stopCluster(cl)
		gc()
	}
	
	res <- list(predictions = predictions, grid = grid0)
	class(res) <- "idwModel"

	return(res)
}


predict.krigModelReml <- function(model, grid){
	config <- attr(model, "config")
	simNum <- config$simNum
	trend <- config$trend
	metaFolder <- attr(model, "metaFolder")
	convexIndex <- grid$convexIndex

	# x should be a data frame with x,y and value
	predSimulation <- matrix(NA,nrow = nrow(grid), ncol = simNum) # n_grid * simNum
	pred <- matrix(NA, nrow = nrow(grid),ncol = 2) 
	# one for prediction, one for variance

	# start predictions
	predictions <- list()
	grid0 <- grid
	grid <- subset(grid,convexIndex==1)
	
	trend.l <- dplyr::select(grid, x,y,bathymetry) %>%  
			as.geodata(covar.col = 3) %>% 
			trend.spatial(trend,.)

	for(i in 1:length(model)){
		if(i == 6){
			print("i=6")
		}
			
		df <- model[[i]] %>% 
			dplyr::select(x, y, value,bathymetry) %>%
			as.geodata(covar.col = c("bathymetry"))

		trend.d <- trend.spatial(trend,df)
		
		semiVariance <- variog(df,trend = trend.d)

		ml <- likfit(df, 
			ini = c(max(semiVariance$v),70),
			fix.nugget = T,  # the nugget is defaulted as zero
			lik.method = "ML",
			cov.model = config$cov.model, # "exponential",
			trend = trend.d,
			fix.psiA = TRUE,
			fix.psiR = TRUE)
			# psiA = pi/4)
			# limits = pars.limits(psiR = c(lower = 1, upper = 10)))
		print(ml)
		if(!is.null(metaFolder)){
			print("saving prediction meta to metaFolder")
			png(paste0(metaFolder,"basis_",i,"_data.png"))
			print(plot(df))
			dev.off()

			png(paste0(metaFolder,"basis_",i, "_reml_vgm",".png"))
			par(mfrow = c(1,2))
			print(plot(variog(df,trend = trend.d,option = "cloud",direction = pi/4)))
			print(lines(ml))
			print(plot(variog(df,trend = trend.d,option = "cloud",direction = 3*pi/4)))
			print(lines(ml))
			dev.off()
		}
		

		modelPred <- krige.conv(df, 
				locations = grid[,c("x","y")], 
				krige = krige.control(obj.m = ml, trend.d = trend.d, trend.l = trend.l),
				output = output.control(n.predictive = simNum))

		predSimulation[convexIndex == 1, ] <- modelPred$simulations
		pred[convexIndex == 1, 1] <- modelPred$predict
		pred[convexIndex == 1, 2] <- modelPred$krige.var

		predictions[[i]] <- list(simulations = predSimulation, pred = pred)

	}
	res <- list(predictions = predictions, grid = grid0, basis = attr(model, "basis"))
	class(res) <- "basisModel"
	return(res)
}





predict.krigModelBaye <- function(model, grid, defaultPrior = FALSE){
	config <- attr(model, "config")
	simNum <- config$simNum
	trend <- config$trend
	convexIndex <- grid$convexIndex
	metaFolder <- attr(model, "metaFolder")

	grid0 <- grid
	predSimulation <- matrix(NA,nrow = nrow(grid), ncol = simNum) # n_grid * simNum
	pred <- matrix(NA, nrow = nrow(grid),ncol = 2)  # one for prediction, one for variance

	# start predictions
	predictions <- list()
	grid <- subset(grid,convexIndex==1)

	trend.l <- dplyr::select(grid, x,y,bathymetry) %>%  
		as.geodata(covar.col = 3) %>%
		trend.spatial(trend,.)
	
	# specify the priors
	PC <- prior.control(
						phi.discrete=seq(20,70,5),  # range is discreted
						beta.prior = "flat",  # beta is flat 
						sigmasq.prior = "reciprocal",
						tausq.rel.prior = "fixed",
						tausq.rel = 0)  # sigma^2 is 

	OC <- output.control(n.pos = simNum, # the number of samples taking from posterior distribution
						n.pred = simNum,   # sample to taken from the predictive distribution
						 signal = FALSE)

	for(i in 1:length(model)){

		df <- model[[i]] %>% 
			dplyr::select(x, y, value, bathymetry) %>%
			as.geodata(covar.col = 4)

		trend.d <- trend.spatial(trend,df)
		
		# change prior distribution if indicated
		if(!defaultPrior){
			print ("using non-default prior")
			if(myPrior == "spherical.cov"){
				print("change to sperical cov.model")
				config$cov.model = "spherical"
			}
			else if(myPrior == "df_1_sc.inv.chisq.sigmasq" | myPrior == "df_3_sc.inv.chisq.sigmasq"){
				semiVariance <- variog(df,trend = trend.d)
				print("change to sc.inv.chisq.sigmasq")
				df.sigmasq <-  as.numeric(substr(myPrior,4,4))

				ml <- likfit(df, 
				ini = c(max(semiVariance$v),70),
				fix.nugget = T,  # the nugget is defaulted as zero
				lik.method = "ML",
				cov.model = config$cov.model, # "exponential",
				trend = trend.d,
				fix.psiA = TRUE,
				fix.psiR = TRUE)

				PC <- prior.control(
						phi.discrete=seq(20,70,5),  # range is discreted
						beta.prior = "flat",  # beta is flat 
						sigmasq.prior = "sc.inv.chisq",
						sigmasq = ml$sigmasq,
						df.sigmasq = df.sigmasq,
						tausq.rel.prior = "fixed",
						tausq.rel = 0) 

			# }else if(myPrior == "fix.sigmasq"){
			# 	print("change to fix sigmasq")
			# 	semiVariance <- variog(df,trend = trend.d)
			# 	ml <- likfit(df, 
			# 		ini = c(max(semiVariance$v),70),
			# 		fix.nugget = T,  # the nugget is defaulted as zero
			# 		lik.method = "ML",
			# 		cov.model = config$cov.model, # "exponential",
			# 		trend = trend.d,
			# 		fix.psiA = TRUE,
			# 		fix.psiR = TRUE)
			# 	PC <- prior.control(
			# 			# phi.discrete=seq(20,70,5),  # range is discreted
			# 			beta.prior = "flat",  # beta is flat 
			# 			sigmasq.prior = "fixed",
			# 			sigmasq = ml$sigmasq,
			# 			tausq.rel.prior = "fixed",
			# 			tausq.rel = 0) 
			# }
			}else{
				print("change to other prior distribution")
				PC = myPrior
			}
		}

		MC <- model.control(
			cov.model = config$cov.model, 
			trend.l = trend.l, 
			trend.d = trend.d) 

		modelRes <- krige.bayes(geodata = df, 
								 locations = grid[,c("x","y")], 
								 model = MC,
								 prior = PC,
								 output = OC)
    
		modelPred <- modelRes$predictive
		
		if(!is.null(metaFolder)){
			print("saving prediction meta to metaFolder")
			png(paste0(metaFolder,"basis_",i,"_data.png"))
			print(plot(df))
			dev.off()
			variogramSummary <- function(x){quantile(x, prob = c(0.05, 0.5, 0.95))}

			png(paste0(metaFolder,"basis_",i, "_bayes_vgm",".png"))
			print(plot(variog(df,option = "cloud",trend = "1st")))
			print(lines(modelRes, summ = variogramSummary, ty="l", lty=c(2,1,2), col=1))
			dev.off()
		}
		
		predSimulation[convexIndex == 1, ] <- modelPred$simulations
		pred[convexIndex == 1, 1] <- modelPred$mean
		pred[convexIndex == 1, 2] <- modelPred$variance

		predictions[[i]] <- list(simulations = predSimulation, 
			pred = pred)
	}
	res <- list(predictions = predictions, grid = grid, basis = attr(model, "basis"))
	
	class(res) <- "basisModel"
	return(res)
}



predict.lakeDO <- function(obj, grid, method, predictType, ...){
	# obj is the lakeDO object
	metaFolder = list(...)$metaFolder
	nmax <- list(...)$nmax
	
	if(!is.null(metaFolder)) createFolder(metaFolder)
	
	if(!method %in% c("idw","Reml","Baye")){
		stop("interpolation method not implemented")
	}

	if(!predictType %in% c("extent","expected","simulations")){
		stop("predictType not implemented")
	}


	if(method == "idw"){
		if(is.null(nmax)){
			stop("nmax is not specify")
		}

		res <- obj %>% 
			idwModel(metaFolder = metaFolder, nmax = nmax) %>% 
			predict(grid, parallel = TRUE)

		if(!is.null(metaFolder)) saveRDS(res, paste0(metaFolder,"trendPredictions.rds"))
		
		if(predictType == "extent"){
		 	res <- summary(res)
		}else if(predictType == "simulations"){
			res <- reConstruct(res)
		}else{
			stop("predictType not implemented")
		}

	}else if(method %in% c("Reml","Baye")){

		trend = list(...)$trend
		r = list(...)$r
		totalSim <- list(...)$totalSim
		randomSeed <- list(...)$randomSeed

		basisModelRes <- obj %>% 
			basisModel(trend, method, r, metaFolder) 
		
		# interpolate residuals
		print("interpolating residuals")
		residualPredictions <- basisModelRes$residuals %>% 
				idwModel(metaFolder = metaFolder, nmax = nmax) %>% 
				predict(grid, parallel = TRUE) %>% reConstruct()	
		
		
		# predict trend model
		print("interpolating trend")
		trendModel  <- basisModelRes$model %>% predict(grid)
		availableSimNum <- ncol(trendModel$predictions[[1]]$simulations)

		if(!is.null(metaFolder)){
			saveRDS(basisModelRes, paste0(metaFolder,"basisModelRes.rds"))
			saveRDS(residualPredictions, paste0(metaFolder,"residualPredictions.rds"))
			saveRDS(trendModel, paste0(metaFolder,"trendModel.rds"))
		}

		basisNum <- r + 1

		set.seed(randomSeed)
		indMatrix <- base::sample(1:availableSimNum, totalSim*basisNum,replace= TRUE) %>%
				matrix(nrow = basisNum)
		
		if(predictType == "extent"){
			res <- summary(
				trendModel,
				residualPredictions, 
				parallel = TRUE, 
				totalSim = totalSim, 
				indMatrix = indMatrix)

		}else if(predictType == "expected"){
			# return all values
			res <- reConstruct(
				trendModel, 
				residualPredictions,
				simulationNum = 0,
				indMatrix = NULL)

		}else if(predictType == "simulations"){
			res <- reConstruct(
				trendModel,
				residualPredictions,
				simulationNum = -totalSim,
				indMatrix = indMatrix,
				parallel = TRUE)
		}else{
			stop("predictType not implemented")
		}
	}else{
		stop("method is not implemented")
	}
	return(res)
}