# backend of interpolation 

spatial_interpolation <- function(df,grid,method = "IDW", simNum = 1,tmpName = "tmp",...){
	# df and grid is a dataframe that contains longitude and latitude and value as columns
	pred <- array(NA,dim=c(nrow(grid),simNum)) # shapes as the [gridNum, simNum]
	convexIndex <- grid$convexIndex
	grid <- subset(grid,convexIndex==1)
	
	# check the input validation
	if(method=="IDW" & simNum!=1){
		throw("IDW shouldn't should have simNum =1")
	}

	if(simNum<1){
		throw("Sim num should be positive")
	}

	if(method == "IDW"){
		coordinates(df) = ~x + y
		coordinates(grid) = ~x + y
		pred[convexIndex == 1,1] <- idw(value ~  1 , df, grid, nmax = 5)$var1.pred    # nearest 5 points to do IDW
	}
	else if(method == "loglik"){
		# Using log likelihood to fit covariance
		df <- df[,c("x","y","value","bathymetry")] %>% as.geodata(covar.col = 4)
		# trend <- as.formula(~bathymetry + I(bathymetry^2))
		# trend <- "cte"
		ml <- likfit(df, ini = c(1000,70), fix.nugget = T, lik.method = "REML",cov.model = "spherical",trend = trend)  # the nugget is defaulted as zero
		
		# plot the variogram and data
		png(paste(metaFolder,tmpName,"loglik_df.png",sep = ""))
		print(plot(df))
		dev.off()
		
		png(paste(metaFolder,tmpName, "loglik_variogram_",trend,".png",sep=""))
		print(plot(variog(df,trend = trend,option = "cloud")))
		print(lines(ml))
		dev.off()

		# kriging control
		KC <- krige.control(obj.m = ml, trend.d = trend, trend.l = trend)
		
		if(simNum ==1){
			pred[convexIndex == 1, 1] <- krige.conv(df, locations = grid[,c("x","y")], krige = KC)$predict
		}else{
			# do conditional simulation
			pred[convexIndex == 1, ] <- krige.conv(df, locations = grid[,c("x","y")], krige = KC, output = output.control(n.predictive = simNum))$simulations
		}
		
	}
	else if(method == "baye"){
		# conduct the bayesian kriging
		#trend = "1st"
		MC <- model.control(cov.model = "exponential", trend.l = trend, trend.d = trend) 
		
		# specify the priors
		PC <- prior.control(phi.discrete=seq(20,60,5),  # range is discreted
							beta.prior = "flat",  # beta is flat 
							sigmasq.prior = "reciprocal",
							tausq.rel.prior = "fixed",
							tausq.rel = 0)  # sigma^2 is 
		 
		# specify the output control
		OC <- output.control(n.pos = simNum, # the number of samples taking from posterior distribution
							n.pred = simNum,   # sample to taken from the predictive distribution
							 signal = FALSE)
		
		df <- df[,c("x","y","value","bathymetry")] %>% as.geodata(covar.col = 4)
		
		predRes <- krige.bayes(geodata = df, 
								 locations = grid[,c("x","y")], 
								 model = MC,
								 prior = PC,
								 output = OC)
		pred[convexIndex == 1,] <- predRes$predictive$simulations
		
		png(paste(metaFolder,tmpName,"baye_df.png",sep = ""))
		print(plot(df))
		dev.off()
		mySummary <- function(x){quantile(x, prob = c(0.05, 0.5, 0.95))}
		
		png(paste(metaFolder,tmpName, "baye_variogram_",trend,".png",sep=""))
		print(plot(variog(df,option = "cloud",trend = "1st")))
		print(lines(predRes, summ = mySummary, ty="l", lty=c(2,1,2), col=1))
		dev.off()
		
	}
	else if(method == "rand"){
		# give random interpolation
	}
	else{
		throw()
	}
	return(pred)
}


basis_interpolation <- function(DOdata, logger_geo, grid, timeIndex = NULL, basis = "fda",intMethod = "likfit", simNum = 1,...){
	# function to decompose the data into basis and then do interpolation
	# DOdata is the dataframe
	# basis can take "fda" or "svd"
	# time index is a vector that specify how long data is used

	args <- list(...)
	saveMeta <- args$saveMeta
	print(saveMeta)
	if(is.null(saveMeta))
		saveMeta <- TRUE
	ID <- as.numeric(colnames(DOdata))
	TimeN <- nrow(DOdata)  # Totally how many timestamps
	
	if(basis == "fda"){
		decompRes <- B_spline(DOdata.args$knots)
	}else{
		decompRes <- SVD_basis(DOdata,args$r)
		nBasis <- dim(decompRes$coef)[1] # columns of the coefficent matrix
		varExpl <- decompRes$varExpl
	}
	print(sprintf("total number of basis: %d",nBasis))
	basis <- decompRes$basis
	
	if(saveMeta)
		saveRDS(basis, file = paste(metaFolder,"basis.rds",sep =""))
	
	# construct complete coefficient matrix
	coef <- decompRes$coef %>% t() %>% as.data.frame()

	coef$ID <- ID
	coef_df <- merge(coef,logger_geo,by.x = "ID",by.y = "loggerID")
	
	if(saveMeta)
		saveRDS(coef_df, file = paste(metaFolder,"coef_df.rds",sep = ""))

	DO_res <- DOdata - decompRes$fit  # a matrix of shape T, nLocation
	
	if(saveMeta)
		saveRDS(DO_res, file = paste(metaFolder,"DO_res.rds", sep = ""))
	# do interpolation on each basis coefficients
	# prediction <- matrix(0,T,nrow(grid)) # a matrix T * nBasis, row is time, columns is locations
	prediction <- array(0, dim = c(simNum, TimeN, nrow(grid)))
	
	for(i in 1:nBasis){
		print(sprintf("doing basis %d",i))

		spData <- coef_df[,c("x","y",paste("X",i,sep = ""),"bathymetry","longitude","latitude")]
		names(spData)[3] <- "value"
		
		pred <- spatial_interpolation(df=spData,grid=grid,method=intMethod, simNum =simNum,tmpName = paste("basis_",i,sep=""))
		# pred is matrix with [# of grid, numSim]
		
		if(saveMeta)
			saveRDS(pred, file = sprintf("%spred_coef_%d.rds",metaFolder,i))

		for(sim in 1:dim(pred)[2]){
			prediction[sim,,] <- prediction[sim,,]+basis[,i] %*% t(pred[,sim])
		}
		rm(pred)
	}
	

	return(list(trend = prediction,res = DO_res,varExpl = varExpl))
}


random_interpolation <- function(data, locationInfo, nsim = 100){
	# function to interpolate the map based on random values. i.e. no correlation
	grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)
	data <- as.numeric(data)
	for(i in 1:nsim){}
}
