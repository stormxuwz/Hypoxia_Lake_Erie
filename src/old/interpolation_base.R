# backend of interpolation 

spatial_interpolation <- function(df,grid,method = "IDW", simNum = 1,tmpName = "tmp",...){
	# function to do spatial interpolation given df and grid
	# df and grid are dataframes that contains longitude and latitude and value as columns
	pred <- array(NA,dim=c(nrow(grid),simNum)) # shapes as the [gridNum, simNum]
	convexIndex <- grid$convexIndex
	grid <- subset(grid,convexIndex==1)
	
	# normalize df and grid coordinates
	# scalePara <-  scale(df[,c("x","y","bathymetry")])
	# df[,c("x","y","bathymetry")] <- scalePara
	# scaleCenter <- attr(scalePara,"scaled:center")
	# scaleVar <- attr(scalePara,"scaled:scale")
	# grid[,c("x","y","bathymetry")] <- scale(grid[,c("x","y","bathymetry")],center = scaleCenter, scale = scaleVar)

	# check the input validation
	if(method=="IDW" & simNum!=1)
		throw("IDW should have simNum =1")

	if(simNum<1)
		throw("Sim num should be positive")

	if(method == "IDW"){
		print("Using IDW")
		coordinates(df) = ~x + y
		coordinates(grid) = ~x + y
		pred[convexIndex == 1,1] <- idw(value ~  1 , df, grid, nmax = 5)$var1.pred    
		# nearest 5 points to do IDW
	}
	else if(method == "localKrig"){
		require(gstat)
		coordinates(df) = ~x+y
		coordinates(grid) = ~x+y
		model_gstat <- gstat(NULL,id = "value", formula = value~x+y+bathymetry+I(bathymetry^2), data = df, nmax = 5, nmin = 1)
		v <- variogram(model_gstat, cressie=T)
	}


	else if(method == "loglik"){
		print("Using log likelihood to fit variogram")
		
		# plot df
		png(paste(metaFolder,tmpName,"df.png",sep = ""))
		newDf <- df
		newDf$lmRes <- lm(value~x+y+bathymetry+I(bathymetry^2), data = newDf)$residuals
		print(qplot(x,y,data = newDf,color = lmRes,size = I(10))+
			scale_color_gradientn(colours = terrain.colors(10))+
			coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE))
		dev.off()

		# Using log likelihood to fit covariance
		df <- df[,c("x","y","value","bathymetry")] %>% as.geodata(covar.col = 4)
		trend.d <- trend.spatial(trend,df)

		invisible(ml <- likfit(df, 
			ini = c(1000,70), 
			fix.nugget = T,  # the nugget is defaulted as zero
			lik.method = "REML",
			cov.model = "exponential",
			trend = trend.d,
			fix.psiA = TRUE, 
			fix.psiR = FALSE,
			psiA = pi/4))  

		trend.l <- trend.spatial(trend,geodata = as.geodata(grid[,c("x","y","bathymetry")], covar.col = 3))

		# plot the variogram and data
		png(paste(metaFolder,tmpName,"loglik_df.png",sep = ""))
		print(plot(df))
		dev.off()

		png(paste(metaFolder,tmpName, "loglik_variogram_",trend,".png",sep=""))
		par(mfrow = c(1,2))
		print(plot(variog(df,trend = trend.d,option = "cloud",direction = pi/4)))
		print(lines(ml))
		print(plot(variog(df,trend = trend.d,option = "cloud",direction = 3*pi/4)))
		print(lines(ml))
		dev.off()

		# kriging control
		KC <- krige.control(obj.m = ml, trend.d = trend.d, trend.l = trend.l)
		
		if(simNum ==1){
			pred[convexIndex == 1, 1] <- krige.conv(df, 
				locations = grid[,c("x","y")], 
				krige = KC)$predict
		}else{
			# do conditional simulation
			pred[convexIndex == 1, ] <- krige.conv(df, 
				locations = grid[,c("x","y")], 
				krige = KC,
				output = output.control(n.predictive = simNum))$simulations
		}
	}
	else if(method == "baye"){
		# conduct the bayesian kriging
		df <- df[,c("x","y","value","bathymetry")] %>% as.geodata(covar.col = 4)
		
		trend.d <- trend.spatial(trend,df)
		trend.l <- trend.spatial(trend,geodata = as.geodata(grid[,c("x","y","bathymetry")], covar.col = 3))												 
		
		MC <- model.control(cov.model = "exponential", trend.l = trend.l, trend.d = trend.d) 
		
		
		 
		# specify the output control
		OC <- output.control(n.pos = simNum, # the number of samples taking from posterior distribution
							n.pred = simNum,   # sample to taken from the predictive distribution
							 signal = FALSE)
		
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
	else{
		throw()
	}
	return(pred)
}


basis_interpolation <- function(DOdata, logger_geo, grid, timeIndex = NULL, basisDecomp = "fda",fitMethod = "likfit", simNum = 1,...){
	# function to decompose the data into basis and then do interpolation
	# DOdata is the dataframe
	# basisDecomp can take "fda" or "svd"
	# time index is a vector that specify how long data is used

	# return list(trendCoeff = coeffPredList,basis = basis, res = DO_res,varExpl = varExpl)

	args <- list(...)
	saveMeta <- args$saveMeta
	print(saveMeta)
	if(is.null(saveMeta))
		saveMeta <- TRUE
	ID <- as.numeric(colnames(DOdata))
	TimeN <- nrow(DOdata)  # Totally how many timestamps
	
	if(basisDecomp == "fda"){
		decompRes <- B_spline(DOdata.args$knots)
	}else{
		decompRes <- SVD_basis(DOdata,args$r) # return a list with "fit","coef","basis","varExpl"
		nBasis <- dim(decompRes$coef)[1] # columns of the coefficent matrix
		varExpl <- decompRes$varExpl
	}
	print(sprintf("total number of basis: %d",nBasis)) # nBasis should = args$r
	basis <- decompRes$basis # the basis matrix
	
	if(saveMeta)
		saveRDS(basis, file = paste(metaFolder,"basis.rds",sep =""))
	
	# construct complete coefficient matrix
	coef <- decompRes$coef %>% t() %>% as.data.frame()

	coef$ID <- ID
	coef_df <- merge(coef,logger_geo,by.x = "ID",by.y = "loggerID") # create the parameter data frame
	
	if(saveMeta)
		saveRDS(coef_df, file = paste(metaFolder,"coef_df.rds",sep = ""))

	DO_res <- DOdata - decompRes$fit  # a matrix of shape [T, nLocation]
	
	if(saveMeta)
		saveRDS(DO_res, file = paste(metaFolder,"DO_res.rds", sep = ""))
	
	# Do interpolation on each basis coefficients
	coeffPredList <- list()
	
	for(i in 1:nBasis){
		print(sprintf("doing basis %d",i))

		spData <- coef_df[,c("x","y",paste("X",i,sep = ""),"bathymetry","longitude","latitude")]
		names(spData)[3] <- "value"
		
		pred <- spatial_interpolation(
			df=spData,
			grid=grid,
			method=fitMethod,
			simNum =simNum,
			tmpName = paste("basis_",i,sep=""))

		coeffPredList[[i]] <- pred
	}
	
	return(list(trendCoeff = coeffPredList,basis = basis, res = DO_res,varExpl = varExpl))
}


constructFromCoeff <- function(coeffPred,basis,gridNum,simNum = 100){
	nBasis <- ncol(basis)
	# generate randomNum
	TimeN <- nrow(basis)
	prediction <- array(0.0, dim = c(simNum, TimeN, gridNum))

	for(i in 1:nBasis){
		pred <- coeffPred[[i]]
		coeff_simNum <- dim(pred)[2]
		sampleIndex <- sample(1:coeff_simNum,size =simNum,replace = TRUE)
		for(sim in 1:simNum){
			prediction[sim,,] <- prediction[sim,,] + basis[,i] %*% t(pred[,sampleIndex[sim]])
		}
	}
	return(prediction)
}


random_interpolation <- function(data, locationInfo, nsim = 100){
	# function to interpolate the map based on random values. i.e. no correlation
	grid <- createGrid(locationInfo,by.x = mapDx, by.y = mapDy)
	data <- as.numeric(data)
	for(i in 1:nsim){}
}
