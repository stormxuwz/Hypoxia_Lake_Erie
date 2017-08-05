# This is the script to represent data with basis functions
# Available basis: B-Spline and temporal basis function through SVD

source("./src/database.R")
source("./src/plot.R")
source("./src/helper.R")
#source("./src/interpolation.R")


require(sp)
require(gstat)
require(dplyr)
require(reshape2)
require(geoR)


coeff2Value <- function(coef, basis){
	return(basis %*% t(coef)) # return a matrix
}

plot_variogram <- function(df, formu = "value~1"){
	coordinates(df) = ~x+y
	print(plot(variogram(as.formula(formu),data =df,cutoff = 120, cloud=TRUE)))
}

NMF_scale <- function(NMFRes, center = "basis"){
	W <- NMFRes$basis
	H <- NMFRes$coef

	for(i in 1:ncol(W)){
		if(center == "basis"){
			print("conducting NMF and scale basis max to 12")
			alpha <- 12 / max(W[,i])
			W[,i] <- W[,i]*alpha
			H[i,] <- H[i,]/alpha
		}else if(center == "basis2"){
			print("conducting NMF and scale basis mean to 6")
			alpha <- 6 / mean(W[,i])
			W[,i] <- W[,i]*alpha
			H[i,] <- H[i,]/alpha
		}

		else if(center == "coef"){
			print("conducting NMF and scale max coefficients to 1")
			alpha <- 1 / max(H[i,])
			# alpha <- 1 / max(H)
			H[i,] <- H[i,] * alpha
			W[,i] <- W[,i] / alpha
		}
	}
	return(list(basis = W, coef = H))
}

NMF_basis <- function(DOdata, r, ...){
	method <- list(...)$method

	require(NMF)
	DOdata <- as.matrix(DOdata)

	if(is.null(method)){
		print("using default NMF fitting method")
		nmfRes <- nmf(DOdata, r, nrun = 60)
	}else{
		print(paste0("using ",method, " fitting gmethod"))
		nmfRes <- nmf(DOdata, r, nrun = 60, method = method, beta = 0.1)
	}
	

	W <- nmfRes@fit@W
	H <- nmfRes@fit@H


	return(list(basis = W, coef = H))
}


NMFReconstruct <- function(year, aggType, basis_r, loggerID, new = FALSE,...){
	erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()
	timeIdx <- index(erieDO$samplingData)

	method <- list(...)$method
	
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
	reconstructed <- data.frame(nmfDecomp$basis %*% nmfDecomp$coef)
	names(reconstructed) <- colnames(nmfDecomp$coef)
	
	res <- data.frame(nmf = reconstructed[,as.character(loggerID)], y = as.numeric(erieDO$samplingData[,loggerID]), time = timeIdx)
	print(qplot(timeIdx, nmf, data = res) + geom_line(aes(timeIdx, y), data = res, color = "red"))
}

NMF_analysis <- function(year, aggType, basis_r = NULL,explore = FALSE, new = FALSE, ...){
	erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()
	timeIdx <- index(erieDO$samplingData)
	baseMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",year)) + labs(x = "Longitude", y = "Latitude")
		 
	if(is.null(basis_r)){
	  if(year == 2015){
	    basis_r = 5
	  }else if(year == 2014){
	    basis_r = 5
	  }else if(year == 2016){
	    basis_r = 6
	  }
	}
	
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
		exploreRank <- nmfEstimateRank(as.matrix(d), c(2:10), method = method)
		saveRDS(exploreRank,sprintf("%s/results/cluster/nmfEstimateRank_%d_%s_%s.rds",outputBaseName, year,aggType, methodAlias))
	}
	# exploreRank <- readRDS(sprintf("%s/results/cluster/nmfEstimateRank_%d_%s.rds",outputBaseName, year,aggType))
	
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



SVD_basis <- function(DOdata, r){
	# r is the column vectors to keep
	DOdata <- as.matrix(DOdata)

	print("Doing Scaling!!")
	svdRes <- svd(DOdata %>% scale()) # is centered necessary? No (3/23) // Yes (the original paper did so)


	basis <- svdRes$u[,1:r]
	t = nrow(basis)

	for(i in 1:r){
		splineRes <- smooth.spline(x = 1:t, y = basis[,i], df = t*0.2)
		basis[,i] <- predict(splineRes, 1:t)$y
	}

	basis <- cbind(basis,1) # add the bias term

	coef <- lsfit(x = basis, y= DOdata, intercept = FALSE)$coefficients
	DO_fit <- basis %*% coef # coef, each column is the coefficents for different basis
	
	print("variance explained")
	print(sum(svdRes$d[1:r]**2)/sum(svdRes$d**2))
	return(list(fit = DO_fit, coef = coef, basis = basis,varExpl = sum(svdRes$d[1:r]**2)/sum(svdRes$d**2)))
}


B_spline <- function(DOdata,knots,norder = 4){
	require(fda)
	T <- nrow(DOdata)
	DOdata <- as.matrix(DOdata)
	nBasis <- length(knots)-2+norder
	bbasis <- create.bspline.basis(range(knots),nBasis,norder,knots)

	basismat <- eval.basis(1:T, bbasis)
	coef <- solve(crossprod(basismat), crossprod(basismat,DOdata))
	
	DOfd <- fd(coef, bbasis)
	DO_fit <- predict(DOfd,1:T)

	return(list(fit = DO_fit, coef = coef, basis = bbasis))  # return fit value, coefficients and basis object
}	
