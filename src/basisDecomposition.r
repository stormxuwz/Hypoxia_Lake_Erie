# This is the script to represent data with basis functions
# Available basis: B-Spline and temporal basis function through SVD

source("./src/database.R")
source("config.R")
source("./src/plot.R")
source("./src/spatialHelper.R")
#source("./src/interpolation.R")

require(fda)
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

SVD_basis <- function(DOdata, r){
	# r is the column vectors to keep
	DOdata <- as.matrix(DOdata)

	svdRes <- svd(DOdata)

	# coef <- (diag(svdRes$d) %*% t(svdRes$v))[1:r,]
	basis <- svdRes$u[,1:r]
	t = nrow(basis)

	for(i in 1:r){
		splineRes <- smooth.spline(x = 1:t,y = basis[,i], df = t*0.2)
		basis[,i] <- predict(splineRes, 1:t)$y
	}
	# basis <- (svdRes$u %*% (diag(svdRes$d)))[,1:r]
	# coef <- t(svdRes$v)[1:r,]
	print(dim(basis))
	basis <- cbind(basis,1)
	coef <- lsfit(x = basis, y= DOdata, intercept = FALSE)$coefficients
	# print(dim(svdRes))
	
	DO_fit <- basis %*% coef # coef, each column is the coefficents for different basis
	
	return(list(fit = DO_fit, coef = coef, basis = basis,varExpl = sum(svdRes$d[1:r])/sum(svdRes$d)))

}


B_spline <- function(DOdata,knots,norder = 4){
	
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



crossValidation <- function(DOdata,locations){
	time <- index(DOdata)
	for(i in 1:length(locations$loggerID)){
		remainingloc <- locations$loggerID[-i] %>% as.character()
		subData <- DOdata[,remainingloc]
		grid <- locations[i,c("longitude","latitude","x","y")]
		grid$convexIndex <- 1
		
		Y <- as.numeric(DOdata[,as.character(locations$loggerID[i])])

		trendPrediction <- basis_interpolation(subData,locations[-i,],grid,basis = "svd", r = length(remainingloc))
		finalPrediction <- cbind(t(trendPrediction$trend),grid)[,1:length(time)] %>% as.numeric()
		print(plot(time,Y,type = "l",main = locations$loggerID[i]))
		print(lines(time,finalPrediction,col="red"))
	}
}


#finalPrediction <- cbind(t(trendPrediction$trend),grid)



# B_spline <- function(DOdata,logger_geo,timeIndex, norder = 4){
# 	# DOdata is a dataframe, each column is a series of data from one sensor
# 	DOdata <- DOdata[1:104,] %>% as.matrix()
	
# 	norder <- 4
# 	knots <- seq(0,110,7)
# 	nBasis <- length(knots)-2+norder
# 	bbasis <- create.bspline.basis(c(0,105),nBasis,norder,knots)

# 	basismat <- eval.basis(1:104, bbasis)
# 	coef <- solve(crossprod(basismat), crossprod(basismat,DOdata))
	
# 	DOfd <- fd(coef, bbasis,list("Day", "Year", "DO (mg/L)"))
# 	# plot(DOfd, lty=1, lwd=2, col=1)

# 	DO_fit <- predict(DOfd,1:104)
# 	DO_res <- (DOdata-DO_fit) %>% as.data.frame()


# 	for(i in 1:ncol(DOdata)){
# 	#	plotfit.fd(DOdata[,i], 1:104, DOfd[i], lty=1, lwd=2)
# 	}
	 
# 	coef_t <- t(coef) %>% as.data.frame()
# 	coef_t$ID <- as.numeric(rownames(coef_t))

# 	coef_df <- merge(coef_t,logger_geo,by.x = "ID",by.y = "loggerID")

	
# 	return(list(fit = DO_fit,res = DO_res,coef = coef_df,basis = bbasis))
# }




# test
# debug(fda)





# fda_res <- fda(data_2014,bottomLogger_2014,time_2014)


# # spatial analysis on the coefficients
# coef_df <- fda_res[[3]]
# coef_df2 <- coef_df
# coordinates(coef_df2) = ~x+y
# bbasis <- fda_res[[4]]

# for(i in 1:18){
# 	formu <- as.formula(sprintf("bspl4.%d~1",i))
# 	print(plot(variogram(formu,data =coef_df2,cutoff = 120)))
# }

# loci <- createGrid(bottomLogger_2014) %>% lonlat2UTM()


# krigingCoeff <- matrix(0,nrow(loci),18)

# for(i in 1:18){
# 	subDF <- coef_df[,c("x","y",paste("bspl4.",i,sep = ""),"ID")] %>% as.geodata()
# 	ml <- likfit(subDF, ini = c(1,0.5), fix.nugget = T, method = "RML") # using loglikelihood to estimate
# 	krigingCoeff[,i] <- krige.conv(subDF, locations = loci[,c("x","y")], krige = krige.control(obj.m = ml))$predict
	
# 	# reconstruct
# 	# plot 
# 	# bin1 <- variog(subDF)
# 	# plot(bin1, main = expression(paste("fixed ")))
# 	# lines(ml)
# }

# # construct 
# fdobj <- fd(t(krigingCoeff),bbasis)
# allInterpolation <- predict(fdobj,1:104)  # a matrix with 104* nGrid

# allDf <- list()
# for(i in 1:104){
# 	loci$pred <- allInterpolation[i,]
# 	attributes(loci)$time <- i
# 	allDf[[i]] <- subset(loci,convexIndex ==1)
# }

# plot_spatial(allDf,bottomLogger_2014,"./output/newRes")


# # bayesian kriging
# MC <- model.control(cov.model = "exponential")  # model control
# PC <- prior.control(phi.discrete=seq(30,100,5)) # priors

# OC <- output.control(n.pos = 500, n.pred = 100, quantile = c(0.25, 0.5, 0.75), thres = 1.5)
# # n.pos: number of samples to be taken from the posterior distribution
# # n.pred: number of samples to be taken from the predictive distribution
# # 

# s100.kb <- krige.bayes(subDF, loc=loci[1:10,c("x","y")], model=MC, prior=PC, out=OC)

# # # spatial temporal analysis
# # DO_res <- fda_res$res
# # DO_res$logger_time <- as.POSIXct(time_2014)
# # DO_res <- melt(DO_res,id.vars = c("logger_time")) %>% 
# # 	arrange(logger_time,variable) %>%
# # 	rename(ID = variable)

# # coordinates(bottomLogger_2014)=~longitude+latitude
# # raster::projection(bottomLogger_2014)=CRS("+init=epsg:4326")

# # # year <- 2015
# # # bottomLogger_2015 <- retriveGeoData(year,"B")
# # # data_2015 <- retriveLoggerData(bottomLogger_2015$loggerID,year,"DO","daily","AVG",transform=TRUE)

# # timeDF <- STFDF(sp=bottomLogger_2014,time=time_2014,data=data.frame(value = DO_res$value))
# vST <- variogramST(value~longitude+latitude+bathymetry+I(bathymetry^2),timeDF,tlags = 0:10,cutoff = 120)
# plot(vST,map=FALSE)

# fda(DOdata,logger)
