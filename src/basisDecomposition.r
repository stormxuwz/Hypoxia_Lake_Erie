# This is the script to represent data with basis functions
# Available basis: B-Spline and temporal basis function through SVD

source("./src/database.R")
source("config.R")
source("./src/plot.R")
source("./src/spatialHelper.R")
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

NMF_basis <- function(DOdata, r){
	require(NMF)
	DOdata <- as.matrix(DOdata)
	nmfRes <- nmf(DOdata, r)

	W <- nmfRes@fit@W
	H <- nmfRes@fit@H
	return(list(basis = W, coef = H))
}


SVD_basis <- function(DOdata, r){
	# r is the column vectors to keep
	DOdata <- as.matrix(DOdata)
	svdRes <- svd(DOdata) # is centered necessary? No (3/23)

	basis <- svdRes$u[,1:r]
	t = nrow(basis)

	for(i in 1:r){
		splineRes <- smooth.spline(x = 1:t,y = basis[,i], df = t*0.5)
		basis[,i] <- predict(splineRes, 1:t)$y
	}
	basis <- cbind(basis,1)
	coef <- lsfit(x = basis, y= DOdata, intercept = FALSE)$coefficients
	DO_fit <- basis %*% coef # coef, each column is the coefficents for different basis
	
	return(list(fit = DO_fit, coef = coef, basis = basis,varExpl = sum(svdRes$d[1:r])/sum(svdRes$d)))
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
