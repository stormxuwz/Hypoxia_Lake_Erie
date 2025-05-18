# This is the script to represent data with basis functions
# Available basis: B-Spline and temporal basis function through SVD
source("./src/database.R")
source("./src/plot.R")
source("./src/helper.R")

require(sp)
require(gstat)
require(dplyr)
require(reshape2)
require(geoR)

coeff2Value <- function(coef, basis) {
  return(basis %*% t(coef)) # return a matrix
}

SVD_basis <- function(DOdata, r) {
  # r is the column vectors to keep
  DOdata <- as.matrix(DOdata)

  svdRes <- svd(DOdata %>% scale()) # Scale matrix before SVD

  basis <- svdRes$u[, 1:r]
  t <- nrow(basis)

  for (i in 1:r) {
    splineRes <- smooth.spline(x = 1:t, y = basis[, i], df = t * 0.2)
    basis[, i] <- predict(splineRes, 1:t)$y
  }

  basis <- cbind(basis, 1) # add the bias term

  coef <- lsfit(x = basis, y = DOdata, intercept = FALSE)$coefficients
  DO_fit <- basis %*% coef # coef, each column is the coefficents for different basis

  print("variance explained by SVD")
  print(sum(svdRes$d[1:r]**2) / sum(svdRes$d**2))

  return(list(fit = DO_fit, coef = coef, basis = basis, varExpl = sum(svdRes$d[1:r]**2) / sum(svdRes$d**2)))
}
