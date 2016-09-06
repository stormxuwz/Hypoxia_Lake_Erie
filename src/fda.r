## script to run functional data anlayis
source("./src/database.R")
source("config.R")
source("./src/plot.R")
source("./src/spatialHelper.R")

require(fda)
require(sp)
require(gstat)
require(dplyr)
require(reshape2)

fda <- function(DOdata,logger_geo,time){
	# DOdata is a dataframe, each column is a series of data from one sensor
	DOdata <- DOdata[1:104,] %>% as.matrix()
	
	norder <- 4
	knots <- seq(0,110,7)
	nbasis <- length(knots)-2+norder
	bbasis <- create.bspline.basis(c(0,105),nbasis,norder,knots)

	basismat <- eval.basis(1:104, bbasis)
	coef <- solve(crossprod(basismat), crossprod(basismat,DOdata))
	
	DOfd <- fd(coef, bbasis,list("Day", "Year", "DO (mg/L)")) 
	plot(DOfd, lty=1, lwd=2, col=1)

	DO_fit <- predict(DOfd,1:104)
	DO_res <- (DOdata-DO_fit) %>% as.data.frame()


	for(i in 1:ncol(DOdata)){
	#	plotfit.fd(DOdata[,i], 1:104, DOfd[i], lty=1, lwd=2)
	}
	 
	coef_t <- t(coef) %>% as.data.frame()
	coef_t$ID <- as.numeric(rownames(coef_t))

	coef_df <- merge(coef_t,logger_geo,by.x = "ID",by.y = "loggerID")

	coordinates(coef_df) = ~x+y

	for(i in 1:18){
		formu <- as.formula(sprintf("bspl4.%d~x+y+bathymetry",i))
		print(plot(variogram(formu,data =coef_df, width = 1)))
	}
	return(list(fit = DO_fit,res = DO_res,coef = coef))
}




# test
# debug(fda)
timeRange <- c("2014-06-23","2014-10-05")


year <- 2014
bottomLogger_2014 <- retriveGeoData(year,"B") %>% arrange(loggerID) %>% lonlat2UTM()
data_2014 <- retriveLoggerData(bottomLogger_2014$loggerID,year,"DO","daily","AVG",timeRange = timeRange)
time_2014 <- index(data_2014)
fda_res <- fda(data_2014,bottomLogger_2014,time_2014)

DO_res <- fda_res$res
DO_res$logger_time <- as.POSIXct(time_2014)
DO_res <- melt(DO_res,id.vars = c("logger_time")) %>% 
	arrange(logger_time,variable) %>%
	rename(ID = variable)

coordinates(bottomLogger_2014)=~longitude+latitude
raster::projection(bottomLogger_2014)=CRS("+init=epsg:4326")

# year <- 2015
# bottomLogger_2015 <- retriveGeoData(year,"B")
# data_2015 <- retriveLoggerData(bottomLogger_2015$loggerID,year,"DO","daily","AVG",transform=TRUE)

timeDF <- STFDF(sp=bottomLogger_2014,time=time_2014,data=data.frame(value = DO_res$value))
vST <- variogramST(value~longitude+latitude+bathymetry+I(bathymetry^2),timeDF,tlags = 0:10)
plot(vST,map=FALSE)

# fda(DOdata,logger)
