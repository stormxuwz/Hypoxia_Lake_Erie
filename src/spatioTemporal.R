library(sp)
library(spacetime)
library(gstat)
source("./src/database.R")
source("config.R")
library(dplyr)
# estimate the thermocline depth

spKriging <- function(loggerIndex,year,var,groupRange,dataType,timeRange=NULL){
	logger_geo <- retriveGeoData(2014,"B",loggerIndex)
	logger_geo <- arrange(logger_geo,loggerID)[,c("longitude","latitude")] # the sp points shouldn't have other features
	coordinates(logger_geo)=~longitude+latitude
	projection(logger_geo)=CRS("+init=epsg:4326")
	
	logger_value <- retriveLoggerData(loggerIndex,2014,"DO","hourly","AVG",c("2014-07-01","2014-07-31"),transform = F)
	logger_value$Time <- as.POSIXct(logger_value$Time)
	logger_value <- arrange(logger_value,Time,logger)
	
	logger_time <- unique(logger_value$Time)
	
	timeDF <- STFDF(logger_geo,	logger_time,logger_value)
	vST<- variogramST(DO~longitude+latitude,data=timeDF,assumeRegular=T,na.omit=T,tlags=1:5,progress=F)
	return(vST)
}



loggerIndex <- c(10523450,10523447,10523437,10523445)
year <- 2014
var <- "DO"
groupRange <- "hourly"
dataType <- "AVG"
timeRange <- c("2014-07-01","2014-08-31")






head(testData)


testData <- arrange(testData,Time)

?STFDF
