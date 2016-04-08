source("./src/database.R")
source("config.R")
library(leaflet)
year <- 2014
var <- "DO"
groupRange <- "hourly"
dataType <- "AVG"
timeRange <- c("2014-07-01","2014-09-31")

allGeoData <- retriveGeoData(2014,"B")
allData <- retriveLoggerData(allGeoData$loggerID,year,var,groupRange,dataType,timeRange=timeRange)

t <- allData[,"10523447"]
dt <- diff(c(t))
m <- tseries::arma(dt,c(2,1))

t <- ts(c(t),start=1,freq=24)
decomposeddata<-stl(t, s.window=24)