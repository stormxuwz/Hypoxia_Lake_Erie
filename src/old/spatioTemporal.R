library(sp)
library(spacetime)
library(gstat)
source("./src/database.R")
source("config.R")
library(dplyr)
library(dygraphs)
library(reshape2)

# estimate the thermocline depth

stKriging <- function(loggerIndex,year,var,groupRange,dataType,timeRange,resultName){
	logger_geo <- retriveGeoData(2014,"B",loggerIndex)
	logger_geo <- arrange(logger_geo,loggerID) # the sp points shouldn't have other features
	coordinates(logger_geo)=~longitude+latitude
	raster::projection(logger_geo)=CRS("+init=epsg:4326")
	
	logger_value <- retriveLoggerData(loggerIndex,2014,var,groupRange,dataType,timeRange,transform = TRUE)

	coeff1=c()
	coeff2=c()
	logger_res <- logger_value
	
	for(i in 1:ncol(logger_value)){
		detrendingResults <- detrending(logger_value[,i])
		coeff1 <- c(coeff1,detrendingResults$coeff[1])
		coeff2 <- c(coeff2,detrendingResults$coeff[2])
	#	logger_res[,i] <- detrendingResults$res
	}
	
	logger_time <- as.POSIXct(time(logger_value))

	logger_res <- as.data.frame(logger_res)
	
	rownames(logger_res) <- NULL
	logger_res$dayTime <- as.character(logger_time)
	logger_res_original <- logger_res
	logger_res <- melt(logger_res,id.vars = c("dayTime"))
	logger_res$dayTime <- as.POSIXct(logger_res$dayTime)
	logger_res <- arrange(logger_res,dayTime,variable)
	yday <- as.numeric(strftime(logger_time, format = "%j"))
	# logger_value <- retriveLoggerData(loggerIndex,2014,var,groupRange,dataType,timeRange,transform = FALSE)
	# logger_value$Time <- as.POSIXct(logger_value$Time)
	# logger_value <- arrange(logger_value,Time,logger)
	# logger_time <- unique(logger_value$Time)
	
	timeDF <- STFDF(sp=logger_geo,time=logger_time,data=data.frame(value = logger_res$value))
	#vST<- variogramST(value~longitude+latitude+time,timeDF,tlags=0:12,cutoff = 100,width = 20,cressie=T)
	vST <- variogramST(value~longitude+latitude+time+bathymetry+I(bathymetry^2),timeDF,tlags=0:4,boundaries=c(0,25,50,75))
	# vST <- variogramST(value~1,timeDF,tlags=0:4,boundaries=seq(0,70,by=20))
	
	# vST2 <- subset(vST,dist>20 | dist<10)

	sumMetric <- vgmST("sumMetric", space = vgm(psill=0.5,"Sph", range=30, nugget=0),
										 time = vgm(psill=0.5,"Sph", range=10, nugget=0), 
										 joint = vgm(psill=0.2,"Sph", range=20, nugget=0),stAni = 10, temporalUnits="days")
	# plot(vST,map=F)
	plot(vST,list(sumMetric),map=FALSE)

	# png(paste("~/Developer/Hypoxia/output/",resultName,".png",sep=""))
	# print(plot(vST,map=F))
	# dev.off()
	return(vST)
}


fftDetrending <- function(i){
	t <- as.numeric(strftime(logger_time, format = "%j"))
	print(names(logger_value)[i])
	plot(logger_value[,i])
	res <- lm(as.vector(logger_value[,i])~t+I(t^2))
	fitValue <- res$fitted.values
	fitValue <- zoo(fitValue,order.by = logger_time)
	lines(fitValue,col="red")
	plot(res$residuals,type = "b")
	fftModel <- fft(res$residuals)
	plot(abs(as.vector(fftModel))[1:as.integer(length(fftModel)/2)],type ="b",main = names(logger_value)[i],ylab = "amplitude")
}





reSummaryVariogram <- function(vST){
	
	distUniq <- unique(vST$dist)
	distT <- unique(vST$timelag)
	summaryList <- list()
	newVST <- data.frame()
	for(dist_ in sort(distUniq)){
		for(timelag_ in sort(distT)){
			subVST <- subset(vST,dist == dist_ & timelag == timelag_)
			# if(dist_ == 0)
				# newGamma <- median(subVST$gamma)/9
			# else if(timelag_ >0)
			summaryList[[paste(dist_,timelag_,sep="__")]] <- summary(subVST$gamma)
			# newGamma <- median(subVST$gamma)
			newGamma <- quantile(subVST$gamma,0.25)

			newVST <- rbind(newVST,data.frame(np=1,dist = dist_,gamma = newGamma,dir.hor = 0, dir.ver = 0,id = paste("lag",timelag_,sep=""),timelag = timelag_,spacelag = 0,avgDist = dist_))
		}
	}
	newVST <- arrange(newVST,timelag,dist)
	class(newVST) <- class(vST)
}


detrendingAll <- function(logger_value,logger_geo,logger_time){
		logger_value <- as.data.frame(logger_value)
		rownames(logger_value) <- NULL
		logger_value$time <- logger_time
		logger_value <- melt(logger_value,id.vars = c("time"))
		logger_value$time <- as.POSIXct(logger_value$time)
		logger_value <- arrange(logger_value,time,variable)
		logger_value$yday <- as.numeric(strftime(logger_time, format = "%j"))
		logger_value <- base::merge(logger_value,as.data.frame(logger_geo),by.x = "variable",by.y = "loggerID",all.x = TRUE)
		
		tpsModel <- fields::Tps(x = as.matrix(logger_value[,c("longitude","latitude","bathymetry")]),Y = as.vector(logger_value$value))
}


detrending <- function(series,base = "linear"){
	data <- data.frame(value=series,time= 1:length(series))

	model <- lm(value~time+I(time^2),data=data)  # the coefficient is the value/second
	trendCoeff <- c(model$coefficients)
	res <- model$residuals

	return(list(coeff = trendCoeff,res = res))
}




totalLoggerIndex <- c(10523446,10523447,10523439,10523450,10523436,10528846,10523443,10523437,10523445)

# loggerIndex <- c(10528846,10523445,10523439)  # 6,9,3
loggerIndex <- totalLoggerIndex
year <- 2014
var <- "DO"
groupRange <- "daily"
dataType <- "AVG"
timeRange <- c("2014-06-15","2014-10-05")
# 
# for(i in 1:length(totalLoggerIndex)){
# 	loggerIndex <- totalLoggerIndex[-i]
# 	resultName <- paste("-",totalLoggerIndex[i],sep="")
# 	vST <- spKriging(loggerIndex,year,var,groupRange,dataType,timeRange,resultName)	
# }
# 
