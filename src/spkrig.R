library(gstat)
library(sp)
library(raster)
library(reshape2)

detrending <- function(data,sp,time){
	library(fields)
	if(ncol(data)!=nrow(sp)){
		print("Wrong")
		return(0)
	}

	myDf <- data.frame(time=rep(time,each=nrow(sp)),
						long = rep(sp$longitude,nrow(data)),
						lat = rep(sp$latitude,nrow(data)),
						bathy= rep(sp$bathymetry,nrow(data)),
						DO = c(t(as.matrix(data))))

	detrendingModel <- lm(DO~long+lat+bathy+I(bathy^2)+time,data=myDf)
	# detrendingModel <- Tps(DO~Long)

	return(list(res=detrendingModel$residuals,detrendingModel=detrendingModel))
}


conditionalSim <- function(data, spMeta, grid){
	data <- DO_bottom_dailyMean
	spMeta <- meta_B[,c("Lat","Long","bathy")]
	spMeta$bathy <- ft2meter(spMeta$bathy)
	x <- as.matrix(spMeta)
	names(grid) <- c("Long","Lat","bathy")

	coordinates(spMeta)=~Long+Lat
	coordinates(grid)=~Long+Lat

	
	krigData <- spMeta

	timeLength <- nrow(data)

	require(gridExtra)
	for(i in 1:nrow(data)){
		print(i)

		
		spMeta$DO <- as.numeric(data[i,1:18])
		tps_model<-Tps(x,as.numeric(spMeta$DO),df=10)
		tps_model_pred=predict(tps_model,as.data.frame(grid)[,c("Lat","Long","bathy")])
        spMeta$res<-tps_model$residuals

		v <- variogram(res~bathy+Long+Lat, spMeta)
		v <- subset(v,np>3)
		# m <- fit.variogram(v)
		png(paste("./variogram/day",i,".png",sep="_"),width=800,height=800)

		p1 <- plot(variogram(res~bathy+Long+Lat, spMeta,cloud=T),plot.numbers = TRUE)
		p2 <- plot(variogram(DO~bathy+Long+Lat, spMeta,cloud=T),plot.numbers = TRUE)
		p3 <- plot(variogram(res~bathy+Long+Lat, spMeta,cressie=T),plot.numbers = TRUE)
		p4 <- plot(variogram(DO~bathy+Long+Lat, spMeta,cressie=T),plot.numbers = TRUE)
		p5 <- plot(variogram(res~bathy+Long+Lat, spMeta,cressie=F),plot.numbers = TRUE)
		p6 <- plot(variogram(DO~bathy+Long+Lat, spMeta,cressie=F),plot.numbers = TRUE)

		grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
		dev.off()
		# sim <- krige(DO~Long+Lat+bathy, spMeta, grid, model = m, nmax = 15, beta = 5.9, nsim = 100)
		# krige(DO,sp_data,grid)$var1.pred

	}


}


spatialTemporalKrig <- function(data,locationData,stationName,endTime=as.POSIXct("2014-10-1"),subsample=TRUE){
	# data is the 
	require(geoR)
	require(raster)
	require(spacetime)
	require(gstat)

	
	data <- subset(data,Time<endTime)
	myTime <- data$Time
	myData <- data[,c(stationName)]

	if(subsample){
		subseq <- seq(1,nrow(myData),by=6)
		myData <- myData[subseq,]
		myTime <- myTime[subseq]
	}

	myLocation <- locationData

	# Create the data
	DOSeries <- c(t(as.matrix(myData)))
	detrendingModel <- detrending(myData,myLocation,myTime)

	res <- detrendingModel$res
	model <- detrendingModel$detrendingModel

	myLocation <- myLocation[,c("longitude","latitude")]
	coordinates(myLocation)=~longitude+latitude
	projection(myLocation)=CRS("+init=epsg:4326")

	myDataFrame <- data.frame(DO=DOSeries,res=res)
	timeDF <- STFDF(myLocation,myTime,data=myDataFrame)

	return(list(timeDF=timeDF,res=res,detrendingModel=model))
}


# mytimeDF <- spatialTemporalKrig(DO_bottom2,meta_B,stationName,endTime=as.POSIXct("2014-10-31"),subsample=TRUE)
# res_vST<- variogramST(res~1,data=mytimeDF,tunit="days",assumeRegular=T,na.omit=T,tlags=0:15,progress=T)
# plot(res_vST,map=F)
# res_vST2<- variogramST(res~1,data=mytimeDF)

# mytimeDF_sp <- as.data.frame(mytimeDF@sp)

# for(i in c(1:2)){
# 	for(j in c(1:nrow(mytimeDF_sp))){
# 		loggerName <- rownames(mytimeDF_sp)[j]
# 			pdf(paste("./plot/",loggerName,ifelse(i==1,"_DO2","_DO_res2"),".pdf",sep=""))
# 			plot(mytimeDF[j,,i],main=loggerName)
# 			dev.off()
# 	}
# }

library(RMySQL)
conn <- dbConnect(MySQL(), dbname = "DO2014", username="root", password="XuWenzhaO", host="127.0.0.1", port=3306)
meta <- dbReadTable(conn,"loggerInfo")
meta_B <- subset(meta,available==1 & loggerPosition=="B")

data  <- list()

for(i in 1:nrow(meta_B)){
	loggerID <- meta_B[i,"loggerID"]
	sql <- paste("select logger, date(Time) as Time, AVG(DO) as DO_avg from loggerData where logger=",loggerID," AND Time < '2014-09-25' AND Time > '2014-07-10' Group by date(Time)")
	res <- dbGetQuery(conn,sql)
	data[[paste("logger_",loggerID,sep="")]] <- res$DO_avg
}

data$Time <- as.POSIXct(res$Time)
DO_bottom <- as.data.frame(data)
stationName=paste("logger",meta_B$loggerID,sep="_")


melt_DO <- melt(DO_bottom,id = "Time")
meta_B$newLogger <- paste("logger_",meta_B$loggerID,sep="")
melt_DO <- base::merge(melt_DO,meta_B,by.x="variable",by.y="newLogger",all.x=TRUE)


# do detrending
tps_X = matrix(0,nrow(melt_DO),3)
tps_X[,1]=as.numeric(melt_DO$Time)
tps_X[,2]=as.numeric(melt_DO$latitude)
tps_X[,3]=as.numeric(melt_DO$longitude)

Tps_model <- Tps(tps_X,melt_DO$value,df=10)
melt_DO$res <- Tps_model$residuals
melt_DO2 <- melt_DO

DO_bottom2 <- dcast(melt_DO2[,c("Time","variable","res")],Time~variable,value.var="res")




allVariogram <- function(data,meta_B){
	spdata <- meta_B
	coordinates(spdata) =~ longitude+latitude
	projection(spdata)=CRS("+init=epsg:4326")
	v <- data.frame()
	for(i in 1:nrow(data)){ 
		spdata$val <- as.numeric(data[i,])
		v <- rbind(v,data.frame(variogram(val~1,data=spdata,cloud=T,cutoff=10000)))
	}
	return(v)
}


allv <- allVariogram(subset(DO_bottom2,select = -Time),meta_B)

m <- ggplot(allv, aes(y = gamma, x = dist,group = round_any(dist, 5)))
m+geom_boxplot()+xlim(c(0,140/2))





mytimeDFList <- spatialTemporalKrig(DO_bottom,meta_B,stationName,endTime=as.POSIXct("2014-10-31"),subsample=FALSE)

detrendingModel <- mytimeDFList$detrendingModel
mytimeDF <- mytimeDFList$timeDF




res_vST<- variogramST(DO~longitude+latitude,data=mytimeDF,assumeRegular=T,na.omit=T,tlags=0:15,progress=F)
res_vST<- variogramST(res~1,data=mytimeDF,tunit="days",assumeRegular=T,na.omit=T,tlags=0:15,progress=T)


# Fit the variogram
pars.l <- c(sill.s = 0, range.s = 10, nugget.s = 0,sill.t = 0, range.t = 1, nugget.t = 0,sill.st = 0, range.st = 10, nugget.st = 0, anis = 0)
pars.u <- c(sill.s = 200, range.s = 1000, nugget.s = 100,sill.t = 200, range.t = 60, nugget.t = 100,sill.st = 200, range.st = 1000, nugget.st = 100,anis = 700)

separable <- vgmST("separable", space = vgm(-60,"Sph", 500, 1),time = vgm(35,"Sph", 500, 1), sill=0.56)
sumMetric <- vgmST("sumMetric", space = vgm(psill=5,"Sph", range=500, nugget=0),time = vgm(psill=500,"Sph", range=500, nugget=0), joint = vgm(psill=1,"Sph", range=500, nugget=10), stAni=500)

separable_Vgm <- fit.StVariogram(res_vST, separable, fit.method=10)
# sumMetric_Vgm <- fit.StVariogram(res_vST, sumMetric, method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="days")
sumMetric_Vgm <- fit.StVariogram(res_vST, sumMetric, method="L-BFGS-B",lower=pars.l,upper=pars.u)
plot(res_vST, list(sumMetric_Vgm),map=FALSE)



# outlier detection

# stKrig_10384436

for(i in 1:nrow(meta_B)){
	print(i)
	loggerID <- meta_B[i,"loggerID"]
	loggerInfo <- meta_B[i,c("longitude","latitude","bathymetry")]

	sql <- paste("select logger, Time, DO from loggerData where logger=",loggerID," AND Time < '2014-10-09' AND Time > '2014-06-25'")
	series <- dbGetQuery(conn,sql)

	timeSeries <- as.POSIXct(series$Time)
	# length(timeSeries)

	myGrid <- data.frame(long=rep(loggerInfo$longitude,length(timeSeries)),lat=rep(loggerInfo$latitude,length(timeSeries)),bathy=rep(loggerInfo$bathymetry),time=timeSeries)
	trend <- predict(detrendingModel$detrendingModel,myGrid)
	cat("trending finished \n")

	spGrid <- data.frame(longitude=loggerInfo$longitude,latitude=loggerInfo$latitude)

	coordinates(spGrid)=~longitude+latitude
	projection(spGrid)=CRS("+init=epsg:4326")
	grid.ST <- STF(spGrid,timeSeries)
	pred <- krigeST(res~1, data=mytimeDF, modelList=sumMetric_Vgm, newdata=grid.ST,computeVar = T)
	pred_res <- pred@data$var1.pred
	pred_var <- pred@data$var1.var
	
	cat("res interpolation finished \n")
	finalPred <- data.frame(time=timeSeries,pred=pred_res+trend,var=pred_var,obs=series$DO)

	finalPred2 <- zoo(finalPred[,-1],order.by=timeSeries)

	finalPred2$upperBoundary <- finalPred2
	p1 <- qplot(time,pred,data=finalPred)+geom_point(aes(time,obs),data=finalPred,color="red")+geom_line(aes(time,pred+pred_var^0.5),data=finalPred)+geom_line(aes(time,pred-pred_var^0.5),data=finalPred)+ggtitle(loggerID)+ylim(c(0,13))
	p2 <- qplot(time,pred,data=finalPred)+geom_point(aes(time,obs),data=finalPred,color="red")+geom_point(aes(time,pred+pred_var^0.5),data=finalPred,color="blue")+geom_point(aes(time,pred-pred_var^0.5),data=finalPred,color="green")+ggtitle(loggerID)+ylim(c(0,13))
	png(paste("stKrig_",loggerID,".png",sep=""),width=1000,height=800)
	print(p1)
	dev.off()
}


# spatial kriging

