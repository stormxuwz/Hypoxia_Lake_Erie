require(rgdal)
require(raster)

findBathy <- function(spData){
	
	coordinates(spData)=~Long+Lat
	# spData must be a sp data class
	bathymetry_raster = raster("~/Downloads/erie_lld/erie_lld.asc")
	bathymetry <- extract(bathymetry_raster,spData)
	return(bathymetry)
}


findIndex <- function(spData){
	require(ncdf)
	# library(abline)
	ncFile = open.ncdf(paste("../Project_IO/DO_LakeErie/OtherData/","e201415100.out1.nc",sep=""))
	lat <- c(get.var.ncdf(ncFile,"lat"))
    lon <- c(get.var.ncdf(ncFile,"lon"))

    grid <- data.frame(Long=lon,Lat=lat)
    coordinates(spData)=~Long+Lat
    coordinates(grid)=~Long+Lat
    projection(spData)=CRS("+init=epsg:4326")
    projection(grid)=CRS("+init=epsg:4326")
    spMatrix <- spDists(spData,grid)
    return(apply(spMatrix,1,which.min))
}


getInfoNCFile <- function(index,var="2Dwave"){
    require(ncdf)
    library(abind)

    ny <- 87
    nx <- 193

    ind_y <- ceiling(index/nx)
    ind_x <- index%%nx
    ind_x <- ifelse(ind_x==0,nx,ind_x)

    startTime <- strptime("2014-6-1 00:00","%Y-%m-%d %H:%M")

    if(var=="2Dwave"){
        fileNameList <- c("e201415100.out1.nc","e201420100.out1.nc","e201425100.out1.nc")
        timeSeries <- seq(from=startTime,by=3600*1,length.out=1200*3)
        myuc <- rep(-1,length.out=1200*3)
        myvc <- myuc 
        mywvd <- myuc 
        
        for(i in 1:3){
        	seriesRange <- ((i-1)*1200+1):((i)*1200)

        	ncFile = open.ncdf(paste("../Project_IO/DO_LakeErie/OtherData/",fileNameList[i],sep=""))
	        uc <- get.var.ncdf(ncFile,"uc")
	        vc <- get.var.ncdf(ncFile,"vc")
	        wvd <- get.var.ncdf(ncFile,"wvd")
	 		
	 		lat <- get.var.ncdf(ncFile,"lat")
        	lon <- get.var.ncdf(ncFile,"lon")	

        	cat(lat[ind_x,ind_y],lon[ind_x,ind_y],"\n")

	 		myuc[seriesRange] <- uc[ind_x,ind_y,]
	 		myvc[seriesRange] <- vc[ind_x,ind_y,]
	 		mywvd[seriesRange] <- wvd[ind_x,ind_y,]

        	close.ncdf(ncFile)
        }
        return(data.frame(time=timeSeries,uc=myuc,vc=myvc,wvd=mywvd))
    }
    else if(var=="Temp"){
    	fileNameList <- c("e201415100.out3.nc","e201420100.out3.nc","e201425100.out3.nc")
    	timeSeries <- seq(from=startTime,by=3600*3,length.out=400*3)
    	myTemp <- rep(-1,length.out=400*3*20)
    	myTemp <- matrix(myTemp,400*3,20)
      	# float d3d[nx,ny,nsigma]  Longname:3D Depth at Nodes Missval:1e+30"

        for(i in 1:3){
        	seriesRange <- ((i-1)*400+1):((i)*400)

        	ncFile = open.ncdf(paste("../Project_IO/DO_LakeErie/OtherData/",fileNameList[i],sep=""))
        	temp <- get.var.ncdf(ncFile,"temp")
        	d3d <- get.var.ncdf(ncFile,"d3d")
        	bathy <- get.var.ncdf(ncFile,"depth")
        	# float temp[nx,ny,nsigma,time]  Longname:Temperature Missval:-99999

        	myTemp[seriesRange,] <- t(temp[ind_x,ind_y,,])
        	myd3d <- d3d[ind_x,ind_y,]
        	mybathy <- bathy[ind_x,ind_y]
        }
        myTemp <- as.data.frame(myTemp)
        # names(myTemp) <- as.character(myd3d)
        myTemp$time <- timeSeries
        return(list(myTemp=myTemp,d3d=myd3d,bathy= mybathy))
    }
}
        
scale2Unit <- function(s){
	return((s-min(s,na.rm=T))/(max(s,na.rm=T)-min(s,na.rm=T)))
}     

find2dWaveData <- function(meta){
	logger_index_inNCFile <- findIndex(meta)

	waveData <- list()
	for(i in 1:length(logger_index_inNCFile) ){
		print(i)
		loggerName <- paste("logger_",meta$No[i],sep="")
		waveData[[loggerName]] <- getInfoNCFile(logger_index_inNCFile[i],var="2Dwave")
	}
	saveRDS(waveData,"2DwaveData.rds")
}

find3dTempData <- function(meta){
	logger_index_inNCFile <- findIndex(meta)

	tempData <- list()
	for(i in 1:length(logger_index_inNCFile) ){
		print(i)
		loggerName <- paste("logger_",meta$No[i],sep="")
		tempData[[loggerName]] <- getInfoNCFile(logger_index_inNCFile[i],var="Temp")
	}

	saveRDS(tempData,"3DtempData.rds")
	return(tempData)
}

plot_wave_DO <- function(meta_available){
	for(i in 1:nrow(meta_available)) {
		loggerName <- paste("logger_",meta_available$No[i],sep="")
		DO_data <- DO_raw_list[[as.character(meta_available$No[i])]]
		wave_data <- waveData[[loggerName]]

		p1=ggplot()+geom_point(aes(Time,scale2Unit(DO)),data=DO_data)+geom_line(aes(time,scale2Unit(uc^2+vc^2)),data=wave_data,color="red")
		p2=ggplot()+geom_point(aes(Time,scale2Unit(DO)),data=DO_data)+geom_line(aes(time,scale2Unit(uc)),data=wave_data,color="yellow")
		p3=ggplot()+geom_point(aes(Time,scale2Unit(DO)),data=DO_data)+geom_line(aes(time,scale2Unit(vc)),data=wave_data,color="green")
		p4=ggplot()+geom_point(aes(Time,scale2Unit(DO)),data=DO_data)+geom_line(aes(time,scale2Unit(wvd)),data=wave_data,color="blue")

		png(paste(loggerName,"wave.png"),width=800,height=1200)
		grid.arrange(p1,p2,p3,p4,nrow=4)
		dev.off()
	}
}

ft2meter <- function(ft){
	return(ft/3.2808399)
}


plotMultipleDO <- function(startTime,endTime,loggerList=c(),colorList=c()){
	p <- ggplot()
	if(length(loggerList)!=length(colorList)){
		cat("!!")
		return()
	}
	startTime <- as.POSIXct(startTime,tz="EST")
	endTime <- as.POSIXct(endTime,tz="EST")

	for(i in 1:length(loggerList)){
		logger <- loggerList[i]
		col <- colorList[i]
		data <- subset(DO_raw_list[[as.character(logger)]],Time>startTime & Time < endTime)
		p <- p+geom_line(aes(Time,DO),data=data,color=col)+geom_point(aes(Time,DO),data=data,color=col)
	}
	return(p)
}


sendtoSQL_loggerData <- function(DO_raw_list){
	alldf <- data.frame(Time=c(),DO=c(),Temp=c(),logger=c())

	for(i in 1:length(DO_raw_list)){
		data <- DO_raw_list[[i]]
		data <- na.omit(data)
		n <- nrow(data)
		data <- data[10:(n-10),]
		# data$Time <- as.character(data$Time)
		data$logger  <- as.numeric(names(DO_raw_list)[i])

		if(max(data$Temp,na.rm = T)>50){
			data$Temp <- (data$Temp - 32) *5/9
		}
		# data <- melt(data[,-1],id=c("Time","logger"))
		alldf<- rbind(alldf,data[,-1])
	}
	alldf$id=1:nrow(alldf)
	dbWriteTable(conn, "loggerData", alldf, overwrite = TRUE,row.names=F)
}


sendtoSQL_waveData <- function(wave){
	alldf <- data.frame(Time=c(),uc=c(),vc=c(),wvd=c(),logger=c())
	# wave <- readRDS("2DwaveData.rds")

	for(i in 1:length(wave)){
		data <- wave[[i]]
		# data$Time <- as.character(data$Time)
		data$logger  <- as.numeric(unlist(strsplit(names(wave)[i],"_"))[2])
		names(data)[1] <- "Time"

		# data <- melt(data[,-1],id=c("Time","logger"))
		alldf<- rbind(alldf,data)
	}
	alldf$id=1:nrow(alldf)
	dbWriteTable(conn, "waveData", alldf, overwrite = TRUE,row.names=F)
}

sendtoSQL_tempData <- function(tempData){
	alldf <- data.frame()

	for(i in 1:length(tempData)){
		data <- tempData[[i]]$myTemp
		# data$Time <- as.character(data$Time)
		data$logger  <- as.numeric(unlist(strsplit(names(tempData)[i],"_"))[2])
		names(data)[21] <- "Time"

		# data <- melt(data[,-1],id=c("Time","logger"))
		alldf<- rbind(alldf,data)
	}
	alldf$id=1:nrow(alldf)
	dbWriteTable(conn, "3DtempData", alldf, overwrite = TRUE,row.names=F)
}


# meta_new <- dbReadTable(conn,"loggerInfo")
# d3dMatrix <- matrix(-1,30,20)
# for(i in 1:length(tempData)){
# 	d3dMatrix[i,] <- tempData[[i]]$d3d
# }
# colnames(d3dMatrix) <- paste("depth",1:20,sep="")
# meta_new2 <- cbind(meta_new,d3dMatrix)
