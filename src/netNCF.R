# function to read ncdf files

library("ncdf")

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


# library("ncdf")
# library("ggmap")
# library(grid)
# fname <- "~/Downloads/e201435100.out1.nc"
# fid <- open.ncdf(fname)
# print(fid)

# uc <- get.var.ncdf(fid,"uc")
# vc <- get.var.ncdf(fid,"vc")

# lat <- get.var.ncdf(fid,"lat")
# lon <- get.var.ncdf(fid,"lon")

# dataFrame <- data.frame(lat=c(lat),lon=c(lon))
# subsampleIndex <- seq(1,nrow(dataFrame),10)
# p0 <- ggmap(get_map(location=c(-83.54,41.35,-78.83,42.92)))
# print(p0)
# # for(i in 1:dim(uc)[3]){

# for(i in 1:10){
#     dataFrame$uc <- c(uc[,,i])
#     dataFrame$vc <- c(vc[,,i])
    
#     png(paste(i,".png",sep=""))
#     p <- ggplot()+geom_segment(aes(x=lon,y=lat,xend=lon+uc,yend=lat+vc),data=dataFrame[subsampleIndex,],arrow=arrow(length = unit(0.1,"cm"))) + geom_point(aes(Long,Lat),data=meta_B,color="red",size=5)
#     print(p)
#     dev.off()
# }

