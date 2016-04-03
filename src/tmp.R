
otherDataPreprocessing <- function(folder="../Project_IO/DO_LakeErie/OtherData/"){
    dataList_45130 <- c("C45132_2014-6","C45132_2014-7","C45132_2014-8","C45132_2014-9","C45132_2014-10")
    data_45132 <- read.csv(paste(folder,dataList_45130[1],".csv",sep=""))
    for(i in 2:length(dataList_45130)){
        data_45132=rbind(data_45132,read.csv(paste(folder,dataList_45130[i],".csv",sep="")))
    }
    data_45132$OBS_DATE <- as.POSIXct(as.character(data_45132$OBS_DATE),"%Y/%m/%d %H:%M",tz="EST")
    saveRDS(data_45132,"data_45132.rds")
    
    
}





readAllNCFile <- function(dataType="2d"){
    require(ncdf)
    library(abind)
    if(dataType=="2d"){
        fileNameList <- c("e201415100.out1.nc","e201420100.out1.nc","e201425100.out1.nc")
        ncFile = open.ncdf(paste("../Project_IO/DO_LakeErie/OtherData/",fileNameList[1],sep=""))
        print(ncFile)
        #     [1] "file ../Project_IO/DO_LakeErie/OtherData/e201420100.out1.nc has 4 dimensions:"
        #     [1] "time   Size: 1200"
        #     [1] "ny   Size: 87"
        #     [1] "nx   Size: 193"
        #     [1] "nsigma   Size: 20"
        #     [1] "------------------------"
        #     [1] "file ../Project_IO/DO_LakeErie/OtherData/e201420100.out1.nc has 14 variables:"
        #     [1] "float ci[nx,ny,time]  Longname:Ice concentration Missval:-99999"
        #     [1] "float depth[nx,ny]  Longname:Bathymetry  Missval:1e+30"
        #     [1] "float eta[nx,ny,time]  Longname:Height Above Model Sea Level Missval:-99999"
        #     [1] "float hi[nx,ny,time]  Longname:Ice thickness Missval:-99999"
        #     [1] "float lat[nx,ny]  Longname:Latitude Missval:1e+30"
        #     [1] "float lon[nx,ny]  Longname:Longitude Missval:1e+30"
        #     [1] "float sigma[nsigma]  Longname:Sigma Stretched Vertical Coordinate at Nodes Missval:1e+30"
        #     [1] "float uc[nx,ny,time]  Longname:Eastward Water Velocity at Surface Missval:-99999"
        #     [1] "float ui[nx,ny,time]  Longname:Ice u-velocity Missval:-99999"
        #     [1] "float vc[nx,ny,time]  Longname:Northward Water Velocity at Surface Missval:-99999"
        #     [1] "float vi[nx,ny,time]  Longname:Ice v-velocity Missval:-99999"
        #     [1] "float wvd[nx,ny,time]  Longname:Wave Direction Missval:-99999"
        #     [1] "float wvh[nx,ny,time]  Longname:Significant Wave Height Missval:-99999"
        #     [1] "float wvp[nx,ny,time]  Longname:Wave Period Missval:-99999"
        
        lat <- get.var.ncdf(ncFile,"lat")
        lon <- get.var.ncdf(ncFile,"lon")
        uc <- get.var.ncdf(ncFile,"uc") # 
        vc <- get.var.ncdf(ncFile,"vc")
        wvd <- get.var.ncdf(ncFile,"uc")
        bathy <- get.var.ncdf(ncFile,"depth")
        
        close.ncdf(ncFile)
        
        for(fname in fileNameList[-1]){
            ncFile = open.ncdf(paste("../Project_IO/DO_LakeErie/OtherData/",fname,sep=""))
            uc <- abind(uc,get.var.ncdf(ncFile,"uc"),along=3)  
            vc <- abind(vc,get.var.ncdf(ncFile,"vc"),along=3)  
            wvd <- abind(wvd,get.var.ncdf(ncFile,"wvd"),along=3)  
        }
        
        return(list(lat=lat,lon=lon,uc=uc,vc=vc,wvd=wvd))
    }
        
        
    else{
        fileNameList <- c("e201415100.out3.nc","e201420100.out3.nc","e201425100.out3.nc")
        ncFile = open.ncdf(paste("../Project_IO/DO_LakeErie/OtherData/",fileNameList[1],sep=""))
        print(ncFile)

        # [1] "file ../Project_IO/DO_LakeErie/OtherData/e201420100.out3.nc has 4 dimensions:"
        # [1] "nsigma   Size: 20"
        # [1] "ny   Size: 87"
        # [1] "nx   Size: 193"
        # [1] "time   Size: 400"  # 3 hours interval
        # [1] "------------------------"
        # [1] "file ../Project_IO/DO_LakeErie/OtherData/e201420100.out3.nc has 8 variables:"
        # [1] "float d3d[nx,ny,nsigma]  Longname:3D Depth at Nodes Missval:1e+30"
        # [1] "float depth[nx,ny]  Longname:Bathymetry  Missval:1e+30"
        # [1] "float lat[nx,ny]  Longname:Latitude Missval:1e+30"
        # [1] "float lon[nx,ny]  Longname:Longitude Missval:1e+30"
        # [1] "float sigma[nsigma]  Longname:Sigma Stretched Vertical Coordinate at Nodes Missval:1e+30"
        # [1] "float temp[nx,ny,nsigma,time]  Longname:Temperature Missval:-99999"
        # [1] "float u[nx,ny,nsigma,time]  Longname:Eastward Water Velocity Missval:-99999"
        # [1] "float v[nx,ny,nsigma,time]  Longname:Northward Water Velocity Missval:-99999"
        
        lat <- get.var.ncdf(ncFile,"lat")
        lon <- get.var.ncdf(ncFile,"lon")
        bath <- get.var.ncdf(ncFile,"depth")
        
        BottomTemperature <- get.var.ncdf(ncFile,"temp")[,,20,]
        BottomU <- get.var.ncdf(ncFile,"u")[,,20,]
        BottomV <- get.var.ncdf(ncFile,"v")[,,20,]
        
        for(fname in fileNameList[-1]){
            ncFile = open.ncdf(paste("../Project_IO/DO_LakeErie/OtherData/",fname,sep=""))
            BottomU <- abind(BottomU,get.var.ncdf(ncFile,"u")[,,20,],along=3)  
            BottomV <- abind(BottomV,get.var.ncdf(ncFile,"v")[,,20,],along=3)  
            BottomTemperature <- abind(BottomTemperature,get.var.ncdf(ncFile,"temp")[,,20,],along=3) 
        }
    
    }
}

startPoint<-as.POSIXlt("2014-6-1 00:00:00",tz="EST")
endPoint<-as.POSIXlt("2015-1-1 00:00:00",tz="EST")
timeIndex<-seq(startPoint,by=3600,length.out=3600)

windData <- readAllNCFile()
saveRDS(windData,"windData.rds")

plotWaveDirection <- function(i){
    require(ggmap)
    yday <- unclass(as.POSIXlt(DO_bottom$Time[i]))$yday
    hour <- unclass(as.POSIXlt(DO_bottom$Time[i]))$hour
    row <- (yday-151)*24+hour
    
    dataFrame <- data.frame(lat=c(windData$lat),lon=c(windData$lon))
    
    # p0 <- ggmap(get_map(location=c(-83.54,41.35,-78.83,42.92)))
    dataFrame$uc <- c(windData$uc[,,row])
    dataFrame$vc <- c(windData$vc[,,row])
    subsampleIndex <- seq(1,nrow(dataFrame),10)
    p <- ggplot()+geom_segment(aes(x=lon,y=lat,xend=lon+uc,yend=lat+vc),data=dataFrame[subsampleIndex,],arrow=arrow(length = unit(0.1,"cm"))) + geom_point(aes(Long,Lat),data=meta_B,color="red",size=5)+ggtitle(as.character(DO_bottom$Time[i]))+xlim(c(-84,-78))+ylim(c(41.25,43.25))
    print(p)
}

trace.animate <- function() {
    lapply(seq(100,10000,6), function(i) {
        plotWaveDirection(i)
    })
}

saveGIF(trace.animate(), interval = .3, movie.name="trace.gif")
