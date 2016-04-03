# Create Bathymetry

bathyGrid <- function(sensorMeta){
	bathymetry_depth<-read.table("erie_central.asc",header=F) # read depths files
	bathymetry_shape<-dim(bathymetry_depth)
	bathymetry_depth<-data.matrix(bathymetry_depth)
	grid<-expand.grid(Lat=seq(from=41.499583333333+(0.000833333333)*1200,by=-0.000833333333*12,length.out=101),Long=seq(from=-82.500416666667,by=0.000833333333*12,length.out=201))
	bathymetry_depth_rough<-bathymetry_depth[seq(1,1201,12),seq(1,2401,12)]  # take a coarser resolution
    grid$depth<-as.vector(bathymetry_depth_rough)
    return(rawGrid=grid,newSensorMeta=sensorMeta)
}

createFaceIndex <- function(rawGrid,convelHullModel){
    # create face index for visualization
    inConvexHull <- function(i){
        if(sum(predict(convexHull,rawGrid[faceVerticeIndex[i,],c("Lat","Long")]))==3)
            return(i)
    }

    nX=unique(rawGrid$Long) # number of columns, 201
    nY=unique(rawGrid$Lat) # number of rows, 101
    
    faceVerticeIndex=matrix(0,(nX-1)*(nY-1)*2,3)
    i <- 0
    for(x in 1:(nX-1))
    {
        for(y in 1:(nY-1)){ 
            a <- (x-1)*nY+y
            b <- (x-1)*nY+y+1
            c <- x*nY+y
            d <- x*nY+y+1
            faceVerticeIndex[i,] <- c(a,b,c)
            faceVerticeIndex[i+1,] <- c(b,c,d)
            i <- i+2
        }
    }

    outOfBoundaryFaceIndex <- apply(1:nrow(raw_face_vertices_index),test_valid)
    
    return(faceVerticeIndex,outOfBoundaryFaceIndex)
}

createBathyGrid <- function(sensorMeta){
    require(dismo,quietly=TRUE)

    convexHullModel<-convHull(sensorMeta[,c("Lat","Long")])
    bathyInfo <- bathyGrid(sensorMeta)
    rawGrid <- bathyInfo$rawGrid
    faces <- createFaceIndex(rawGrid,convexHullModel)
    
    rawGrid$convexIndex <- predict(convexHull,rawGrid[,c("Lat","Long")])
    rawGrid$remainIndex <- 0
    rawGrid[rawGrid$convexIndex>0,"remainIndex"] <- c(1:sum(rawGrid$convexIndex))
   
    # Create Faces
    faces <- createFaceIndex(rawGrid,convexHullModel)
    faceVerticeIndex <- faces$faceVerticeIndex
    omitFaceIndex <- faces$outOfBoundaryFaceIndex
    faceVerticeIndex <- faceVerticeIndex[-omitFaceIndex,]

    apply(1:nrow(faceVerticeIndex),function(i) {faceVerticeIndex[i,]<<-rawGrid[faceVerticeIndex[i,],"remainIndex"]} )

    return(list(grid=subset(rawGrid,convexIndex>0),faceVerticeIndex=faceVerticeIndex))
}

# write.table(face_vertices_index_final,file="face_index.csv",col.names=F,row.names=F,sep=",")
# write.table(grid_final[,c("Lat","Long","depth")],file="grid_final.csv",col.names=F,row.names=F,sep=",")






