library("ncdf")
library("ggmap")
library(grid)
fname <- "~/Downloads/e201435100.out1.nc"
fid <- open.ncdf(fname)
print(fid)

uc <- get.var.ncdf(fid,"uc")
vc <- get.var.ncdf(fid,"vc")

lat <- get.var.ncdf(fid,"lat")
lon <- get.var.ncdf(fid,"lon")

dataFrame <- data.frame(lat=c(lat),lon=c(lon))
subsampleIndex <- seq(1,nrow(dataFrame),10)
p0 <- ggmap(get_map(location=c(-83.54,41.35,-78.83,42.92)))
print(p0)
# for(i in 1:dim(uc)[3]){

for(i in 1:10){
	dataFrame$uc <- c(uc[,,i])
	dataFrame$vc <- c(vc[,,i])
	
	png(paste(i,".png",sep=""))
	p <- ggplot()+geom_segment(aes(x=lon,y=lat,xend=lon+uc,yend=lat+vc),data=dataFrame[subsampleIndex,],arrow=arrow(length = unit(0.1,"cm"))) + geom_point(aes(Long,Lat),data=meta_B,color="red",size=5)
	print(p)
	dev.off()
}

