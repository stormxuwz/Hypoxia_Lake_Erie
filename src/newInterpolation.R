
interpolation <- function(year,method="idw",position="B",area="center"){
	conn <- dbConnect(MySQL(), dbname = paste("DO",year,sep=""), username="root", password="XuWenzhaO", host="127.0.0.1", port=3306)
	loggerInfo <- dbReadTable(conn,"loggerInfo")
	loggerInfo$loggerID <- as.character(loggerInfo$loggerID)
	loggerInfo <- subset(loggerInfo,available==1 & loggerPosition==position)
	
	if(area=="center"){
		centerLogger <- c("10523437","10523442","10528846","10523438","10528847","10523446","10523449","10523450","10523439")
		loggerInfo  <- subset(loggerInfo,loggerID %in% centerLogger)
		sql <- "select avg(DO) as DO,logger as loggerID,DATE_FORMAT(Time,'%Y-%m-%d') as Time from loggerData where logger in (select loggerID from loggerInfo where loggerPosition= 'B' and available=1) group by logger, date(Time)"
		
		DO_data <- dbGetQuery(conn,sql)
		DO_data <- subset(DO_data,loggerID %in% centerLogger)
	}
	else{
		sql <- "select avg(DO) as DO,logger as loggerID,DATE_FORMAT(Time,'%Y-%m-%d') as Time from loggerData where logger in (select loggerID from loggerInfo where loggerPosition= 'B' and available=1) AND Time > '2015-7-23' AND Time < '2015-9-20' group by logger, date(Time)"

		DO_data <- dbGetQuery(conn,sql)
	}


	folder <- paste("./",year,"_",area,sep="")
	dir.create(folder, showWarnings = FALSE)

	grid <- createGrid(loggerInfo)

	if(method=="idw"){
		idwInterpolation(DO_data,loggerInfo,grid,folder)
	}
	
}


idwInterpolation <- function(DO_data, loggerInfo, grid, folder="./"){
	require(gstat)
	require(ggmap)
	times <- unique(DO_data$Time)
	map <- qmplot(longitude,latitude,data=loggerInfo,source="google",zoom=8)

	coordinates(grid)=~longitude+latitude
	
	for(time in times){
		subData <- subset(DO_data,Time==time)
		subData <- merge(subData,loggerInfo)
		coordinates(subData)=~longitude+latitude

		grid$DO <- idw(DO~1,subData,grid)$var1.pred

		# plot data
		p<-map+geom_tile(aes(longitude,latitude,fill=DO),data=as.data.frame(grid))+scale_fill_gradient2(low="red",high="cyan",limits=c(0, 15),mid="yellow",midpoint=7,space = "Lab")+ggtitle(time)

		png(paste(folder,"/",time,".png",sep=""))
		print(p)
		dev.off()
	}

}


createGrid <- function(loggerInfo){
	require("dismo")

	longitudeRange <- range(loggerInfo$longitude)
	latitudeRange <- range(loggerInfo$latitude)
	grid <- expand.grid(longitude=seq(longitudeRange[1],longitudeRange[2],by=0.01),latitude=seq(latitudeRange[1],latitudeRange[2],by=0.02))

	convexHullModel<-convHull(loggerInfo[,c("longitude","latitude")])
	convexIndex <- predict(convexHullModel,grid)
	grid <- subset(grid,convexIndex==1)
	return(grid)
}

interpolation(2015,area="whole")
