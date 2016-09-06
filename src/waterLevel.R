# function to read water level data
# not useful

readWaterLevel <- function(stationList = c(9063020,9063028,9063038,9063053,9063063,9063079),year = 2014){
	res <- data.frame()
	for(stationNo in stationList){
		tmp <- read.csv(sprintf("../input/OtherData/waterLevel_%d/CO-OPS__%d__hr.csv",year,stationNo))[,1:2]
		tmp$station <- paste(stationNo,"_station",sep="")
		res <- rbind(res,tmp)

		# names(tmp)[2] <- paste(stationNo,"_station",sep="")
		# data <- merge(data,tmp[,1:2],by = "Date.Time",all = TRUE)
	}
	t <- res[,"Date.Time"]
	res[,"Date.Time"]<- as.POSIXct(t,format = "%Y-%m-%d %H:%M",tz = "GMT")
	
	res <- dcast(res,Date.Time~station,value.var = "Water.Level")
	return(res)
}

combineWaterLevelWithDOData <- function(DOdata,waterLevelData){
	n <- nrow(DOdata)
	waterLevelFirstIndex <- which(waterLevelData$Date.Time == index(DOdata)[1])
	subWaterLevel <- waterLevelData[waterLevelFirstIndex:(waterLevelFirstIndex+n-1),]
	# print(waterLevelFirstIndex)
	# subWaterLevel <- zoo(subWaterLevel,order.by = subWaterLevel$Date.Time)
	time <- index(DOdata)
	DOdata <- as.data.frame(DOdata)
	
	DOdata <- cbind(DOdata,subWaterLevel[,-1])
	DOdata <- zoo(DOdata,order.by = time)
	return(DOdata)
}