require(RMySQL)
require(dplyr)
require(reshape2)
require(zoo)

retriveLoggerData_DO_Temp <- function(loggerIndex,year,groupRange,dataType,timeRange = NULL){
	DOData <- retriveLoggerData(loggerIndex,year,"DO",groupRange,dataType,timeRange = timeRange)
	names(DOData) <- paste("DO",names(DOData),sep="_")

	tempData <- retriveLoggerData(loggerIndex,year,"Temp",groupRange,dataType,timeRange = timeRange)
	names(tempData) <- paste("Temp",names(tempData),sep="_")

	return(cbind(DOData,tempData))
}


retriveLoggerData <- function(loggerIndex,year,var,groupRange,dataType,timeRange=NULL){
	data <- data.frame()
	sqlRes <- getTimeSeriesSQL(loggerIndex,year,var,groupRange,dataType,timeRange = timeRange)
	sql <- sqlRes$sql
	timeFormat <- sqlRes$timeFormat
	# print(timeFormat)
	data <- sqlQuery(sql) %>% dcast(Time~logger,value.var=var)
	data <- zoo(subset(data,select=-Time),order.by=strptime(data$Time,format=timeFormat,tz="GMT"))
	return(data)
}

retriveGeoData <- function(year,position){
	sql <- sprintf("select loggerID,longitude,latitude,bathymetry from loggerInfo where available=1 and loggerPosition='%s' and year = %s",position, year)
	return(sqlQuery(sql))
}


retrivePairLogger <- function(loggerIndex){
	condition <- sprintf("(%s)",paste(loggerIndex,collapse=","))
	sql <- sprintf("select bottom,upper from loggerBottomUpper where bottom in %s or upper in %s",condition,condition)
	return(sqlQuery(sql))
}


sqlQuery <- function (sql) {
	if(nchar(sql)<1){
		stop("wrong sql")  	
  	}
  	conn <- dbConnect(MySQL(), dbname = dbConfig$dbname, username=dbConfig$username, password=dbConfig$password, host=dbConfig$host, port=3306)

  	result <- dbGetQuery(conn,sql)
  	dbDisconnect(conn)
  	return(result)
}

getTimeSeriesSQL <- function(loggerIndex,year,var,groupRange,dataType,timeRange=NULL){
	loggerCondition <- sprintf("(%s)",paste( paste("logger =",loggerIndex),collapse=" OR "))
	if(is.null(timeRange)){
		sqlCondition <- loggerCondition
	}else{
		sqlCondition <- sprintf("%s AND Time >= '%s' AND Time <= '%s' ",loggerCondition,timeRange[1],timeRange[2])
	}

	if(dataType == "STD"){

		if(groupRange == "daily"){
			sql <- sprintf("Select date(Time) as Time, STD(%s) as %s, logger from loggerData_%s where %s Group by date(Time),logger",var,var,year,sqlCondition)
			timeFormat="%Y-%m-%d"
		}else if(groupRange == "hourly"){
			sql <- sprintf("Select DATE_FORMAT(Time,'%%Y-%%m-%%d %%H') as Time, STD(%s) as %s, logger from loggerData_%s where %s Group by DATE_FORMAT(Time,'%%Y-%%m-%%d %%H'),logger",var,var,year,sqlCondition)
			timeFormat <- "%Y-%m-%d %H"
		
		}else{
			 stop("Invalid groupRange type")
		}
	}

	else if(dataType == "AVG"){

		if(groupRange == "daily"){
			sql <- sprintf("Select date(Time) as Time, AVG(%s) as %s, logger from loggerData_%s where %s Group by date(Time),logger",var,var,year,sqlCondition)
			timeFormat <- "%Y-%m-%d"
		}else if(groupRange == "hourly"){
			sql <- sprintf("Select DATE_FORMAT(Time,'%%Y-%%m-%%d %%H') as Time, AVG(%s) as %s, logger from loggerData_%s where %s Group by DATE_FORMAT(Time,'%%Y-%%m-%%d %%H'),logger",var,var,year,sqlCondition)
			timeFormat <- "%Y-%m-%d %H"
		}else{
			stop("Invalid groupRange type")
		}
	}else{
		sql <- sprintf("Select Time, %s, logger from loggerData_%s where %s",var,year,sqlCondition)
		timeFormat="%Y-%m-%d %H:%M:%S"
	}

	return(list(sql=sql,timeFormat=timeFormat))
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
