require(ggplot2)
require(dygraphs)
require(reshape2)
require(leaflet)
require(ggmap)


plot_value<- function(data,label="value",type="dygrphs",outlierSeries = NULL){
	# data is a zoo dataframe
	
	if(type == "ggplot"){
		data <- as.data.frame(data)
		data$time <- strptime(rownames(data))
		data <- melt(data,id.vars="time")
		# ggplot(data)+geom_point(aes(time,value,color=))
	}
	else if (type=="dygrphs"){
		
		if(!is.null(outlierSeries)){
			# adding outlier series
			print("Plot with outliers")
			
			data <- combindWithOutlier(data,outlierSeries)
			p <- dygraph(data) %>% 
				dyRangeSelector(retainDateWindow=TRUE) %>% 
				dyAxis("y", label = label) %>% 
				dyHighlight(highlightSeriesOpts = list(strokeWidth = 3))
				#dySeries("outlier", color = "red")
		}else{
			p <- dygraph(data) %>% 
				dyRangeSelector(retainDateWindow=TRUE) %>% 
				dyAxis("y", label = label)
		}

		
	}
	else{
		stop ("invalid plot type")
	}
}

plot_spatial_matrix <- function(resultMatrix,outputFolder,prefix = NULL){
	# resultMatrix is a matrix, each raw represent a point on the map

	lonRange <- range(resultMatrix$longitude)
	latRange <- range(resultMatrix$latitude)

	bbox <- make_bbox(lonRange,latRange,f = 0.1)
	myMap <- get_map(location=bbox, source="osm",crop=FALSE)
	p <- ggmap(myMap)

	for(i in 1:104){
		subDf <- resultMatrix[,c(as.character(i),"longitude","latitude")]
		names(subDf)[1] <- "pred"
		q <- p+geom_tile(aes(longitude,latitude,fill = pred),data = subDf)+scale_fill_gradient(name = "DO (mg/L)",low = "red",high = "cyan",limit = c(0,15))
		# +geom_point(aes(longitude,latitude),data =locationInfo,size = I(3))+ggtitle(df_time)
		png(paste(outputFolder,prefix,"_",i,"_",".png",sep=""))
		print(q)
		dev.off()
	}
}


plot_spatial <- function(dataList,locationInfo,outputFolder){
	# dataList is a list which contains dataframes that contains latitude, longitude and interpolation values for each time
	lonRange <- range(unlist(lapply(dataList,function(i) range(i$longitude))))
	latRange <- range(unlist(lapply(dataList,function(i) range(i$latitude))))

	bbox <- make_bbox(lonRange,latRange,f = 0.1)
	myMap <- get_map(location=bbox, source="osm",crop=FALSE)
	p <- ggmap(myMap)
	for(i in 1:length(dataList)){
		print(i)
		df <- dataList[[i]]
		df_time <- attributes(df)$time
		q <- p+geom_tile(aes(longitude,latitude,fill = pred),data = df)+scale_fill_gradient(name = "DO (mg/L)",low = "red",high = "cyan",limit = c(0,15))+geom_point(aes(longitude,latitude),data =locationInfo,size = I(3))+ggtitle(df_time)
		
		png(paste(outputFolder,i,"_",df_time,".png",sep=""))
		print(q)
		dev.off()
	}

}


plot_DO_temp <- function(data){
	TempIndex <- grepl("Temp",names(data))
	plotStr <- c("dygraph(data)")
	print(TempIndex)
	for(i in 1:ncol(data)){
		if(TempIndex[i]){
			plotStr <- c(plotStr,sprintf("dySeries(names(data)[%d], axis = 'y2')",i))
		}
	}
	plotStr <- paste(plotStr,collapse="%>%")
	p <- eval(parse(text=plotStr))
	p <- p %>% dyAxis("y", label = varUnit[["DO"]]) %>% dyAxis("y2", label = varUnit[["Temp"]], independentTicks = TRUE) %>% dyRangeSelector(retainDateWindow=TRUE)
	return(p)
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

