require(ggplot2)
require(dygraphs)
require(reshape2)
require(leaflet)

plot_value<- function(data,type="ggplot"){
	# data is a zoo dataframe
	# change data type
	if(type == "ggplot"){
		data <- as.data.frame(data)
		data$time <- strptime(rownames(data))
		data <- melt(data,id.vars="time")
		# ggplot(data)+geom_point(aes(time,value,color=))
	}
	else if (type=="dygrphs"){
		dygraph(data)
	}
	else{
		stop ("invalid plot type")
	}

}


plot_spatial <- function(data){

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
	p <- p %>% dyAxis("y", label = "DO(mg/L)") %>% dyAxis("y2", label = "Temperature (C)", independentTicks = TRUE) %>% dyRangeSelector()
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

