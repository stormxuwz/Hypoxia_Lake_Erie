library(RMySQL)

conn <- dbConnect(MySQL(), dbname = "DO", username="root", password="XuWenzhaO", host="127.0.0.1", port=3306)

sql_readData <- "Select DATE_FORMAT(Time,'%Y-%m-%d %H') as Time, VARIANCE(DO) as DO_var, AVG(DO) as DO_avg, logger as logger from loggerData where logger=10528849 or logger=10523447 Group by date(Time),hour(Time),logger"
res <- dbGetQuery(conn,sql_readData)
res$logger <- as.factor(res$logger)
res$Time <- as.POSIXct(res$Time,format="%Y-%m-%d %H")
qplot(Time,DO_var,data=res,color=logger,geom="line")




sql_readData1 <- "Select Time as Time, DO as DO, logger from loggerData where logger=10528849 or logger=10523447"
res1 <- dbGetQuery(conn,sql_readData1)
res1$Time <- as.POSIXct(res1$Time)
res1$logger <- as.factor(res1$logger)

p <- ggplot()+geom_point(aes(Time,DO,color=logger),data=res1)

sql_readData2 <- "Select Time as Time, uc, vc, wvd, logger from waveData where logger=10528849 or logger=10523447"
res2 <- dbGetQuery(conn,sql_readData2)
res2$Time <- as.POSIXct(res2$Time)
res2$logger <- as.factor(res2$logger)


sql_readData3 <- "select Time as Time, V20, logger from 3DtempData where logger=10528849"
res3 <- dbGetQuery(conn,sql_readData3)
res3$Time <- as.POSIXct(res3$Time)
qplot(Time,V20,data=res3)



loggerSummary <- function(logger){
	# function to give plot summary
	sql_readData_logger <- paste("Select Time, DO, Temp, logger from loggerData where logger=",logger," AND Time < '2014-10-09'  AND Time > '2014-06-25'")
	res1 <- dbGetQuery(conn,sql_readData_logger)
	res1$Time <- as.POSIXct(res1$Time)
	p1 <- ggplot()+geom_point(aes(Time,DO),data=res1)
	p2 <- ggplot()+geom_point(aes(Time,Temp),data=res1)

	sql_readData_wave <- paste("Select uc, vc, Time, logger from waveData where logger=",logger,"AND Time < '2014-10-09' AND Time > '2014-06-25'")
	res2 <- dbGetQuery(conn,sql_readData_wave)
	res2 <- melt(res2,id=c("Time","logger"))
	res2$Time <- as.POSIXct(res2$Time)
	p3 <- ggplot()+geom_point(aes(Time,value,color=variable),data=res2)

	sql_readData_temp <- paste("Select V20, V17, V1,Time,logger from 3DtempData where logger=",logger,"AND Time < '2014-10-09' AND Time > '2014-06-25'")
	res3 <- dbGetQuery(conn,sql_readData_temp)
	res3 <- melt(res3,id=c("Time","logger"))
	res3$Time <- as.POSIXct(res3$Time)
	p4 <- ggplot()+geom_point(aes(Time,value,color=variable),data=res3)

	png(paste(logger,"_summary.png",sep=""),height=800,width=800)
	grid.arrange(p1,p3,p2,p4,nrow=2)
	dev.off()
}

loggerSummary_2015 <- function(logger){
	sql_readData_logger <- paste("Select Time, DO, Temp, logger from loggerData where logger=",logger)
	res1 <- dbGetQuery(conn,sql_readData_logger)
	res1$Time <- as.POSIXct(res1$Time)
	p1 <- ggplot()+geom_point(aes(Time,DO),data=res1)
	p2 <- ggplot()+geom_point(aes(Time,Temp),data=res1)

	png(paste("./2015/",logger,"_summary.png",sep=""),height=800,width=800)
	grid.arrange(p1,p2,nrow=2,top =paste("logger",logger))
	dev.off()
}

plotLogger <- function(conn,dataType,loggerList){
	conn <- dbConnect(MySQL(), dbname = "DO2015", username="root", password="XuWenzhaO", host="127.0.0.1", port=3306)
	tmp <- paste("logger =",loggerList)
	tmp <- paste(tmp,collapse=" OR ")
	sql <- paste("Select",dataType,",Time,logger from loggerData where",tmp)
	res <- dbGetQuery(conn,sql)
	res$Time <- as.POSIXct(res$Time)
	res$logger <- as.factor(res$logger)
	names(res)[1]="value"
	qplot(Time,value,data=res,color=logger)
}


library(RMySQL)
conn <- dbConnect(MySQL(), dbname = "DO2014", username="root", password="XuWenzhaO", host="127.0.0.1", port=3306)
meta <- dbReadTable(conn,"loggerInfo")

for(i in 1:nrow(meta)){
	if(meta[i,"available"]==1){
		loggerSummary(meta[i,"loggerID"])
	}
}

