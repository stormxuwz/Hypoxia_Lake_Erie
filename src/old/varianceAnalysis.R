library(RMySQL)
library(gstat)
library(fields)
library(sp)
meta <- dbReadTable(conn,"loggerInfo")

conn <- dbConnect(MySQL(), dbname = "DO2014", username="root", password="XuWenzhaO", host="127.0.0.1", port=3306)
meta <- dbReadTable(conn,"loggerInfo")

getSeries <- function(sql){
	res <- dbGetQuery(conn,sql)
	return(res)
}

meta <- subset(meta,available==1 & loggerPosition=="B")
coordinates(meta)=~longitude+latitude


distanceMatrix=matrix(0,nrow(meta),nrow(meta))
rownames(distanceMatrix) <- meta$loggerID
colnames(distanceMatrix) <- meta$loggerID

# distanceMatrix <- sapply(spDistsN1,meta,pts=meta)
allDist <- c()
loggerWithinDistance <- list()
for(i in 1:nrow(meta)){
    distanceMatrix[i,] <- spDistsN1(meta,meta[i,],TRUE) # kilometers will return
    distanceToOthers <- distanceMatrix[i,]
    allDist <- c(allDist,distanceToOthers)

	loggerWithinDistance[[as.character(meta$loggerID[i])]] <- meta$loggerID[which(distanceToOthers>30 & distanceToOthers<40)]
}

submeta <- subset(meta,loggerID %in% unlist(loggerWithinDistance))
leaflet() %>% addMarkers(data=submeta) %>% addTiles() %>% addPopups(data=submeta,popup = ~paste(loggerID))


tmp <- "logger = "
var <- "temp"
sql <- sprintf("Select DATE_FORMAT(Time,'%%Y-%%m-%%d %%H') as Time, AVG(%s) as %s, logger from loggerData where %s Group by DATE_FORMAT(Time,'%%Y-%%m-%%d %%H'),logger",var,var,tmp)


pairComb <- combn(1:nrow(meta),2)


getHourlyAgg<- function(conn,loggerID,aggType,startTime="2014-07-01",endTime="2014-10-01"){
	sql <- sprintf("Select DATE_FORMAT(Time,'%%Y-%%m-%%d %%H') as Time, %s(DO) as DO, logger from loggerData where logger = %s and Time > '%s' and Time < '%s' Group by DATE_FORMAT(Time,'%%Y-%%m-%%d %%H'),logger",aggType,as.character(loggerID),as.character(startTime),as.character(endTime))
	print(sql)
	res <- dbGetQuery(conn,sql)
	return(res)
}


corTable <- data.frame()

for(i in 1:ncol(pairComb)){
	j <- pairComb[1,i]
	k <- pairComb[2,i]
	pairdist <- spDistsN1(meta[j,],meta[k,],TRUE)
	loggerID_j <- meta$loggerID[j]
	loggerID_k <- meta$loggerID[k]

	j_data_avg <- getHourlyAgg(conn,loggerID_j,"AVG")$DO
	k_data_avg <- getHourlyAgg(conn,loggerID_k,"AVG")$DO

	j_data_std <- getHourlyAgg(conn,loggerID_j,"STD")$DO
	k_data_std <- getHourlyAgg(conn,loggerID_k,"STD")$DO

	corTable <- rbind(corTable,data.frame(id_1=j,id_2=k,loggerID_1=loggerID_j,loggerID_2=loggerID_k,dist=pairdist,avgCor=cov(j_data_avg,k_data_avg),stdCor=cov(j_data_std,k_data_std)))
}

# head(corTable)
saveRDS(corTable,"covTable.rds")
library(dplyr)
intervalIndex <- findInterval(corTable$dist,seq(0,max(corTable$dist),10))
corTable$bins <- intervalIndex
d <- group_by(corTable,bins) %>% summarise(dist=mean(dist),avgCor=mean(avgCor),count=n())




