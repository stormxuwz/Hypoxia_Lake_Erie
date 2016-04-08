source("../config.R")



# setwd("~/Developer/DO_Sensor_data/")
# require("ggplot2")
# require("gridExtra")
# require("reshape2")
# source("./settings.R")
# source("dataPreprocess.R")
# source("chaotic.R")
# source("dynamic_time_warping.R")
# source("dataStatistics.R")
# library("dplyr")
# require("sp")
# source("findBathy.R")

# # Read data
# meta <- readRDS(paste(dataSource,"meta.rds",sep=""))
# DO_raw_list <- readRDS(paste(dataSource,"rawFile.rds",sep=""))

# grid <- read.table("grid_final.csv",sep=",")

# badStationNo <- c("10384448","10523440","10523438","10523444","10384450")

# meta_available <- meta[!(as.character(meta$No) %in% badStationNo),]

# # meta_B <- subset(meta_available,position=="B")  # find the bottom site
# meta_B <- as.data.frame(readRDS("meta_bottom2.rds"))

# # meta_B3 <- subset(meta_available,position=="B3")  # find the bottom site

# stationName=paste("logger",meta_B$No,sep="_")
# rownames(meta_B) <- stationName

# # Rearrange all data
# syncDataList <- dataPreprocessing(FALSE)

# DO_upper <- syncDataList$upperDO
# Temp_upper <- syncDataList$upperTemp
# DO_upper <- additionTime(DO_upper)

# DO_bottom <- syncDataList$bottomDO
# Temp_bottom <- syncDataList$bottomTemp
# DO_bottom <- additionTime(DO_bottom)

# diff_temp <- Temp_upper$logger_10523448-Temp_bottom$logger_10523446
# diff_DO <- DO_upper$logger_10523448-DO_bottom$logger_10523446

# Time <- DO_bottom$Time
# qplot(Time,scale(diff_temp))+geom_point(aes(Time,scale(diff_DO)),color="red")
# # write.csv2(DO_bottom,paste(dataSource,"DO_bottom.csv"))

# # DO_bottom_yday <- aggregateByTime(DO_bottom,"yday")
# # DO_bottom_hour <- aggregateByTime(DO_bottom,"hour")
# # DO_bottom_cumHour <- aggregateByTime(DO_bottom,"cumHour")


# DO_bottom_dailyMean <- aggregateByTime(DO_bottom,"yday",method="mean")
# DO_bottom_dailyVar <- aggregateByTime(DO_bottom,"yday",method="var")

# DO_bottom_hourlyMean <- aggregateByTime(DO_bottom,"cumHour",method="mean")
# DO_bottom_hourlyVar <- aggregateByTime(DO_bottom,"cumHour",method="var")


# qplot(cumHour,value,data=DO_bottom_aggregated_reshaped,color= variable)+geom_line()




# # DO_bottom_hour_reshaped <- as.data.frame(scale(DO_bottom_hour,scale=FALSE,center=FALSE))
# # DO_bottom_hour_reshaped$hour <- DO_bottom_hour$hour

# # DO_bottom_hour_reshaped <- melt(DO_bottom_hour_reshaped,id.vars=c("hour"))


# #### Interpolation ####

# # Function to calculate the distance between each loggers
# distanceMatrix=matrix(0,nrow(meta_B),nrow(meta_B))
# rownames(distanceMatrix) <- stationName
# colnames(distanceMatrix) <- stationName

# for(i in 1:nrow(meta_B)){
#     distanceMatrix[i,] <- spDistsN1(as.matrix(meta_B[,c("Long","Lat")]),as.matrix(meta_B[i,c("Long","Lat")]),TRUE)
#     # Unit is "km"
# }


# for(i in 1:length(stationName)){
#     print(stationName[i])
#     # find the nearest station 
#     nearestDis <- min(distanceMatrix[i,-i])
#     j <- which(distanceMatrix[i,]==nearestDis)
#     # plotStationPlot(DO_bottom_aggregated,"yday",stationName[c(i,j)])
    
#     dtw_alignment <- dynamic_time_warping(DO_bottom_aggregated,stationName[i],stationName[j])
# }

#  plotStationPlot <- function(dataTable,aggrOnName,stationName){
#     dataTable_reshaped <- melt(dataTable[,c(stationName,aggrOnName)],id.vars=c(aggrOnName))
#     qplot(dataTable_reshaped[,aggrOnName],value,data=dataTable_reshaped,color=variable)+geom_line()

#  }

# # read wind direction
#   plotDOTemp <- function(loggerName){
#      if(loggerName %in% names(DO_bottom)){
#          dataTable_DO <- DO_bottom
#          dataTable_Temp <- Temp_bottom
#      }
#      else if(loggerName %in% names(DO_upper)){
#          dataTable_DO <- DO_upper
#          dataTable_Temp <- Temp_upper
#      }
#      else{
#          return(0)
#      }
#      Time <- dataTable_DO$Time
#      DO <- scale(dataTable_DO[,loggerName])
#      Temp <- scale(dataTable_Temp[,loggerName])
     
#      dataTable <- data.frame(Time=Time,DO=DO,Temp=Temp)
     
#      p <- ggplot(data=dataTable)+geom_point(aes(Time,DO),label="DO")+geom_point(aes(Time,Temp),label="Temp")+geom_text()
#      print(p)
#  }

#  changeToCels <- function(dataTable){

#      for logger in c("logger_10384449","logger_")
#      tempSeries <- 
#  }

 

# dataTable <- data.frame(lat=c(lat),lon=c(lon),depth=c(depth))

# retrivePoint <- function(fileName,lon,lat){
#     data=0
#     return(data)
# }

# # qplot(V2,V1,data=grid,z=V3,color=V3)+stat_contour(bins=25)


# # Visualization





# # idw_result_matrix <- idwInterpolation(DO_bottom,grid_final,meta_B)
# # Plot interpolation results
# # DO_bottom_t=meta_B





# for(i in 1:dim(a)[2]){
#     grid$DO=a[,i]
#     DO_bottom_t$DO=as.numeric(DO_bottom[i,-c(1,2)])
#     timeStr=as.character(time_line[i])
#     colorRange=range(range(grid$DO),range(subset(DO_bottom_t,DO>0)$DO))
#     colorRange[1]=colorRange[1]-colorRange[1]%%0.2
#     colorRange[2]=colorRange[2]-colorRange[2]%%0.2+0.2
    
#     p=qmplot(Long,Lat,data=grid,color=DO)+scale_color_gradientn(name="DO",colours = topo.colors(10),limit=colorRange)+geom_point(aes(Long,Lat),data=DO_bottom_t,color="black",alpha=0.5,size=3)
    
#     q=qmplot(Long,Lat,data=DO_bottom_t,color=DO,size=I(5),label=sprintf("%.3f",DO))+scale_color_gradientn(name="DO",colours = topo.colors(10),limit=colorRange)+geom_text(size=4,vjust=-1.5,color="black")
    
#     jpeg(paste("~/Developer/Project_IO/DO_LakeErie/output/Tps/",i,"_DO.jpg",sep=""),height=1.6*200,width=6*200)
#     grid.arrange(p,q,ncol=2,main = timeStr)
#     dev.off()
# }


# require("gridExtra")
# loggerNameList <- c("10384449","10523447","10384443")

# for(loggerName in loggerNameList){
#     df <- DO_raw_list[[loggerName]]
#     df <- df[10:(nrow(df)-10),]
#     df$Time <- as.POSIXct(df$Time)

#     if(max(df$Temp)>50){
#         df$Temp  <- (df$Temp - 32)/1.8
#     }
#     df <- melt(df[,-1],id.vars=c("Time"))
#     qplot(Time,value,data=df)+ facet_wrap(~ variable, ncol= 2,scale="free")
# }


# DOdf <- DO_bottom[,c("Time","logger_10384449","logger_10523445")]
# tempDF <- Temp_bottom[,c("Time","logger_10384449","logger_10523445")]

# tempDF$logger_10384449 <- (tempDF$logger_10384449-32)/1.8

# DOdf<- melt(DOdf,id.vars=c("Time"))
# tempDF<- melt(tempDF,id.vars=c("Time"))

# DOdf$type <- "DO"
# tempDF$type <- "Temp"

# final <- rbind(DOdf,tempDF)

# qplot(Time,value,data=DOdf,color=variable)+facet_wrap(. ~ type,nrow=2)

# p1 <- qplot(Time,value,data=DOdf,color=variable,size=I(1))+labs(y = "DO (mg/L)",colour = "loggers")
# p2 <- qplot(Time,value,data=tempDF,color=variable,size=I(1))+labs(y = "Temp (C)",colour = "loggers")

# grid.arrange(p1,p2,nrow=2)





