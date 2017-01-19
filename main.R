rm(list = ls())
setwd("/Users/wenzhaoxu/Developer/Hypoxia/Hypoxia_Lake_Erie")
source("./src/database.R")
source("config.R")
source("./src/plot.R")
source("./src/interpolation.R")

# analysis <- function(year){
#year <- 2014
#bottomLogger <- retriveGeoData(year,"B")
# upperLogger <- retriveGeoData(year,"B3")

# loggerData_DO_Temp <- retriveLoggerData_DO_Temp(upperLogger$loggerID,year,"hourly","AVG",c("2014-06-15","2014-10-06"))

# plot_DO_temp(loggerData_DO_Temp)


analysis <- function(year,timeAggType, r){
	# function to do interpolation
	# =============================

	# the logger info
	loggerInfo <- retriveGeoData(year,"B") 

	# the data used to interpolation
	data <- retriveLoggerData(loggerInfo$loggerID,year,"DO",timeAggType,"AVG",transform = TRUE) %>% na.omit()  # remove the data

	# assign the output folder
	metaFolder <<- sprintf("../meta_%d_%s_%d/", year,timeAggType,r)
	outputFolder <<- sprintf("../output_%d_%s_%d/", year,timeAggType,r)

	# create the output folder
	dir.create(file.path("..", sprintf("meta_%d_%s_%d", year,timeAggType,r)), showWarnings = FALSE)
	dir.create(file.path("..", sprintf("output_%d_%s_%d", year,timeAggType,r)), showWarnings = FALSE)
	
	crossValidation(data,loggerInfo, "svd",rList=c(r))

	basis_hypoxiaExtent <- interpolation_main(data,loggerInfo,"basis",basis_method = "svd", simNum = 20, intMethod = "baye", r = r)

	saveRDS(basis_hypoxiaExtent,paste(metaFolder,"basis_hypoxiaExtent.rds",sep = ""))
	
}

summary_plot <- function(year,timeAggType,r){
	metaFolder <<- sprintf("../meta_%d_%s_%d/", year,timeAggType,r)
	outputFolder <<- sprintf("../output_%d_%s_%d/", year,timeAggType,r)

	loggerInfo <- retriveGeoData(year,"B")
	data <- retriveLoggerData(loggerInfo$loggerID,year,"DO",timeAggType,"AVG",transform = TRUE) %>% na.omit()  # remove the data
	
	timeIndex  <- index(data)
	finalPrediction_basis <- readRDS(paste(outputFolder,"final_prediction_basis.rds",sep =""))
	finalPrediction_idw <- readRDS(paste(outputFolder,"final_prediction_IDW.rds",sep = ""))

	basis_hypoxiaExtent <- summaryHypoxia(finalPrediction_basis,timeIndex)
	IDW_hypoxiaExtent <- summaryHypoxia(finalPrediction_idw,timeIndex)

	trendPrediction <- readRDS(paste(metaFolder,"trend_prediction_basis.rds",sep = ""))
	coef_df <- readRDS(paste(metaFolder,"coef_df.rds",sep = ""))
	DO_res <- readRDS(paste(metaFolder,"DO_res.rds",sep = ""))
	basis <- readRDS(paste(metaFolder,"basis.rds",sep = "")) %>% zoo(order.by = timeIndex)

	
	# plot cross validation cvResults
	# cvResults <- readRDS(file = paste(metaFolder,"/cvResults.rds",sep = ""))
	# for(loggerID in loggerInfo$loggerID){
	# 	colName <- paste(r,loggerID,sep = "_")
	# 	cvRes <- cvResults[[colName]]
	# 	names(cvRes) <- c("CV Prediction","CV Prediction_trend","Observation")
	# 	cvRes$time <- as.POSIXct(timeIndex)
	# 	cvRes <- melt(cvRes,id.vars = "time")
		
	# 	png(paste(outputFolder,loggerID,"cv_res.png",sep = ""),width=500,height=350,res=72)
	# 	print(qplot(time,value,data = cvRes,color = variable)+ylab("DO(mg/L)")+scale_color_manual(loggerID,values=c("red","black")))
	# 	dev.off()
	# }

	grid <- trendPrediction$grid
	# plot interpolation range
	# bbox <- make_bbox(range(grid$longitude),range(grid$latitude),f = 0.1)
	# myMap <- get_map(location=bbox, source="osm",crop=FALSE)
	# p <- ggmap(myMap)
	# p <- p+geom_tile(aes(longitude,latitude),data = subset(grid,convexIndex == 1), alpha = I(0.5))
	# p <- p+geom_point(aes(longitude,latitude),data = loggerInfo)
	
	# png(paste(outputFolder,"interpolationArea.png",sep = ""))
	# print(p)
	# dev.off()
	# calculate the spatio-temporal variogram
	# stvgm <- st_variogram(DO_res,loggerInfo)
	# print(plot(stvgm$vgmModel, map=FALSE))

	area <- attr(grid,"totalArea")
	n <- sum(grid$convexIndex)
	names(basis_hypoxiaExtent)  <- paste(names(basis_hypoxiaExtent),"_basis",sep = "")
	names(IDW_hypoxiaExtent)  <- paste(names(IDW_hypoxiaExtent),"_idw",sep = "")
	

	hypoxiaExtent <- cbind(data.frame(basis_hypoxiaExtent), data.frame(IDW_hypoxiaExtent)) %>% zoo(order.by = timeIndex)
	
	for(hypoxiaType in c("less0","less2","less4")){
		krigNames <- paste(hypoxiaType,c("Upper_basis","Mean_basis","Lower_basis"),sep = "")
		IDWNames <- paste(hypoxiaType,c("_idw"),sep = "")

		p <- dygraph(hypoxiaExtent[,c(krigNames,IDWNames)]/n*area, main = sprintf("Total Area %.2f km^2",area)) %>%
		dyAxis("y", label = "Hypoxia Area (km^2)",valueRange = c(0, area+10)) %>%
		dySeries(krigNames, label = "Spatiao-Temporal Kriging") %>% 
		dySeries(IDWNames, label = "IDW Interpolation") %>% 
		dyRangeSelector(height = 20)
		saveWidget(p, file = paste(outputFolder,hypoxiaType,"_hypoxia.html",sep = ""))
	}

	# saveWidget(dygraph(hypoxiaExtent[,c("less0_idw","less0Mean_basis","less0Upper_basis","less0Lower_basis")]/n, main = sprintf("Total Area Less than 0.01",area)) %>%  dySeries(c("less0Lower_basis", "less0Mean_basis", "less0Upper_basis"), label = "Basis Interpolation") %>% dyRangeSelector(height = 20))
	

	# print(dygraph(hypoxiaExtent[,c("less2_idw","less2Mean_basis","less2Upper_basis","less2Lower_basis")]/n, main = "Less than 2") %>%  dySeries(c("less2Lower_basis", "less2Mean_basis", "less2Upper_basis"), label = "Basis Interpolation") %>% dyRangeSelector(height = 20))
	

	# print(dygraph(hypoxiaExtent[,c("less4_idw","less4Mean_basis","less4Upper_basis","less4Lower_basis")]/n, main = "Less than 4") %>%  dySeries(c("less4Lower_basis", "less4Mean_basis", "less4Upper_basis"), label = "Basis Interpolation") %>% dyRangeSelector(height = 20))
}

main <- function(){
	for(timeAggType in c("daily")){
		for(r in c(10:10)){
		 print(system.time(analysis(2014,timeAggType,r)))
		# analysis(2015,timeAggType,r)
		#	analysis(2016, timeAggType, r)
		}
	}
	# summary_plot(2014, "daily",10)
}

metaFolder <- NULL
outputFolder <- NULL
main()

# 
# library(htmlwidgets)
# for(r in c(5:12)){
#   print(r)
#   summary_plot(2014,"daily",r)
#   summary_plot(2014,"hourly",r)
#   
#   summary_plot(2015,"daily",r)
#   summary_plot(2015,"hourly",r)
# }
# debug(summary_plot)

# # read prediction


# timeIndex <- index(data)



# library(ggmap)

# # t <- 2331
# t <- 999
# sim <- 100
# threshold <- 4
# loggerNames  <- as.numeric(colnames(data))
# subData <- data.frame(logger = loggerNames,DO = as.numeric(data[t,]))
# subData <- merge(subData, loggerInfo,by.x = "logger",by.y = "loggerID") %>% rename(value = DO)

# qmplot(longitude,latitude, data = subData, color = ifelse(value<threshold,0,1), size = I(5)) + scale_color_gradient(low = "red",high = "blue", limits = c(0,1),name = "DO")
# qmplot(longitude,latitude,data = grid, color = ifelse(finalPrediction_basis[sim,t,]<threshold,0,1)) + scale_color_gradient(low = "red",high = "blue", limits = c(0,1),name = "DO")
# qmplot(longitude,latitude,data = grid, color = ifelse(finalPrediction_idw[1,t,]<threshold,0,1)) + scale_color_gradient(low = "red",high = "blue", limits = c(0,1),name = "DO")



# dygraph(hypoxiaExtent[,c("less0_idw","less0Mean_basis","less0Upper_basis","less0Lower_basis")]/n*area, main = "Less than 0.01") %>%  dySeries(c("less0Lower_basis", "less0Mean_basis", "less0Upper_basis"), label = "Basis Interpolation") %>% dyRangeSelector(height = 20)
# dygraph(hypoxiaExtent[,c("less2_idw","less2Mean_basis","less2Upper_basis","less2Lower_basis")]/n*area, main = "Less than 2") %>%  dySeries(c("less2Lower_basis", "less2Mean_basis", "less2Upper_basis"), label = "Basis Interpolation") %>% dyRangeSelector(height = 20)
# dygraph(hypoxiaExtent[,c("less4_idw","less4Mean_basis","less4Upper_basis","less4Lower_basis")]/n*area, main = "Less than 4") %>%  dySeries(c("less4Lower_basis", "less4Mean_basis", "less4Upper_basis"), label = "Basis Interpolation") %>% dyRangeSelector(height = 20)




# method : loglik, baye, IDW
# plot(basis_hypoxiaExtent)



# allPlot <- list()
# for(logger in upperLogger$loggerID){
#     p <- plot_DO_temp(retriveLoggerData_DO_Temp(c(logger),year,"hourly","raw",NULL))
#     allPlot[[paste("logger_",logger,sep="")]] <- p
# }
# 
# 
# saveRDS(allPlot,"../output/allPlot_dygraphs.rds")

# loggerData_DO <- retriveLoggerData(upperLogger$loggerID,year,"DO","hourly","AVG",c("2014-06-15","2014-10-06"))
# plot_value(loggerData,"dygrphs")



# }

