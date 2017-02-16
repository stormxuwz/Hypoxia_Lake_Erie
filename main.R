rm(list = ls())
setwd("/Users/wenzhaoxu/Developer/Hypoxia/Hypoxia_Lake_Erie")
source("./src/database.R")
source("config.R")
source("./src/plot.R")
source("./src/interpolation.R")
source("./src/crossValidation.R")

dataTransformation <- function(df, method = "sqrt"){
	# data transformation
	# Args:
	# 	df: a zoo data frame
	time <- index(df)
	df <- as.matrix(df)

	for(i in 1:ncol(df)){
		df[,i] <- ifelse(df[,i]<-0.1, 0.1, df[,i])
		if(method=="sqrt"){
			df[,i] <- sqrt(df[,i])
		}
	}

	df <- zoo(df, order.by = time)
	return(df)
}


analysis <- function(year,timeAggType, method, withCV = FALSE, ...){
	# function to do interpolation
	# =============================
	# the logger info
	loggerInfo <- retriveGeoData(year,"B") %>% arrange(loggerID)
	
	# the data used to interpolation
	#timeRange = c("2014-07-01","2014-09-01")
	timeRange = NULL # no specific time range
	data <- retriveLoggerData(loggerInfo$loggerID,year,"DO",timeAggType,"AVG",timeRange = timeRange,transform = TRUE) %>% na.omit()  # remove the data
	
	time <- index(data)
	
	# Create the results folder
	metaFolder <<- sprintf("../meta_%d_%s_%s/", year,timeAggType,method)
	outputFolder <<- sprintf("../output_%d_%s/", year,timeAggType,method)
	
	if(method == "basis"){
		metaFolder <<- sprintf("../meta_%d_%s_%s_%s_%d/", 
			year,timeAggType,method,list(...)$fitMethod,list(...)$r)
		outputFolder <<- sprintf("../output_%d_%s_%s_%s_%d/", 
			year,timeAggType,method,list(...)$fitMethod, list(...)$r)
	}else{
		metaFolder <<- sprintf("../meta_%d_%s_%s/", year,timeAggType,method)
		outputFolder <<- sprintf("../output_%d_%s/", year,timeAggType,method)
	}
	dir.create(file.path("..", substring(metaFolder,4)), showWarnings = FALSE)
	dir.create(file.path("..", substring(outputFolder,4)), showWarnings = FALSE)
	
	if(withCV)
		crossValidation(data,loggerInfo, "svd",rList=c(r),method = "basis")
	
	grid <- createGrid(loggerInfo,by.x = mapDx, by.y = mapDy)

	hypoxiaExtent <- interpolation_main(
		data = data,
		locationInfo = loggerInfo,
		method = method, 
		grid = grid, ...)
	
	saveRDS(list(res=hypoxiaExtent,timeIndex = time, grid = grid),
		paste0(outputFolder,"hypoxiaExtent.rds"))

	# read results and plot
	summary_plot(year,timeAggType, method,...)

}

main <- function(year){
	for(timeAggType in c("daily","hourly")){
		print(system.time(analysis(year,timeAggType,"IDW")))
		for(r in c(5,10,15)){
		 print(system.time(analysis(year,timeAggType,"basis",
		 	basisDecomp = "svd",
		 	simNum = 100,
		 	fitMethod = "loglik",
		 	r = r,
		 	residualMethod = "IDW")))

		 print(system.time(analysis(year,timeAggType,"basis",
		 	basisDecomp = "svd",
		 	simNum = 100,
		 	fitMethod = "baye",
		 	r = r,
		 	residualMethod = "IDW")))
		}
	}
}

metaFolder <- NULL
outputFolder <- NULL
main(2014)

