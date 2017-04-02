rm(list = ls())
setwd("/Users/wenzhaoxu/Developer/Hypoxia/Hypoxia_Lake_Erie")
source("./src/database.R")
source("config.R")
source("./src/plot.R")
source("./src/interpolation.R")
source("./src/crossValidation.R")

dataTransformation <- function(df, method = "sqrt"){
	# Do data transformation if necessary
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


analysis <- function(year,timeAggType, method, withCV = TRUE, ...){
	# function to do interpolation
	# =============================
	# the logger info
	loggerInfo <- retriveGeoData(year,"B") %>% arrange(loggerID)
	
	# the data used to interpolation
	#timeRange = c("2014-07-01","2014-09-01")
	timeRange = NULL # no specific time range
	data <- retriveLoggerData(loggerInfo$loggerID,year,"DO",timeAggType,"AVG",timeRange = timeRange,transform = TRUE) %>% na.omit()  # remove the data
	time <- index(data)
	
	#### Create the results folder ####
	prefix_meta <- "meta"
	prefix_output <- "output"
	
	if(withCV){
		 prefix_meta <-  paste(prefix_meta,"CV",sep = "_")
		 prefix_output <-  paste(prefix_output,"CV",sep = "_")
	}

	if(method == "basis"){
		metaFolder <<- sprintf("../%s_%d_%s_%s_%s_%d/", prefix_meta,
			year,timeAggType,method,list(...)$fitMethod,list(...)$r)
		outputFolder <<- sprintf("../%s_%d_%s_%s_%s_%d/", prefix_output,
			year,timeAggType,method,list(...)$fitMethod, list(...)$r)
	}else{
		metaFolder <<- sprintf("../%s_%d_%s_%s/", prefix_meta,year,timeAggType,method)
		outputFolder <<- sprintf("../%s_%d_%s_%s/", prefix_output, year,timeAggType,method)
	}
	dir.create(file.path("..", substring(metaFolder,4)), showWarnings = FALSE)
	dir.create(file.path("..", substring(outputFolder,4)), showWarnings = FALSE)
	#### Done ####


	if(withCV){
		# Do cross validation
		crossValidation(data,loggerInfo,method = method,...)
	}else{
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
main(2015)
