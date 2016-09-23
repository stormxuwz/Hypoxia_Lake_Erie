# this script is used to detect outlier
# Including (1) a simple moving average


tsDecomposition <- function(x,season){
	#stl functions to do time series decomposition
	res <- ts(x,frequency = season) %>% stl(s.window = 7, robust = TRUE) # s.window = 7 follows the suggestions
	print(plot(res))
	return(res$time.series)
	# time.series has three columns: seasonal, trend and remainder
}

remainderOutlier <- function(x,threshold,window= 3*24){
	ut <- function(y) {m = median(y); median(y) + threshold * median(abs(y - m))}
	lt <- function(y) {m = median(y); median(y) - threshold * median(abs(y - m))}
	
	z_ut <- rollapply(zoo(x), window, ut, align="right")
	z_lt <- rollapply(zoo(x), window, lt, align="right")
	
	z_ut <- c(rep(z_ut[1], window-1), z_ut)
	z_lt <- c(rep(z_lt[1], window-1), z_lt)
	
	#print(z_lt)
	
	outliers <- (x > z_ut) | (x < z_lt)
	
	return(outliers)
}


stlOutlierDetection <- function(x, threshold,testSeasons = c(8,17,24)){
	remainder <- c()
	outlier <- c()
	for(season in testSeasons){
		decomp <- tsDecomposition(x,season)
		remainder <- c(remainder,decomp[,3])
		outlier <- c(outlier,remainderOutlier(decomp[,3],threshold))
	}
	remainder <- matrix(remainder, ncol = 3)
	#print(head(outlier))
	outlier <- matrix(outlier,ncol = 3)
	outlier <- rowSums(outlier)
	# print(outlier)
	return(outlier)
}


outlierTest <- function(){
	library(tsoutliers)
	year <- 2014
	loggerInfo <- retriveGeoData(year,"B")
	data <- retriveLoggerData(loggerInfo$loggerID,year,"DO","hourly","AVG",transform = TRUE) %>% na.omit()
	series <- data[,10]
	series2 <- series
	ol <- stlOutlierDetection(series,3)
	
	#series2[ol<3] <- NA # 3 means all three tests are and
	#plot(series)
	#points(series2,col = "red")
	#return(ol)
	#dat.ts <- ts(series, frequency = 8)
	#a <- stl(x = dat.ts, s.window = 7,robust = TRUE)
	#plot(a)
	#plot(dat.ts)
	#lines(a$time.series[,1]+a$time.series[,2],col = "red")
	
	
	#data.ts.outliers <- tso(dat.ts)
	
	#stl(x = data.ts.outliers)
	
	#plot(data.ts.outliers)
}

#ol <- outlierTest()
