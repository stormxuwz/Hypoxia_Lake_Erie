# This is to estimate the changes of thermocline


series <- retriveLoggerData(10384439,year,"Temp","hourly","raw",c("2014-07-18","2014-10-06"))


onsetDetection <- function(series,window = 24){
	diff_series <- diff(c(series))
	z <- c(rep(NA,window+1))
	for(i in (window+1):length(diff_series) ){
		subDiff <- diff_series[(i-window):i]
		subDiff <- subDiff*(subDiff>0)
		z <- c(z,sum(subDiff))
	} 
	return(z)
}


z <- onsetDetection(series)

series$sumGradient <- z
dygraph(series) %>% dySeries("sumGradient",axis="y2")