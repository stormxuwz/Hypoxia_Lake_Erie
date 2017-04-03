

plot_hypoxia <- function(){

	# trendPrediction <- readRDS(paste(metaFolder,"trend_prediction.rds",sep = ""))
	# coef_df <- readRDS(paste(metaFolder,"coef_df.rds",sep = ""))
	# DO_res <- readRDS(paste(metaFolder,"DO_res.rds",sep = ""))
	# basis <- readRDS(paste(metaFolder,"basis.rds",sep = "")) %>% zoo(order.by = timeIndex)


	
	# plot interpolation range
	
	# calculate the spatio-temporal variogram
	# stvgm <- st_variogram(DO_res,loggerInfo)
	# print(plot(stvgm$vgmModel, map=FALSE))

	hypoxiaExtent <- cbind(data.frame(basis_hypoxiaExtent), data.frame(IDW_hypoxiaExtent)) %>% zoo(order.by = timeIndex)

	for(hypoxiaType in c("less0","less2","less4")){
		p <- dygraph(hypoxiaExtent[,c(krigNames,IDWNames)]/n*area, main = sprintf("Total Area %.2f km^2",area)) %>%
		dyAxis("y", label = "Hypoxia Area (km^2)",valueRange = c(0, area+10)) %>%
		dySeries(krigNames, label = "Spatiao-Temporal Kriging") %>% 
		dySeries(IDWNames, label = "IDW Interpolation") %>% 
		dyRangeSelector(height = 20)
		saveWidget(p, file = paste(outputFolder,hypoxiaType,"_hypoxia.html",sep = ""))
	}
}



RMSE <- function(x,y){
	return(round(mean((x-y)^2),digits = 3))
}

nash_coef <- function(pred,obs){
	return(round(1-sum((obs-pred)^2)/sum( (obs-mean(obs))^2 ),2))
}

comparisonCV <- function(year,ID,timeInterval){
	res = list()
	
	# store all the results
	for(r in c(5,10,15)){
		savedFile = sprintf("../meta_%d_%s_%d/cv_r%d_%d.rds",year,timeInterval,r,r,ID)
		d <- readRDS(savedFile)
		info <- strsplit(strsplit(basename(savedFile),"[.]")[[1]][1],"_")[[1]]
		nsim = length(d$trend)
		
		for(i in 1:nsim){
			d$trend[,i] <- d$trend[,i]+d$pred_res
			d$trend[,i] <- d$trend[,i]*(d$trend[,i]>0)
		}

		m <- apply(d$trend,1,quantile,probs = c(0.5))
		upper = apply(d$trend, 1, quantile,probs = c(0.95))
		lower = apply(d$trend, 1, quantile,probs = c(0.05))
		
		plotDF <- data.frame(m = m, upper = upper, lower = lower,true = d$label,time = 1:length(m))
		res[[as.character(r)]] <- plotDF
	}
	
	yRange <- range(m,upper,lower,d$label)
	
	# plot
	library(gridExtra)
	
	p5 = ggplot(data = res[["5"]]) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,m),color = "black") + 
		geom_line(aes(time,true),color = "red",alpha = 0.8) + 
		ylim(yRange)+ylab("DO")+
		ggtitle(paste(ID,5,timeInterval,RMSE(res[["5"]]$m,res[["5"]]$true),sep = "_"))
	
	p10 = ggplot(data = res[["10"]]) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,m),color = "black") + 
		geom_line(aes(time,true),color = "red",alpha = 0.8) + 
		ylim(yRange)+ylab("DO")+
		ggtitle(paste(ID,10,timeInterval,RMSE(res[["10"]]$m,res[["10"]]$true),sep = "_"))
	
	p15 = ggplot(data = res[["15"]]) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,m),color = "black") + 
		geom_line(aes(time,true),color = "red",alpha = 0.8) + 
		ylim(yRange)+ylab("DO")+
		ggtitle(paste(ID,15,timeInterval,RMSE(res[["15"]]$m,res[["15"]]$true),sep = "_"))
	
	png(paste("~/Desktop/CV",year,ID,timeInterval,".png",sep = "_"),width = 1500)	
	print(grid.arrange(p5, p10, p15, ncol=3))
	dev.off()
}






year = 2014
loggerInfo <- retriveGeoData(year,"B") %>% arrange(loggerID)

for(loggerID in loggerInfo$loggerID){
	comparisonCV(year,loggerID,"hourly")
}



year = 2015
loggerInfo <- retriveGeoData(year,"B") %>% arrange(loggerID)

for(loggerID in loggerInfo$loggerID){
	comparisonCV(year,loggerID,"hourly")
}






