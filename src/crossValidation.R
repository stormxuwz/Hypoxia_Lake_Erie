
crossValidation <- function(data, locationInfo, method, ...){
	
	crossPred <- function(j,...){
		# j is the left out logger index
		trainData <- data[,-j]
		testData <- data[,j]
		
		trainLocation <- locationInfo[-j,]
		targetLocation <- locationInfo[j,c("longitude","latitude","x","y","bathymetry")]
		targetLocation$convexIndex <- 1
		
		if(method == "basis"){
			res <- basis_interpolation_step1(
				data = trainData, 
				locationInfo = trainLocation,
				basisDecomp = list(...)$basisDecomp,
				simNum = list(...)$simNum, 
				fitMethod = list(...)$fitMethod,
				r = list(...)$r, 
				residualMethod = list(...)$residualMethod,
				grid = targetLocation,
				saveMeta = FALSE)

			pred <- basis_interpolation_step2(res[[1]],res[[2]],nSim = 1000,parallel = TRUE,returnHypoxia = FALSE, saveMeta = FALSE)
			nsim <- length(pred) # total simulations

			# extract the first column which should only have the one column
			stopifnot(ncol(pred[[1]])==1)

			pred  <- data.frame(lapply(pred,function(x) return(x[,1]))) 
			colnames(pred) <- paste0("sim_",as.character(1:nsim))
			pred_res = res[[2]][[1]]
		}
		else if(method == "IDW"){
			pred <- idw_interpolation_main(
				data = trainData, 
				locationInfo = trainLocation, 
				grid = targetLocation)[[1]]
			pred <- data.frame(sim_0 = pred)
			pred_res = rep(0,nrow(trainData))
		}
		
		return(list(trend = pred, pred_res = pred_res,label = as.numeric(testData)))
	}

	for(j in 1:nrow(locationInfo)){
		cp <- crossPred(j,...)
		fileName <- sprintf("%scv_%s_%s.rds",outputFolder,locationInfo[j,"loggerID"],method)
		saveRDS(cp,fileName)
		crossValidation_summary(fileName,withResidual = TRUE)
	}
}


crossValidation_summary <- function(savedFile,withResidual= TRUE){
	d <- readRDS(savedFile)
	info <- strsplit(strsplit(basename(savedFile),"[.]")[[1]][1],"_")[[1]]
	nsim  <- length(d$trend)
	
	if(withResidual){
		for(i in 1:nsim){
			d$trend[,i] <- d$trend[,i]+d$pred_res
			d$trend[,i] <- d$trend[,i]*(d$trend[,i]>0)
		}
		#m = rowMeans(d$trend)
		m <- apply(d$trend,1,quantile,probs = c(0.5))
		upper = apply(d$trend, 1, quantile,probs = c(0.95))
		lower = apply(d$trend, 1, quantile,probs = c(0.05))
	}
	# ifelse(m-std<0,0,m-std)
	plotDF <- data.frame(m = m, upper = upper, lower = lower,true = d$label,time = 1:length(m))
	
	yRange <- range(m,upper,lower,d$label)

	p <- ggplot(data = plotDF) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,m),color = "black") + 
		geom_line(aes(time,true),color = "red",alpha = 0.8) + 
		ylim(c(-1,15))+ylab("DO")+
		ggtitle(paste0("RMSE: ",RMSE(plotDF$m,plotDF$true)))
	
	pdf(sprintf("%splot_%s.pdf",outputFolder,paste(info,collapse = "_")),width = 6, height = 4)
	print(p)
	dev.off()
}

RMSE <- function(x,y){
	return(round(mean((x-y)^2),digits = 3))
}

nash_coef <- function(pred,obs){
	return(round(1-sum((obs-pred)^2)/sum( (obs-mean(obs))^2 ),2))
}


library(dplyr)
for(year in 2014:2015){
	loggerInfo <- retriveGeoData(year,"B") %>% arrange(loggerID)
	for(ID in loggerInfo$loggerID){
		for(fitMethod in c("loglik","baye"))
		{
			compareCV_IDW(year, ID, "hourly", 15, fitMethod)
		}
	}
}



compareCV_IDW <- function(year, ID, timeAggType,r,fitMethod){

	# timeRange = NULL # no specific time range
	# data <- retriveLoggerData(loggerInfo$loggerID,year,"DO",timeAggType,"AVG",timeRange = timeRange,transform = TRUE) %>% na.omit()  # remove the data
	# time <- index(data)
	basis_file <- sprintf("../output_CV_%d_%s_basis_%s_%d/cv_%d_basis.rds",year,timeAggType,fitMethod,r,ID)
	IDW_file <- sprintf("../output_CV_%d_%s_IDW/cv_%d_IDW.rds",year,timeAggType,ID)

	d <- readRDS(basis_file)
	nsim  <- length(d$trend)
	for(i in 1:nsim){
		d$trend[,i] <- d$trend[,i]+d$pred_res
		d$trend[,i] <- d$trend[,i]*(d$trend[,i]>0)
	}
	m <- apply(d$trend,1,quantile,probs = c(0.5))
	upper = apply(d$trend, 1, quantile,probs = c(0.95))
	lower = apply(d$trend, 1, quantile,probs = c(0.05))
	idw <- readRDS(IDW_file)$trend[,1]
	
	plotDF <- data.frame(m = m, upper = upper, lower = lower,true = d$label,idw = idw, time = 1:length(m))

	p <- ggplot(data = plotDF) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,m),color = "black") + 
		geom_line(aes(time,true),color = "red",alpha = 0.8) + 
		geom_line(aes(time,idw),color = "yellow",alpha = 0.8,linetype = 2) + 
		ylim(c(-1,17))+ylab("DO")+
		ggtitle(paste0("RMSE basis/IDW: ",RMSE(plotDF$m,plotDF$true),"/",RMSE(plotDF$idw,plotDF$true)))
	
	png(paste("../finalResults/",year,ID,timeAggType,r,fitMethod,".png",sep = "_"),width = 1500)	
	print(p)
	dev.off()
}



comparisonCV_basis <- function(year,ID,timeInterval){
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


