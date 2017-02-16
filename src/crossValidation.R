
crossValidation <- function(data, locationInfo, basis, rList = c(5,10,15),method = "IDW"){
	
	crossPred <- function(r,j){
		# r is the number of basis, j is the left out logger index
		trainData <- data[,-j]
		testData <- data[,j]
		
		trainLocation <- locationInfo[-j,]
		targetLocation <- locationInfo[j,c("longitude","latitude","x","y","bathymetry")]
		targetLocation$convexIndex <- 1
		
		if(method == "basis"){
			res <- basis_interpolation_step1(
				data = trainData, 
				locationInfo = trainLocation,
				basis = basis,
				simNum = 100, 
				intMethod = "loglik", 
				r = r, 
				residualMethod = "IDW",
				saveMeta = FALSE, 
				grid = targetLocation)

			pred <- basis_interpolation_step2(res[[1]],res[[2]],nSim = 1000,TRUE,FALSE)
			nsim <- length(pred) # total simulations

			# extract the first column which should only have the one column
			pred  <- data.frame(lapply(pred,function(x) return(x[,1]))) 
			colnames(pred) <- paste0("sim_",as.character(1:nsim))
			pred_res = res[[2]][,,1]
		}
		else if(method == "IDW"){
			pred <- idw_interpolation_main(
				data = trainData, 
				locationInfo = trainLocation, 
				grid = targetLocation)[1,,]
			pred <- data.frame(sim_0 = pred)
			pred_res = rep(0,nrow(trainData))
		}
		
		return(list(trend = pred, pred_res = pred_res,label = as.numeric(testData)))
	}

	for(r in rList){
		for(j in 1:nrow(locationInfo)){
			cp <- crossPred(r,j)
			fileName <- sprintf("%scv_r%d_%s_%s.rds",metaFolder,r,locationInfo[j,"loggerID"],method)
			saveRDS(cp,fileName)
			crossValidation_summary(fileName,withResidual = TRUE)
		}
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

	p = ggplot(data = plotDF) + 
		geom_ribbon(aes(time,ymin = lower, ymax = upper), fill = "grey70") +
		geom_line(aes(time,m),color = "black") + 
		geom_line(aes(time,true),color = "red",alpha = 0.8) + 
		ylim(c(-1,15))+ylab("DO")+
		ggtitle(paste0("RMSE: ",RMSE(plotDF$m,plotDF$true)))
	
	pdf(sprintf("%splot_%s.pdf",metaFolder,paste(info,collapse = "_")),width = 6, height = 4)
	print(p)
	dev.off()
}

RMSE <- function(x,y){
	return(round(mean((x-y)^2),digits = 3))
}

nash_coef <- function(pred,obs){
	return(round(1-sum((obs-pred)^2)/sum( (obs-mean(obs))^2 ),2))
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


