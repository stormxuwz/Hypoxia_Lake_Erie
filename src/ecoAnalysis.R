
require(chron)

buoyDataFolder <- "../input/"
EPASite <- c("ER30","ER31","ER32","ER36","ER37","ER38","ER42","ER43","ER73","ER78")

#buoyName <- "45005"
#year <- 2014
#dtype <- "h"

#V1. V2.V3 V4 V5 V6.  V7.   V8.  V9.   V10.  V11 V12.   V13.  V14.  V15. V16.  V17.  V18
#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE

# evernote:///view/1223127095/s45/1819a5a5-3513-4ce4-9731-eb8927322eb7/1819a5a5-3513-4ce4-9731-eb8927322eb7/

# WDIR
# Wind direction (the direction the wind is coming from in degreesclockwise from true N) during the same period used for WSPD. See Wind Averaging Methods
# WSPD
# Wind speed (m/s) averaged over an eight-minute period for buoys anda two-minute period for land stations.  Reported Hourly.  See Wind Averaging Methods.
# GST
# Peak 5 or 8 second gust speed (m/s) measured during the eight-minuteor two-minute period.  The 5 or 8 second period can be determined by payload, See the Sensor Reporting, Sampling, andAccuracy section.
# WVHT
# Significant wave height (meters) is calculated as the average of thehighest one-third of all of the wave heights during the 20-minutesampling period.  See the Wave Measurementssection.
# DPD
# Dominant wave period (seconds) is the period with the maximum waveenergy. See the Wave Measurements section.
# APD
# Average wave period (seconds) of all waves during the 20-minuteperiod. See the Wave Measurements section.
# MWD
# The direction from which the waves at the dominant period (DPD) are coming. The units are degrees from true North, increasing clockwise, with North as 0 (zero) degrees and East as 90 degrees.  See the Wave Measurements section.
# PRES
# Sea level pressure (hPa).  For C-MAN sites and Great Lakes buoys,the recorded pressure is reduced to sea level using the methoddescribed in NWS Technical Procedures Bulletin 291 (11/14/80). (labeled BAR in Historical files)
# ATMP
# Air temperature (Celsius).  For sensor heights on buoys, seeHull Descriptions.  For sensor heights atC-MAN stations, see C-MAN Sensor Locations
# WTMP
# Sea surface temperature (Celsius). For buoys the depth is referencedto the hull's waterline. For fixed platforms it varies with tide, but isreferenced to, or near Mean Lower Low Water (MLLW).
# DEWP
# Dewpoint temperature taken at the same height as the airtemperature measurement.
# VIS
# Station visibility (nautical miles). Note that buoy stations are limited to reports from 0 to 1.6 nmi.
# PTDY
# Pressure Tendency is the direction (plus or minus) and the amount ofpressure change (hPa)for a three hour period ending at the time ofobservation. (not in Historical files)
# TIDE
# The water level in feet above or below Mean Lower Low Water (MLLW).


readBuoyData <- function(year, buoyName, dtype = "h"){
	filename <- paste0(buoyDataFolder,"/otherData/", buoyName,dtype,year,".txt")
	data <- read.table(filename,skip = 1,sep = "") 

	data$time <- paste0(data[,1],"-",data[,2],"-",data[,3]," ",data[,4]) 

	if(dtype == "h"){
		data  <- data %>% 
			rename(WDIR=V6,WTMP = V15,ATMP =V14,WVHT =V9, WSPD = V7) %>% 
			group_by(time) %>% 
			summarise_all(funs(mean)) %>% 
			data.frame()

		data$WDIR = ifelse(data$WDIR == 999, NA, data$WDIR)
		data$WVHT = ifelse(data$WVHT >= 40, NA, data$WVHT)
		data$WTMP = ifelse(data$WTMP >= 150, NA, data$WTMP)
		data$ATMP = ifelse(data$ATMP > 150, NA, data$ATMP)
		data$WSPD = ifelse(data$WSPD == 99, NA, data$WSPD)

		data$southeastShoreWSPD <- (cos(data$WDIR - 315)/180*pi)*data$WSPD
		data$southwestShoreWSPD <- (sin(data$WDIR - 315)/180*pi)*data$WSPD

		time <- data$time %>% as.POSIXlt(format = "%Y-%m-%d %H",tz = "GMT")
		data <- zoo(data[,c("WDIR","WTMP","ATMP","WVHT","WSPD","southeastShoreWSPD","southwestShoreWSPD")],order.by = time)

	}else{
		data <- data %>% rename(WDIR=V6,WSPD = V7,GDR =V8,GST =V9) %>% 
			group_by(time) %>% 
			summarise_all(funs(mean)) %>% 
			data.frame()

		data$WDIR = ifelse(data$WDIR == 999, NA, data$WDIR)
		data$GST = ifelse(data$GST == 99, NA, data$GST)
		data$GDR = ifelse(data$GDR == 999, NA, data$GDR)
		data$WSPD = ifelse(data$WSPD == 99, NA, data$WSPD)

		data$southeastShoreWSPD <- (cos(data$WDIR - 315)/180*pi)*data$WSPD
		data$southwestShoreWSPD <- (sin(data$WDIR - 315)/180*pi)*data$WSPD

		time <- data$time %>% as.POSIXlt(format = "%Y-%m-%d %H",tz = "GMT")
		data <- zoo(data[,c("WDIR","WSPD","GDR","GST","southeastShoreWSPD","southwestShoreWSPD")],order.by = time)
	}

	data <- subset(data, time > paste(year, "06-15",sep = "-") & time < paste(year, "10-15",sep = "-"))
	names_buoyData <- names(data)
	names(data) <- paste(names_buoyData,buoyName,dtype,sep="_")

	return(data)
}


convertWindDirection <- function(WDIR, side = "south"){
	if(side == "south"){
		degree_projection = rep(0, 360)

		for(d in 225:315){
			degree_projection[d + 1] = 1 - abs((315-d)/90)
		}
		
		for(d in 315:359){

			degree_projection[d + 1] = 1 - abs((315 - d)/90)
		}

		for(d in 0:45){
			degree_projection[d + 1] = 0.5 - abs((0 - d)/90)
		}

		newWDIR <- sapply(WDIR, function(i){degree_projection[i+1]})
	}


	return(newWDIR)
}


read45132 <- function(year){
	filename <- paste0(buoyDataFolder,"/otherData/", "c45132.csv")
	data <- read.csv(filename)
	time <- as.POSIXlt(data$DATE, format = "%m/%d/%Y %H:%M",tz = "GMT") %>% 
		format("%m/%d/%Y %H") %>% 
		as.POSIXlt(format = "%m/%d/%Y %H",tz = "GMT")
	
	data <- data %>% 
		select(WDIR, ATMS, VCAR, WSPD, SSTP,VWH.) %>% 
		rename(WVHT = VCAR, WTMP = SSTP, WVHT2 = VWH.) # VWH: Characteristic significant wave height (reported by the buoy) (m)


	data <- zoo(data, order.by = time)

	data <- subset(data, time > paste(year, "06-15",sep = "-") & time < paste(year, "10-15",sep = "-")) 

	tmp <- names(data)
	names(data) <- paste(tmp,"c45132",sep="_")

	return(data)

}

scale_0_1 <- function(x){
	min_x = min(x, na.rm = T)
	max_x = max(x, na.rm = T)

	return( (x-min_x)/(max_x - min_x))
}

readAllData <- function(year){
	erieDO <- getLakeDO(year, "B", "hourly") %>% na.omit()
	loggerInfo <- erieDO$loggerInfo
	
	buoyData_45005_h <- readBuoyData( year, "45005","h")
	buoyData_45164_h <- readBuoyData( year, "45164","h")
	buoyData_cndo1_h <- readBuoyData( year, "cndo1","h")
	buoyData_faio1_h <- readBuoyData(year, "faio1","h")	# only wind data, temperature data available 
	buoyData_gelo1_h <- readBuoyData(year, "gelo1","h")	# no data in sampling time in 2015
	buoyData_cblo1_h <- readBuoyData(year, "cblo1","h") # 2015 data, no southeastShoreWSPD_cblo1_h in sampling time

	buoyData_45169_h <- readBuoyData(year, "45169","h")
	buoyData_45132_h <- read45132(year)


	samplingData <- erieDO$samplingData %>% na.omit()
	t <- index(samplingData)

	basis_r <- ifelse(year == 2015, 5, 3)

	nmfDecomp <- readRDS(sprintf("%s/results/cluster/nmfBasis_%d_%s_%d.rds",outputBaseName, year,"hourly",basis_r))
	nmfDecomp <- nmfDecomp$basis %>% data.frame()
	# nmfDecomp <- NMF_basis(samplingData,5)$basis %>% data.frame()

	names(nmfDecomp) <- paste("NMF",1:ncol(nmfDecomp),sep = "_")
	
	nmfDecomp <- zoo(nmfDecomp, order.by = t)
	
	# data <- zoo(data[,c("WDIR","WTMP","ATMP","WVHT","WSPD",southeastShoreWSPD","southwestShoreWSPD")],order.by = time)
	
	if(year == 2015){
		mergedData <- merge(nmfDecomp, buoyData_45169_h)[,c("NMF_5","WVHT_45169_h")]
		newt <- index(mergedData) %>% as.POSIXct()
		mergedData <- mergedData %>%
			data.frame() %>% 
			rename(Basis_5 = NMF_5) %>% 
			mutate(time = newt) %>% 
			subset(newt >= t[1] & newt <= t[length(t)]) %>% 
			mutate(Basis_5 = scale_0_1(Basis_5),WVHT_45169_h = scale_0_1(WVHT_45169_h)) %>%
			melt(id.vars = c("time"))

		pdf(sprintf("%s/results/%d_nmfBasis_%d_WVHT_%s.pdf",outputBaseName, year, basis_r, "45169h"), width = 6.25, height = 3)
		print(ggplot(mergedData) + geom_line(aes(time, value, color = variable), alpha = 0.85) + 
			theme_bw() + ylab("Scaled value") + xlab("") + theme(legend.position="top"))
		dev.off()
			# geom_line(aes(newt, scale_0_1(WVHT_45169_h)),color = "blue", alpha = 0.5) 


		mergedData <- merge(nmfDecomp, buoyData_45132_h)[,c("NMF_4","WVHT_c45132")]
		newt <- index(mergedData) %>% as.POSIXct()
		mergedData <- mergedData %>% 
			data.frame() %>% 
			rename(Basis_4 = NMF_4) %>% 
			mutate(time = newt) %>% 
			subset(newt >= t[1] & newt <= t[length(t)]) %>% 
			mutate(Basis_4 = scale_0_1(Basis_4),WVHT_c45132 = scale_0_1(WVHT_c45132)) %>%
			melt(id.vars = c("time"))

		pdf(sprintf("%s/results/%d_nmfBasis_%d_WVHT_%s.pdf",outputBaseName, year, basis_r, "c45132"), width = 6.25, height = 3)
		print(ggplot(mergedData) + geom_line(aes(time, value, color = variable), alpha = 0.85) + 
			theme_bw() + ylab("Scaled value") + xlab("") + theme(legend.position="top"))
		dev.off()
	}else{
		mergedData <- merge(nmfDecomp, buoyData_45164_h)[,c("NMF_2","WVHT_45164_h")]
		newt <- index(mergedData) %>% as.POSIXct()
		mergedData <- mergedData %>%
			data.frame() %>% 
			rename(Basis_2 = NMF_2) %>% 
			mutate(time = newt) %>% 
			subset(newt >= t[1] & newt <= t[length(t)]) %>% 
			mutate(Basis_2 = scale_0_1(Basis_2),WVHT_45164_h = scale_0_1(WVHT_45164_h)) %>%
			melt(id.vars = c("time"))

		pdf(sprintf("%s/results/%d_nmfBasis_%d_WVHT_%s.pdf",outputBaseName, year, basis_r, "45164h"), width = 6.25, height = 3)
		ggplot(mergedData) + geom_line(aes(time, value, color = variable), alpha = 0.85) + 
			theme_bw() + ylab("Scaled value") + xlab("") + theme(legend.position="top")
		dev.off()
	}



	fullData <- merge(nmfDecomp, buoyData_45005_h) %>% 
				merge(buoyData_45164_h) %>% merge(buoyData_cndo1_h) %>% merge(buoyData_faio1_h)

	merged_noApprox <- merge(nmfDecomp,buoyData_45005_h) %>% na.omit()

	merged <- merge(nmfDecomp,buoyData_45005_h, all = TRUE)  %>% na.approx() %>% na.omit()
	spectrum(merged$ATMP_45005_h)
	spectrum(merged$NMF_3)

	plot(merge(nmfDecomp,buoyData_45132_h)[,c("NMF_4","WVHT_c45132","WVHT2_c45132")] %>% na.omit())

	plot(merge(nmfDecomp,buoyData_45132_h)[,c("NMF_2","WVHT_45169_h")] %>% na.omit())
	# plot(fullData[,c("NMF_1","WDIR_45005_h")])

	ggplot(data = data.frame(fullData))+geom_line(aes(t,scale_0_1(NMF_3))) + geom_line(aes(t, newWDIR),color = "red")
}


# function to analyze the EPA loggers

filterEPA <- function(lakeDO, changeColumnName = TRUE){
	lakeDO$loggerInfo <- subset(lakeDO$loggerInfo, site %in% EPASite) %>%
				arrange(site)
	lakeDO$samplingData <- lakeDO$samplingData[,lakeDO$loggerInfo$loggerID]
	
	if(changeColumnName) names(lakeDO$samplingData) <- lakeDO$loggerInfo$site
	return(lakeDO)
}


getSiteSeries <- function(allData, site){
	d_2014 <- allData[["data2014"]]$samplingData[-c(1:100),site]
	d_2016 <- allData[["data2016"]]$samplingData[-c(1:200),site]

	if(site!="ER78"){
		d_2015 <- allData[["data2015"]]$samplingData[-c(1:100),site]
		d <- rbind(d_2014,d_2015,d_2016)
	}else{
		d <- rbind(d_2014,d_2016)
	}

	return(d)
}

getSiteDO <- function(allDO, allTemp, site_){

	d_do <- getSiteSeries(allDO, site_)
	d_temp <- getSiteSeries(allTemp, site_)
	
	

	if(site_!="ER78"){
		e <- sapply(allDO, function(x) subset(x$loggerInfo,site == site_)[,"bathymetry"]) + 173
		d_elevation <- data.frame(year = c(2014,2015,2016), elevation = e)
	}else{
		e <- sapply(allDO[c(1,3)], function(x) subset(x$loggerInfo,site == site_)[,"bathymetry"]) + 173
		d_elevation <- data.frame(year = c(2014,2016), elevation = e)
	}


	d <- cbind(d_do, d_temp)

	timeIdx <- index(d)

	newTime <- strftime(timeIdx, format = "%m-%d %H:%M:%S") %>% as.POSIXct(format = "%m-%d %H:%M:%S", tz = "GMT")

	newd <- data.frame(d) %>% rename(DO = d_do, Temp = d_temp) %>%
			mutate(Year = as.factor(lubridate::year(timeIdx)), 
					Time = newTime) %>% 
			merge(d_elevation, by.x = "Year", by.y = "year", all.x = TRUE) %>%
			mutate(DOsat= rMR::DO.saturation(DO.mgl = DO, temp.C = Temp, elevation.m = elevation) * 100)

	pdf(paste0("DO_",site_,".pdf"),width = 6, height = 3.5)
	print(qplot(Time, DO, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("DO (mg/L)"))
	dev.off()

	pdf(paste0("Temp_",site_,".pdf"),width = 6, height = 3.5)
	print(qplot(Time, Temp, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("Temp (C)"))
	dev.off()

	pdf(paste0("DOsat_",site_,".pdf"),width = 6, height = 3.5)
	print(qplot(Time, DOsat, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("DOsat (%)"))
	dev.off()

	return(newd)
}

EPAloggerAnalysis <- function(){

	allDO <- list(
		data2014 = getLakeDO(2014, "B", "hourly") %>% filterEPA(),
		data2015 = getLakeDO(2015, "B", "hourly") %>% filterEPA(),
		data2016 = getLakeDO(2016, "B", "hourly") %>% filterEPA()
	)

	allTemp <- list(
		data2014 = getLakeDO(2014, "B", "hourly","Temp") %>% filterEPA(),
		data2015 = getLakeDO(2015, "B", "hourly","Temp") %>% filterEPA(),
		data2016 = getLakeDO(2016, "B", "hourly","Temp") %>% filterEPA()
	)

	for(site in EPASite){
		newd <- getSiteDO(allDO, allTemp, site)
	}

}

