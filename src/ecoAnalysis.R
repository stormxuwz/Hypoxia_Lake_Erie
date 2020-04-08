
require(chron)

buoyDataFolder <- "../input/"
EPASite <- c("ER30","ER31","ER32","ER36","ER37","ER38","ER42","ER43","ER73","ER78")

#buoyName <- "45005"
#year <- 2014
#dtype <- "h"

#V1. V2.V3 V4 V5 V6.  V7.   V8.  V9.   V10.  V11 V12.   V13.  V14.  V15. V16.  V17.  V18
#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE

# https://www.ndbc.noaa.gov/measdes.shtml

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


plotBuoyData <- function(year){
	myMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",2014)) + labs(x = "Longitude", y = "Latitude")

	erieDO <- getLakeDO(year, "B", "hourly")
	
	erieEPA <- getLakeDO(2014, "B", "hourly") %>% filterSites(EPASite)
	EPAloggerInfo <- erieEPA$loggerInfo

	nearShoreSite <- list(
		nearShoreSite_2014 = c("ASH003","CLE002","CLE001","GEN001","ESL001","CND001","CND002","CND003"),
		nearShoreSite_2015 = c("CLE003","CLE001","ASH003","PRG001","PST006","WTY004","WTY002","LRG001"),
		nearShoreSite_2016 = c("CBG_83","CBG_43","CBG_55","CBG_94","CBG_11","CBG_52","CBG_73","HER_01","HER_02"))

	
	loggerInfo <- erieDO$loggerInfo

	buoy <- list(pos_45164 = c(41.732, -81.694), 
			pos_45169 = c(41.615, -81.821),
			pos_45005 = c(41.677, -82.398),
			pos_45132 = c(42.460, -81.220),
			pos_45167 = c(42.186, -80.137),
			pos_45176 = c(41.550, -81.765))
	buoyName <- names(buoy)
	dfBuoy <- do.call(rbind.data.frame, buoy)
	names(dfBuoy) <- c("latitude","longitude")
	dfBuoy$site <- names(buoy)
	allB <- sqlQuery( "Select longitude, latitude from loggerInfo")

	lonRange <- range(c(allB$longitude, dfBuoy$longitude))
	latRange <- range(c(allB$latitude, dfBuoy$latitude))
	bbox <- make_bbox(lonRange,latRange,f = 0.2)
	myMap <- get_map(location=bbox, source="google",crop=FALSE) %>% ggmap()
	
	saveRDS(myMap,"erieGoogleMap.rds")


	pdf("EPA_sites.pdf")
	myMap + geom_point(aes(longitude, latitude, color = bathymetry), data = EPAloggerInfo, size = I(4)) + 
	scale_color_gradientn(colours=terrain.colors(10))
		# geom_text(aes(longitude, latitude, label = site), data = EPAloggerInfo)
	dev.off()

	myMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",2014)) + labs(x = "Longitude", y = "Latitude")
	nonEPA <- getLakeDO(2014, "B", "hourly") %>% filterSites(nearShoreSite[[paste0("nearShoreSite_", 2014)]], FALSE)
		myMap + geom_point(aes(longitude, latitude, color = bathymetry), data = nonEPA$loggerInfo, size = I(4)) + 
	scale_color_gradientn(colours=terrain.colors(10), limit = c(-21,-13))
	# geom_point(aes(longitude, latitude), shape = I(4), size = I(4), data = subset(dfBuoy, site %in% c("pos_45164","pos_45167")))

	myMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",2015)) + labs(x = "Longitude", y = "Latitude")
	nonEPA <- getLakeDO(2015, "B", "hourly") %>% filterSites(nearShoreSite[[paste0("nearShoreSite_", 2015)]], FALSE)
		myMap + geom_point(aes(longitude, latitude, color = bathymetry), data = nonEPA$loggerInfo, size = I(4)) + 
	scale_color_gradientn(colours=terrain.colors(10), limit = c(-21,-13))
	# geom_point(aes(longitude, latitude), shape = I(4), size = I(4), data = subset(dfBuoy, site %in% c("pos_45169","pos_45167")))

	myMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",2016)) + labs(x = "Longitude", y = "Latitude")
	nonEPA <- getLakeDO(2016, "B", "hourly") %>% filterSites(nearShoreSite[[paste0("nearShoreSite_", 2016)]], FALSE)
		myMap + geom_point(aes(longitude, latitude, color = bathymetry), data = nonEPA$loggerInfo, size = I(4)) + 
	scale_color_gradientn(colours=terrain.colors(10), limit = c(-21,-13)) 
	# geom_point(aes(longitude, latitude), shape = I(4), size = I(4), data = subset(dfBuoy, site %in% c("pos_45169","pos_45167")))
}



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
		data$WVHT = ifelse(data$WVHT >= 20, NA, data$WVHT)
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


plot_wave_DO <- function(mergedData, variables){
	t <- index(mergedData) %>%  as.POSIXct()
	mergedData <- mergedData[,variables]
	tRange <- na.omit(mergedData) %>% index() %>% range()

	mergedData <- data.frame(mergedData) %>% 
		mutate_all(funs(scale_0_1)) %>%
		mutate(time = t) %>% 
		filter(time > tRange[1] & time < tRange[2]) %>% 
		melt(id.vars = c("time"))

	print(ggplot(mergedData) + geom_line(aes(time, value, color = variable), alpha = 0.85) + 
		theme_bw() + ylab("Scaled value") + xlab("") + theme(legend.position="none"))
}

plot_wave_DO2 <- function(mergedData, variables){
	t <- index(mergedData) %>%  as.POSIXct()
	mergedData <- mergedData[,variables]
	tRange <- na.omit(mergedData) %>% index() %>% range()

	mergedData <- data.frame(mergedData) %>% 
		mutate_all(funs(scale_0_1)) %>%
		mutate(time = t) %>% 
		filter(time > tRange[1] & time < tRange[2])
	names(mergedData) <- c("wave","DO", "time")

	print(ggplot(mergedData) + geom_line(aes(time, DO), color = "cyan4") + geom_bar(aes(time, wave), stat = "identity", alpha = I(0.5)) + 
		theme_bw() + ylab("Scaled value") + xlab("") + theme(legend.position="top") + theme(legend.position="none"))
}


readAllData <- function(year, createNew = TRUE){
	erieDO <- getLakeDO(year, "B", "hourly")
	loggerInfo <- erieDO$loggerInfo

	buoy <- list(pos_45164 = c(41.732, -81.694), 
			pos_45169 = c(41.615, -81.821),
			pos_45005 = c(41.677, -82.398),
			pos_45132 = c(42.460, -81.220),
			pos_45167 = c(42.186, -80.137),
			pos_45176 = c(41.550, -81.765))
	buoyName <- names(buoy)

	dfBuoy <- do.call(rbind.data.frame, buoy)
	names(dfBuoy) <- c("latitude","longitude")
	dfBuoy$buoyName <- buoyName

	if(createNew){
		myMap <- readRDS(sprintf("./resources/erieGoogleMap_%d.rds",year)) + labs(x = "Longitude", y = "Latitude")
		p <- myMap + geom_point(aes(longitude, latitude), data = loggerInfo) + geom_text(aes(longitude, latitude,label = loggerID), data = loggerInfo)
		p <- p + geom_point(aes(longitude, latitude), data = dfBuoy, color = "red") + geom_text(aes(longitude, latitude,label = buoyName), data = dfBuoy)

		pdf(paste0(year, "_plot_buoy_DO.pdf"))
		print(p)
		dev.off()
	# buoyData_cndo1_h <- readBuoyData(year, "cndo1","h")
	# buoyData_faio1_h <- readBuoyData(year, "faio1","h")	# only wind data, temperature data available 
	# buoyData_gelo1_h <- readBuoyData(year, "gelo1","h")	# no data in sampling time in 2015
	# buoyData_cblo1_h <- readBuoyData(year, "cblo1","h") # 2015 data, no southeastShoreWSPD_cblo1_h in sampling time

		buoyData_45005_h <- readBuoyData(year, "45005","h")  # 2014 to 2016
		buoyData_45164_h <- readBuoyData(year, "45164","h")  # 2014 to 2016
		buoyData_45167_h <- readBuoyData(year, "45167","h")  # 2014 to 2016  # eastsouth shore
		buoyData_45132_h <- read45132(year) 

		mergedAll <- merge(erieDO$samplingData, buoyData_45005_h, buoyData_45164_h, buoyData_45167_h, buoyData_45132_h)

		if(year == 2015){
			buoyData_45169_h <- readBuoyData(year, "45169","h") # no 2014 data
			mergedAll <- merge(mergedAll,buoyData_45169_h)	

		}else if(year == 2016){
			buoyData_45169_h <- readBuoyData(year, "45169","h")  # no 2014 data
			buoyData_45176_h <- readBuoyData(year, "45176","h")  # only 2016 exist
			mergedAll <- merge(mergedAll,buoyData_45169_h, buoyData_45176_h)	
		}

		saveRDS(mergedAll, paste0(year,"_wave_DO.rds"))
	}else{
		mergedAll <- readRDS(paste0(year,"_wave_DO.rds"))
	}

	if(year == 2014){
		# SouthWest
		pdf("2014_45164_10384443_SW.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll, c("WVHT_45164_h","10384443"))
		dev.off()
		
		pdf("2014_45164_10384445_SW.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll, c("WVHT_45164_h","10384445"))	
		dev.off()

		pdf("2014_45164_10384437_SW.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll, c("WVHT_45164_h","10384437"))	
		dev.off()


		# SouthEast
		pdf("2014_45167_10384438_SE.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45167_h","10384438"))
		dev.off()

		pdf("2014_45167_10384436_SE.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45167_h","10384436"))
		dev.off()

		#NorthEast

		plot_wave_DO(mergedAll, c("WVHT_c45132","10523446"))
		plot_wave_DO(mergedAll, c("WVHT_c45132","10523443"))
		plot_wave_DO(mergedAll, c("WVHT_c45132","10384439"))
		dev.off()

	}else if(year == 2015){
		
		# plot_wave_DO(mergedAll,c("WVHT_45164_h","10534118")) # no 45164 data
		# plot_wave_DO(mergedAll,c("WVHT_45164_h","10384439"))

		# SouthWest
		pdf("2015_45169_10534118_SW.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45169_h","10534118"))
		dev.off()

		pdf("2015_45169_10384439_SW.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45169_h","10384439"))
		dev.off()
		# mergedAll2 <- data.frame(mergedAll)
		# mergedAll2$time <- index(mergedAll)
		# tRange <- na.omit(mergedAll[,c("WVHT_45169_h","10534118")]) %>% index() %>% range()
		# ggplot(subset(mergedAll2, time >tRange[1] & time < tRange[2])) + 
		# 	geom_bar(aes(time, scale_0_1(WVHT_45169_h)), stat = "identity", color = I("cyan"), alpha = I(0.5)) + 
		# 	geom_line(aes(time, scale_0_1(X10534118)))
		# plot_wave_DO(mergedAll,c("WVHT_45169_h","10384439"))  # not significant

		# West
		pdf("2015_45005_10534123_W.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45005_h","10534123"))
		dev.off()
		
		# SouthEast
		pdf("2015_45167_10461951_SE.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45167_h","10461951"))
		dev.off()

		# NorthEast
		plot_wave_DO(mergedAll,c("WVHT_c45132","10523442"))
		plot_wave_DO(mergedAll,c("WVHT_c45132","10384449"))
		dev.off()
	}else if(year == 2016){
		# SouthWest
		# plot_wave_DO(mergedAll,c("WVHT_45164_h","10534123"))
		pdf("2016_45169_10534123_SW.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45169_h","10534123"))
		dev.off()

		pdf("2016_45169_10384449_SW.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45169_h","10384449"))
		dev.off()


		pdf("2016_45169_10534118_SW.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45169_h","10534118"))
		dev.off()
		
		plot_wave_DO(mergedAll,c("WVHT_45176_h","10534123"))
		
		# SouthEast
		pdf("2016_45167_10534122_SE.pdf",width = 6.25, height = 2)
			plot_wave_DO(mergedAll,c("WVHT_45167_h","10534122"))
		dev.off()
		# NorthEast
		plot_wave_DO(mergedAll,c("WVHT_c45132","10523441"))
		dev.off()
	}
	


	# basis_r <- ifelse(year == 2015, 5, 3)
	# nmfDecomp <- readRDS(sprintf("%s/results/cluster/nmfBasis_%d_%s_%d.rds",outputBaseName, year,"hourly",basis_r))
	# nmfDecomp <- nmfDecomp$basis %>% data.frame()
	# nmfDecomp <- NMF_basis(samplingData,5)$basis %>% data.frame()
	# names(nmfDecomp) <- paste("NMF",1:ncol(nmfDecomp),sep = "_")
	# nmfDecomp <- zoo(nmfDecomp, order.by = t)
	# data <- zoo(data[,c("WDIR","WTMP","ATMP","WVHT","WSPD",southeastShoreWSPD","southwestShoreWSPD")],order.by = time)
	
	# if(year == 2015){
	# 	mergedData <- merge(nmfDecomp, buoyData_45169_h)[,c("NMF_5","WVHT_45169_h")]
	# 	newt <- index(mergedData) %>% as.POSIXct()
	# 	mergedData <- mergedData %>%
	# 		data.frame() %>% 
	# 		rename(Basis_5 = NMF_5) %>% 
	# 		mutate(time = newt) %>% 
	# 		subset(newt >= t[1] & newt <= t[length(t)]) %>% 
	# 		mutate(Basis_5 = scale_0_1(Basis_5),WVHT_45169_h = scale_0_1(WVHT_45169_h)) %>%
	# 		melt(id.vars = c("time"))

	# 	pdf(sprintf("%s/results/%d_nmfBasis_%d_WVHT_%s.pdf",outputBaseName, year, basis_r, "45169h"), width = 6.25, height = 3)
	# 	print(ggplot(mergedData) + geom_line(aes(time, value, color = variable), alpha = 0.85) + 
	# 		theme_bw() + ylab("Scaled value") + xlab("") + theme(legend.position="top"))
	# 	dev.off()
	# 		# geom_line(aes(newt, scale_0_1(WVHT_45169_h)),color = "blue", alpha = 0.5) 


	# 	mergedData <- merge(nmfDecomp, buoyData_45132_h)[,c("NMF_4","WVHT_c45132")]
	# 	newt <- index(mergedData) %>% as.POSIXct()
	# 	mergedData <- mergedData %>% 
	# 		data.frame() %>% 
	# 		rename(Basis_4 = NMF_4) %>% 
	# 		mutate(time = newt) %>% 
	# 		subset(newt >= t[1] & newt <= t[length(t)]) %>% 
	# 		mutate(Basis_4 = scale_0_1(Basis_4),WVHT_c45132 = scale_0_1(WVHT_c45132)) %>%
	# 		melt(id.vars = c("time"))

	# 	pdf(sprintf("%s/results/%d_nmfBasis_%d_WVHT_%s.pdf",outputBaseName, year, basis_r, "c45132"), width = 6.25, height = 3)
	# 	print(ggplot(mergedData) + geom_line(aes(time, value, color = variable), alpha = 0.85) + 
	# 		theme_bw() + ylab("Scaled value") + xlab("") + theme(legend.position="top"))
	# 	dev.off()
	# }else{
	# 	mergedData <- merge(nmfDecomp, buoyData_45164_h)[,c("NMF_2","WVHT_45164_h")]
	# 	newt <- index(mergedData) %>% as.POSIXct()
	# 	mergedData <- mergedData %>%
	# 		data.frame() %>% 
	# 		rename(Basis_2 = NMF_2) %>% 
	# 		mutate(time = newt) %>% 
	# 		subset(newt >= t[1] & newt <= t[length(t)]) %>% 
	# 		mutate(Basis_2 = scale_0_1(Basis_2),WVHT_45164_h = scale_0_1(WVHT_45164_h)) %>%
	# 		melt(id.vars = c("time"))

	# 	pdf(sprintf("%s/results/%d_nmfBasis_%d_WVHT_%s.pdf",outputBaseName, year, basis_r, "45164h"), width = 6.25, height = 3)
	# 	ggplot(mergedData) + geom_line(aes(time, value, color = variable), alpha = 0.85) + 
	# 		theme_bw() + ylab("Scaled value") + xlab("") + theme(legend.position="top")
	# 	dev.off()
	# }



	# fullData <- merge(nmfDecomp, buoyData_45005_h) %>% 
	# 			merge(buoyData_45164_h) %>% merge(buoyData_cndo1_h) %>% merge(buoyData_faio1_h)

	# merged_noApprox <- merge(nmfDecomp,buoyData_45005_h) %>% na.omit()

	# merged <- merge(nmfDecomp,buoyData_45005_h, all = TRUE)  %>% na.approx() %>% na.omit()
	# spectrum(merged$ATMP_45005_h)
	# spectrum(merged$NMF_3)

	# plot(merge(nmfDecomp,buoyData_45132_h)[,c("NMF_4","WVHT_c45132","WVHT2_c45132")] %>% na.omit())

	# plot(merge(nmfDecomp,buoyData_45132_h)[,c("NMF_2","WVHT_45169_h")] %>% na.omit())
	# # plot(fullData[,c("NMF_1","WDIR_45005_h")])

	# ggplot(data = data.frame(fullData))+geom_line(aes(t,scale_0_1(NMF_3))) + geom_line(aes(t, newWDIR),color = "red")
}


# function to analyze the EPA loggers

filterSites <- function(lakeDO, siteList, changeColumnName = TRUE){
	lakeDO$loggerInfo <- subset(lakeDO$loggerInfo, site %in% siteList) %>%
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
	require("cowplot")
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


	p_DO <- qplot(Time, DO, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("DO (mg/L)")
	p_Temp <- qplot(Time, Temp, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("Temperature(C)")

	prow <- plot_grid(
		p_DO+theme(legend.position="none"), 
		p_Temp+theme(legend.position = "none"))
	
	legend_b <- get_legend(p_DO + theme(legend.position="bottom"))
	pdf(paste0("DO_Temp",site_,".pdf"),width = 6, height = 2.5)
	print(plot_grid(legend_b,prow,ncol = 1, rel_heights = c(0.05, 1)))
	dev.off()
	# pdf(paste0("DO_",site_,".pdf"),width = 6, height = 3.5)
	# print(qplot(Time, DO, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("DO (mg/L)"))
	# dev.off()

	# pdf(paste0("Temp_",site_,".pdf"),width = 6, height = 3.5)
	# print(qplot(Time, Temp, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("Temp (C)"))
	# dev.off()

	# pdf(paste0("DOsat_",site_,".pdf"),width = 6, height = 3.5)
	# print(qplot(Time, DOsat, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("DOsat (%)"))
	# dev.off()

	return(newd)
}

EPAloggerAnalysis <- function(){

	allDO <- list(
		data2014 = getLakeDO(2014, "B", "hourly") %>% filterSites(EPASite),
		data2015 = getLakeDO(2015, "B", "hourly") %>% filterSites(EPASite),
		data2016 = getLakeDO(2016, "B", "hourly") %>% filterSites(EPASite)
	)

	allTemp <- list(
		data2014 = getLakeDO(2014, "B", "hourly","Temp") %>% filterEPA(),
		data2015 = getLakeDO(2015, "B", "hourly","Temp") %>% filterEPA(),
		data2016 = getLakeDO(2016, "B", "hourly","Temp") %>% filterEPA()
	)

	for(site in EPASite){
		newd <- getSiteDO(allDO, allTemp, site)
	}

	# for(year in c(2014,2015,2016)){
		tmp <- allDO[[paste0("data",2014)]]$samplingData[-c(1:200),c("ER43","ER32")] %>% na.omit()
		t <- index(tmp) %>% as.POSIXct()
		tmp <- data.frame(tmp) %>% mutate(time = t) %>% melt(id.vars = c("time"))
		p_2014 <- ggplot(data = tmp) + geom_line(aes(time, value,color = variable), alpha = 0.85) + theme_bw() + ylab("DO (mg/L)") + xlab("") +
		 ggtitle(2014) + scale_color_discrete(name = "Station")

		
		tmp <- allDO[[paste0("data",2015)]]$samplingData[-c(1:200),c("ER43","ER32")] %>% na.omit()
		t <- index(tmp) %>% as.POSIXct()
		tmp <- data.frame(tmp) %>% mutate(time = t) %>% melt(id.vars = c("time"))
		p_2015 <- ggplot(data = tmp) + geom_line(aes(time, value,color = variable), alpha = 0.85) + theme_bw() + ylab("DO (mg/L)") + xlab("") +
		 ggtitle(2015) + scale_color_discrete(name = "Station")

			
		tmp <- allDO[[paste0("data",2016)]]$samplingData[-c(1:200),c("ER43","ER32")] %>% na.omit()
		t <- index(tmp) %>% as.POSIXct()
		tmp <- data.frame(tmp) %>% mutate(time = t) %>% melt(id.vars = c("time"))
		p_2016 <- ggplot(data = tmp) + geom_line(aes(time, value,color = variable), alpha = 0.85) + theme_bw() + ylab("DO (mg/L)") + xlab("") +
		 ggtitle(2016) + scale_color_discrete(name = "Station")



		prow <- plot_grid(
		p_2014+theme(legend.position="none"), 
		p_2015+theme(legend.position = "none"),
		p_2016+theme(legend.position = "none"), ncol = 1)

		legend_b <- get_legend(p_2014 + theme(legend.position="bottom"))

		pdf(paste0("ER43_ER32.pdf"), height = 6, width = 5)
		print(plot_grid(prow,legend_b, ncol = 1, rel_heights = c(1,0.05)))
		dev.off()
	# }

	for(year in c(2014,2015,2016)){
		tmp <- allDO[[paste0("data",year)]]$samplingData[-c(1:200),c("ER42","ER30")] %>% na.omit()
		t <- index(tmp) %>% as.POSIXct()
		tmp <- data.frame(tmp) %>% mutate(time = t) %>% melt(id.vars = c("time"))
		p <- ggplot(data = tmp) + geom_line(aes(time, value,color = variable), alpha = 0.85) + theme_bw() + ylab("DO(mg/L)") + xlab("") +
			theme(legend.position="none") + ggtitle(year)

		pdf(paste0(year,"_ER42_ER30.pdf"), height = 2, width = 4)
		print(p)
		dev.off()
	}
	

}




# data2014 = getLakeDO(2014, "B", "hourly")
# data2015 = getLakeDO(2015, "B", "hourly")
# data2016 = getLakeDO(2016, "B", "hourly")

# # SW
# d <- rbind(data2014$samplingData[,"10384437"],data2015$samplingData[,"10534118"],data2016$samplingData[,"10534123"])
# timeIdx <- index(d)
# newTime <- strftime(timeIdx, format = "%m-%d %H:%M:%S") %>% as.POSIXct(format = "%m-%d %H:%M:%S", tz = "GMT")

# newd <- data.frame(d) %>% rename(DO = d) %>%
# 			mutate(Year = as.factor(lubridate::year(timeIdx)), 
# 					Time = newTime) %>% na.omit()
# qplot(Time, DO, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("DO (mg/L)")

# # SE
# d <- rbind(data2014$samplingData[,"10384436"],data2015$samplingData[,"10461951"],data2016$samplingData[,"10534122"])
# timeIdx <- index(d)
# newTime <- strftime(timeIdx, format = "%m-%d %H:%M:%S") %>% as.POSIXct(format = "%m-%d %H:%M:%S", tz = "GMT")

# newd <- data.frame(d) %>% rename(DO = d) %>%
# 			mutate(Year = as.factor(lubridate::year(timeIdx)), 
# 					Time = newTime) %>% na.omit()
# qplot(Time, DO, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("DO (mg/L)")


# # N
# d <- rbind(data2014$samplingData[,"10461951"],data2015$samplingData[,"10534121"],data2016$samplingData[,"10384450"])
# timeIdx <- index(d)
# newTime <- strftime(timeIdx, format = "%m-%d %H:%M:%S") %>% as.POSIXct(format = "%m-%d %H:%M:%S", tz = "GMT")

# newd <- data.frame(d) %>% rename(DO = d) %>%
# 			mutate(Year = as.factor(lubridate::year(timeIdx)), 
# 					Time = newTime) %>% na.omit()
# qplot(Time, DO, color = Year, data = newd,geom = "line",alpha = I(0.85)) + theme_bw() + ylab("DO (mg/L)")



# # data2014 = getLakeDO(2014, "B", "hourly","Temp") 
# data2015 = getLakeDO(2015, "B", "hourly","Temp")
# data2016 = getLakeDO(2016, "B", "hourly","Temp")

