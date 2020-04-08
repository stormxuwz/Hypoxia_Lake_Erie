library(shiny)
library(dygraphs)
library(zoo)
library(reshape2)
library(leaflet)
library(RColorBrewer)
library(plotly)
library(sp)
library(gstat)
library(raster)
source("./config.R")

source("src/database.R")
source("src/plot.R")
source("src/classPrediction.R")
source("src/classDef.R")
source("src/classReConstruct.R")
source("src/classSummary.R")
source("src/helper.R")
source("src/basisDecomposition.R")
source("src/outlierDetection.R")

# Define global variables for shiny session
emptyData <- zoo(c(rep(NA, 4)),order.by=as.Date(c("2014-1-1","2014-1-2")))
trend <- ~coords[,"x"]+ coords[,"y"] + bathymetry + I(bathymetry^2)
mapDx <- 0.025
mapDy <- 0.025

varUnit <- list(DO="DO(mg/L)",Temp="Temperature(C)")
method <- "basis_baye" # "basis_MLE"

parseDataTypeInput <- function(dataType){
	if(dataType == "Raw"){
		return(list(groupRange = "", dataType= "Raw"))
	}else{
		dataType <- unlist(strsplit(dataType,"_"))
		return(list(groupRange = dataType[1], dataType= dataType[2]))
	}
}

shinyServer(
	function(input,output,session){
	# Plot the map
	output$mymap <- renderLeaflet({
       	leaflet("mymap") %>% clearShapes() %>% addTiles() %>% fitBounds(-82.41, 41.59,-80.75, 42.43)
    })
	
	observe({
    	updateSelectizeInput(session, 'selectedID', choices = unique(geoData()[,"loggerID"]), selected=NULL,server = FALSE)
	})
	
	# Action -- ClearAll
  observe({
    input$ClearAll
    leafletProxy("mymap") %>% clearPopups()
    updateSelectizeInput(session, 'selectedID', choices = unique(geoData()[,"loggerID"]), selected=NULL,server = FALSE)
  })
  
  # Action -- Interpolation
	interpolation <- reactive({
    spdata <- spatialDataAll()
    myGeoData <- geoData()
 
    if(input$var == "All"){
    	# can only interpolate one varialbe at once
    	return(NULL)
    }

    if(input$mapData != "Interpolated"){
    	return(NULL)
    }

	if(is.null(spdata)){
		# do nothing
	}else{
		names(spdata)[3] = "value"
    	myGeoData <- subset(myGeoData,loggerID %in% spdata$logger)
		
		grid <- createGrid(myGeoData,by.x = 0.01, by.y = 0.01)
		convexIndex <- grid$convexIndex

		coordinates(spdata) = ~x + y
		coordinates(grid) = ~x + y

		grid$pred <- NA
		grid$pred[convexIndex == 1] <-  idw(value~1 , spdata, subset(grid, convexIndex == 1), nmax = 5)$var1.pred

		grid <- grid %>% data.frame() %>% select(longitude,latitude,pred) %>% rasterFromXYZ()
		raster::projection(grid)=CRS("+init=epsg:4326")
		return(grid)
		}
  })

	calHypoxiaExtent <- reactive({
  	year <- input$year
  	print("start calculate hypoxia extent")

		loggerInfo <- retriveGeoData(year,"B")

		erieDO <- getLakeDO(year, "B", "daily") %>% na.omit()
		grid <- createGrid(erieDO$loggerInfo, mapDx, mapDy)
		timeIdx <- index(erieDO$samplingData)

		
		
		stopifnot(method %in% c("idw","basis_MLE","basis_baye"))
		print(paste0("using ", method, " to estimate hypoxia extent"))
		if(method == "idw"){
			hypoxiaExtent <- predict(obj = erieDO, 
				grid = grid, 
				method = "idw", 
				predictType = "extent",
				metaFolder = NULL, 
				nmax = 5)
			hypoxiaExtent <- zoo(hypoxiaExtent, order.by = timeIdx)

		}else if(method == "basis_MLE"){
			hypoxiaExtent <- predict(obj = erieDO,
					grid = grid, 
					method = "Reml", 
					predictType = "extent", 
					trend = trend, 
					r = 10, 
					totalSim = 1000,
					nmax = 5,
					metaFolder = NULL)

		}else if(method == "basis_baye"){
			hypoxiaExtent <- predict(obj = erieDO,
					grid = grid, 
					method = "Baye", 
					predictType = "extent", 
					trend = trend, 
					r = 10, 
					totalSim = 1000,
					nmax = 5,
					metaFolder = NULL)
		}

		if(method != "idw"){
			hypoxiaExtent <- do.call(cbind, hypoxiaExtent) %>% zoo(order.by = timeIdx)
		}

		return(list(hypoxiaExtent = hypoxiaExtent, grid = grid, method = method))
  	})

	# plot the hypoxia extent
	output$hypoxiaExtentPlot <- renderDygraph({
		hypoxiaRes  <- calHypoxiaExtent()
		hypoxia <- hypoxiaRes$hypoxiaExtent
		grid <- hypoxiaRes$grid
		method <- hypoxiaRes$method
		totalArea <- attr(grid, "totalArea")
		interpolatedNum <- sum(grid$convexIndex==1)
		
		# show interpolation area on the map
		grid$pred <- NA
		grid$pred[grid$convexIndex==1] <- 1
		pal <- colorNumeric("black", domain = NULL)(grid$pred)
		grid <- grid[,c("longitude","latitude","pred")] %>% rasterFromXYZ()
		raster::projection(grid)=CRS("+init=epsg:4326")
		
		leafletProxy("mymap") %>% 
			clearControls() %>% 
			clearImages() %>%
			addRasterImage(grid, colors = pal, opacity = 0.2)

		# show hypoxia extent
		hypoxia <- hypoxia / interpolatedNum
		label <- "Hypoxia extent ratio (km^2)"

		if(input$showArea){
			hypoxia <- hypoxia * totalArea
			label <- "Hypoxia extent (km^2)"
		}

		if(method == "idw"){
			return(dygraph(hypoxia) %>% dyRangeSelector(retainDateWindow=TRUE) %>% dyAxis("y", label = label))
			
		}else{
			tmpName <- c("lower","median","upper")
			p <- dygraph(hypoxia) %>% 
				dySeries(paste("less0",tmpName,sep = "."), label = "less0") %>% 
				dySeries(paste("less2",tmpName,sep = "."), label = "less2") %>% 
				dySeries(paste("less4",tmpName,sep = "."), label = "less4")

			return(p %>% dyRangeSelector(retainDateWindow=TRUE) %>% dyAxis("y", label = label) )
		}
	})

	visData <- reactive({
	  	year <- isolate(input$year)
	  	tmpDataType <- parseDataTypeInput(input$dataType)
	  	dataType <- tmpDataType$dataType
	  	groupRange <- tmpDataType$groupRange
		
		loggerIndex <- input$selectedID

		if(input$withUpperLogger){
			loggerIndex <- c(loggerIndex,retrivePairLogger(loggerIndex,year)$upper)
		}
		var <- input$var
		
		if(var =="All"){
			retriveLoggerData_DO_Temp(loggerIndex,year,groupRange,dataType,timeRange=NULL)
		}else{
			retriveLoggerData(loggerIndex,year,var,groupRange,dataType,timeRange=NULL)
		}
	})

	# get geoData
	geoData <- reactive({
		year <- input$year
   		retriveGeoData(year,"B")
 	})

	# color scheme
	colorpal <- reactive({
    if(input$mapData == "Bathy"){
      colorNumeric(input$colors, geoData()$bathymetry)
    }
    else if(input$mapData == "logData"){
      colorNumeric(input$colors, spatialDataAll()[,3])
    }
    else{
    	grid <- interpolation()
    	print(grid)
    	if(is.null(grid)){
    		return(colorNumeric(input$colors, spatialDataAll()[,3]))
    	}else{
    		domain <- range(c(range(values(grid),na.rm = TRUE),range(spatialDataAll()[,3],na.rm = TRUE)))
    		colorNumeric(input$colors, domain, na.color = "transparent")
    	}
    }
  })

	# Map ploting functions
	observe({
	   pal <- colorpal()
		
	  # plot data on the map
	  if(input$mapData=="Bathy"){
		   	mygeodata <- geoData()
		    leafletProxy("mymap", data = mygeodata) %>% 
		    	clearImages() %>% 
		    	clearShapes()%>% 
		    	clearControls() %>% 
		    	addCircles(layerId=~loggerID,lng=~longitude,lat=~latitude,radius = 3000, weight = 1, color = "#777777",fillColor = ~pal(bathymetry), fillOpacity = 0.8) %>% 
		    	addLegend(position = "bottomright",pal = pal, values = ~bathymetry, title = "Bathymetry")
	  }else if(input$mapData == "logData"){
	    spdata <- spatialDataAll()
	    if(is.null(spdata)){
	    	print("no spdata")
	      leafletProxy("mymap", data = spdata) %>% clearImages() %>% clearShapes() %>% clearControls()
	    }else{
	      names(spdata)[3]="val"
	      leafletProxy("mymap", data = spdata) %>% clearImages() %>% clearShapes() %>% clearControls() %>% addCircles(layerId=~logger,lng=~longitude,lat=~latitude,radius = 3000, weight = 1, color = "#777777",fillColor = ~pal(val), fillOpacity = 0.8) %>% addLegend(position = "bottomright",pal = pal, values = ~val, title = "avg")
	    }
	  }

	  else{
	  	# do the interpolation
			spdata <- spatialDataAll()
			grid <- interpolation()

			if(is.null(grid)){
	
			}
			else{
				# adding raster
				names(spdata)[3]="val"
				leafletProxy("mymap") %>% clearControls() %>% clearImages() %>% 
				addRasterImage(grid, colors = pal, opacity = 0.8) %>% 
				addLegend(position = "bottomright",pal = pal, values = values(grid), title = "avg")
	
				# adding the logger locations 
				leafletProxy("mymap", data = spdata) %>% clearShapes() %>%
				addCircles(layerId=~logger,lng=~longitude,lat=~latitude,radius = 3000, weight = 1, color = "#777777",fillColor = ~pal(val), fillOpacity = 0.8)
			}
	  }

  	})

	observe({
		click<-input$mymap_shape_click
   	 	if(is.null(click))
          return()

		  leafletProxy("mymap")%>%addPopups(click$lng,click$lat, paste(click$id),options=popupOptions(maxHeight=20,zoomAnimation=FALSE))

     	# use isolate to avoid repeat call to input$selectedID
    	ID <- isolate(input$selectedID)
    	updateSelectizeInput(session, 'selectedID', choices = unique(isolate(geoData())[,"loggerID"]), selected= c(ID,click$id), server = FALSE)
	})

	calOutliers <- reactive({
		data <- visData()
		outliers <- data
		for(i in 1:ncol(data)){
			print(i)
			outliersIndex <- stlOutlierDetection(data[,i],3) == 3
			outliers[!outliersIndex,i] <- NA
		}
		return(outliers)
	})

	output$timeSeriesPlot <- renderDygraph({
		if(is.null(input$selectedID)){
			leafletProxy("mymap")%>%clearPopups()
			return(dygraph(emptyData) %>% dyRangeSelector())
		}
		
		selectedGeoData <- subset(geoData(),loggerID %in% as.numeric(input$selectedID))
		leafletProxy("mymap") %>% 
			clearPopups()%>%
			addPopups(data=selectedGeoData,lng=~longitude,lat=~latitude,paste(selectedGeoData$loggerID,"(",round(selectedGeoData$bathymetry,1),")"),options=popupOptions(maxHeight=20,zoomAnimation=FALSE))
		
		data <- visData()
		if(input$var == "All"){
			plot_DO_temp(data)
		}else{
			outliers = NULL
			if(input$outlier){
				# show outliers
				outliers <- calOutliers()
			}
			plot_value(data,isolate(varUnit[[input$var]]),"dygrphs", outlierSeries = outliers)
		}

	})

	output$corr <- renderTable({
		if(is.null(input$selectedID))
			return(NULL)
		if(input$dataType=="Raw")
			return(NULL)
		return(cor(visData(),use="pairwise.complete.obs"))
	},rownames = TRUE)

	observe({
	  input$year
	  leafletProxy("mymap") %>% clearPopups()
	})
	
	spatialDataAll <- reactive({
	  # get the spatial data at one certain time point with sp locations
	  var <- input$var
	  if(var == "All")
	  	var <- "DO"
	  QueryDay <- input$myDate
	  print(input$dataType)
	  if(input$dataType == "Raw"){
	  	return(NULL) # Raw data can't visualize spatially on the map since the time is different
	  }
	  groupRange <- parseDataTypeInput(input$dataType)$groupRange

	  if(groupRange == "daily"){
			QueryHour <- NULL
		}else{
			QueryHour <- input$myHour
	  }
	  year <- input$year # change to non-isolate
	  
	  data <- retriveSnapShot(var,"AVG",year, QueryDay,QueryHour, NULL,"America/New_York")
	  return(data)
	})
	
	output$downloadData <- downloadHandler(
		filename = function() {
    		paste('data-', isolate(input$year), '.csv', sep='')
  		},
		content = function(con) {
  			
  			if(length(input$selectedID)==0){
  				write.csv(NULL, con)
  			}else{

  				data <- visData() # data in UTC time
				if(input$outlier){
  					outliers <- calOutliers()
  					data <- combindWithOutlier(data,outliers)
  				}
  				write.zoo(data, con,sep = ",")
  			}
  		}
	)

	output$downloadHypoxia <- downloadHandler(
		filename = function() {
    	paste('Hypoxia Area', isolate(input$year), '.csv', sep='')
		},
		content = function(con) {
			data <- calHypoxiaExtent()[[1]]
			print(data)
			data <- data * attr(data,"totalArea")

			write.zoo(data, con,sep = ",")
		}
	)

	output$Variogram <- renderPlotly({
		if(is.null(input$selectedID)){
			return()
		}
		if(length(input$selectedID)<3){
			return()
		}
		spdata <- spatialDataAll() %>% subset(logger %in% input$selectedID)
		
		#  logger,Time,Temp,longitude,latitude,bathymetry,id
		spdata[,3] <- ifelse(spdata[,3]>0.01,spdata[,3],0.01)
		if(is.null(spdata))
		  return()
		names(spdata)[3]="var"
		coordinates(spdata) = ~longitude+latitude
		projection(spdata) = CRS("+init=epsg:4326")
		eq <- paste("var",input$equation)
		v <- data.frame(variogram(as.formula(eq),data=spdata,cloud=T,cutoff=10000))
		
		v$leftLogger <- spdata$logger[v$left]
		v$rightLogger <- spdata$logger[v$right]

		v$leftValue <- spdata$var[v$left]
		v$rightValue <- spdata$var[v$right]
		print(v)
		p <- plot_ly(v, x = ~dist, y=~gamma, mode="markers",hoverinfo = "text",
          text = ~paste(leftLogger,"(",round(leftValue,2),")--",rightLogger,"(",round(rightValue,2),")",sep=""))

	  p
	})
})