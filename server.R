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
source("./src/database.R")
source("./config.R")
source("./src/plot.R")
source("./src/interpolation.R")
source("./src/outlierDetection.R")

emptyData <- zoo(c(rep(NA, 4)),order.by=as.Date(c("2014-1-1","2014-1-2")))

boxcox <- function(x,lambda){return((x^lambda-1)/lambda)}

# ID <- input$selectedID

shinyServer(function(input,output,session)
{

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
    input$Interpolation
    spdata <- spatialDataAll()
  
    myGeoData <- geoData()
 
    if(input$withOther){
    	# can only interpolate one varialbe at once
    	return(NULL)
    }

    if(input$mapData != "Interpolated"){
    	return(NULL)
    }

	if(is.null(spdata)){
		# if no data is returned, plot the original data
		# leafletProxy("mymap", data = spdata) %>% clearShapes()
		# return()
	}else{
		names(spdata)[3] = "value"
		
    	myGeoData <- subset(myGeoData,loggerID %in% spdata$logger)

		# do interpolation
		grid <- createGrid(myGeoData)
		#print(head(grid))
		#print(head(spdata))

		grid$pred <- spatial_interpolation(spdata,grid,method = "IDW")
		grid <- grid[,c("longitude","latitude","pred")] %>% rasterFromXYZ()
		raster::projection(grid)=CRS("+init=epsg:4326")
		return(grid)
	}
  })

  	calHypoxiaExtent <- reactive({
  		# input$calHypoxiaButtom
  		year <- input$year
  		print("start calculate hypoxia extent")
  		# daily <- input$dataType

		loggerInfo <- retriveGeoData(year,"B")
		data <- retriveLoggerData(loggerInfo$loggerID,year,"DO","daily","AVG",transform = TRUE) %>% na.omit()
		
		## hard code for progress bar
		data <- na.omit(data)  # only remain the time where all data are available
		times <- index(data)
		grid <- createGrid(loggerInfo)  # will also return an area
		
		
		interpolationRes <- matrix(0,nrow = nrow(grid),ncol = nrow(data))
		loggerNames  <- as.numeric(colnames(data))
		print("# of interpolation:")
		print(nrow(data))

		withProgress(message = 'Calculate Hypoxia Extent', value = 0, {
		for(i in 1:nrow(data)){
		# for(i in 1:3){
			subData <- data.frame(logger = loggerNames,DO = as.numeric(data[i,]))
			subData <- merge(subData, loggerInfo,by.x = "logger",by.y = "loggerID") %>% rename(value = DO)

			grid$pred <- spatial_interpolation(subData,grid)
			interpolationRes[,i] <- ifelse(grid$pred<0,0,grid$pred)
			incProgress(1/nrow(data), detail = paste("Calculate for time", times[i]))
		}
		})

		interpolationRes <- interpolationRes[grid$convexIndex == 1,] # remove the locations that are NA

		totalPx <- nrow(interpolationRes)

		hypoxia_2 <- colSums(interpolationRes<2)/totalPx
		hypoxia_0 <- colSums(interpolationRes<0.01)/totalPx
		hypoxia_4 <- colSums(interpolationRes<4)/totalPx

		hypoxiaExtent <- zoo(data.frame(below_0.01 = hypoxia_0, below_2 = hypoxia_2, below_4 = hypoxia_4),order.by = times)

		# attr(hypoxiaExtent,"pixSize") <- attr(grid,"pixSize")
		attr(hypoxiaExtent,"totalArea") <- attr(grid,"totalArea")

		# hExtent <- calulateHypoxiaExtent(data,loggerInfo) # calculate hypoxia extent
		


		return(list(hypoxiaExtent = hypoxiaExtent, grid = grid))
  	})

  	output$hypoxiaExtentPlot <- renderDygraph({

		hypoxiaRes  <- calHypoxiaExtent()
		hypoxia <- hypoxiaRes[[1]]
		grid <- hypoxiaRes[[2]]
		
		grid$pred <- NA
		grid$pred[grid$convexIndex==1] <- 1
		
		pal <- colorNumeric("black", domain = NULL)(grid$pred)

		grid <- grid[,c("longitude","latitude","pred")] %>% rasterFromXYZ()
		raster::projection(grid)=CRS("+init=epsg:4326")
		
		leafletProxy("mymap") %>% clearControls() %>% clearImages() %>% clearShapes()
			addRasterImage(grid, colors = pal, opacity = 0.2)

		if(input$showArea){
			hypoxia <- hypoxia*attr(hypoxia,"totalArea")
			label = "Hypoxia extent (km^2)"
		}else{
			label = "Hypoxia extent ratio (km^2)"
		}
		
		return(dygraph(hypoxia) %>% dyRangeSelector(retainDateWindow=TRUE)) %>% dyAxis("y", label = label)

	})



	visData <- reactive({
	  	year <- isolate(input$year)
	  	dataType <- input$dataType
		loggerIndex <- input$selectedID

		if(input$withUpperLogger){
			loggerIndex <- c(loggerIndex,retrivePairLogger(loggerIndex,year)$upper)
		}
		var <- input$var
		groupRange <- input$GroupRange


		if(input$withOther){
			retriveLoggerData_DO_Temp(loggerIndex,year,groupRange,dataType,timeRange=NULL)
		}else{
			retriveLoggerData(loggerIndex,year,var,groupRange,dataType,timeRange=NULL)
		}
	})

	geoData <- reactive({
		year <- input$year
		# print(year)
   		retriveGeoData(year,"B")
 	})

	colorpal <- reactive({
		# to be extracted into plot.R files
      if(input$mapData == "Bathy"){
        colorNumeric(input$colors, geoData()$bathymetry)
      }
	    else if(input$mapData == "logData"){
	      colorNumeric(input$colors, spatialDataAll()[,3])
	    }
	    else{
	    	grid <- interpolation()

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

	  if(input$mapData=="Bathy"){
		   	mygeodata <- geoData()
		    leafletProxy("mymap", data = mygeodata) %>% clearImages() %>% clearShapes()%>% clearControls() %>% addCircles(layerId=~loggerID,lng=~longitude,lat=~latitude,radius = 3000, weight = 1, color = "#777777",fillColor = ~pal(bathymetry), fillOpacity = 0.8) %>% addLegend(position = "bottomright",pal = pal, values = ~bathymetry, title = "Bathymetry")
	   }
	  else if(input$mapData == "logData"){
	    spdata <- spatialDataAll()
	    # names(spdata[,3]) = "avg"
	    if(is.null(spdata)){
	      leafletProxy("mymap", data = spdata) %>% clearImages() %>% clearShapes() %>% clearControls()
	    }else{
	      names(spdata)[3]="val"
	      leafletProxy("mymap", data = spdata) %>% clearImages() %>% clearShapes() %>% clearControls() %>% addCircles(layerId=~logger,lng=~longitude,lat=~latitude,radius = 3000, weight = 1, color = "#777777",fillColor = ~pal(val), fillOpacity = 0.8) %>% addLegend(position = "bottomright",pal = pal, values = ~val, title = "avg")
	    }
	  }

	  else{
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
		if(input$withOther){
			plot_DO_temp(data)
		}else{
			outliers = NULL
			if(input$outlier){
				# show outliers
				#print(head(data))
				outliersIndex <- stlOutlierDetection(data[,1],3) == 3
				#print(outliersIndex)
				outliers <- data[,1]
				outliers[!outliersIndex] <- NA
				# print(outliers)
			}
			plot_value(data,isolate(varUnit[[input$var]]),"dygrphs", outlierSeries = outliers)
		}

	})

	output$corr <- renderTable({
		if(is.null(input$selectedID))
			return(NULL)
		if(input$dataType=="Raw")
			return(NULL)
		# if(input$var==""){}
		return(cor(visData(),use="pairwise.complete.obs"))
	})

	observe({
	  input$year
	  leafletProxy("mymap") %>% clearPopups()
	})
	

	spatialDataAll <- reactive({
		# get the spatial data at one certain time point with sp locations
	  var <- input$var
	  QueryDay <- input$myDate
	  
	  if(input$GroupRange == "daily"){
			QueryHour <- NULL
		}else{
			QueryHour <- input$myHour
	  }
	  year <- isolate(input$year)
	  
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
  				# print(index(data))
  				write.zoo(data, con,sep = ",")
  			}
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
		# spdata[,3] <- boxcox(spdata[,3],0.2)
		#print(head(spdata))
		# print(spdata$Time)
		if(is.null(spdata))
		  return()
		names(spdata)[3]="var"
		coordinates(spdata)= ~longitude+latitude
		projection(spdata)=CRS("+init=epsg:4326")
		#print(spdata)
		eq <- paste("var",input$equation)
		#print(eq)
		v <- data.frame(variogram(as.formula(eq),data=spdata,cloud=T,cutoff=10000))
		# print(v)
		v$leftLogger <- spdata$logger[v$left]
		v$rightLogger <- spdata$logger[v$right]

		v$leftValue <- spdata$var[v$left]
		v$rightValue <- spdata$var[v$right]
		# saveRDS(v,"test.rds")

		p <- plot_ly(v, x = dist, y=gamma, mode="markers",hoverinfo = "text",
          text = paste(v$leftLogger,"(",round(v$leftValue,2),")--",v$rightLogger,"(",round(v$rightValue,2),")",sep=""))
		p

	})

})