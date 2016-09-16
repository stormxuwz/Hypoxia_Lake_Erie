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
    #print("spData")
    #print(head(spdata))

    myGeoData <- geoData()
    
    #print(head(spdata))
    #print(head(myGeoData))

    if(input$withDO){
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
		 if("DO" %in% names(spdata)){
    		spdata <- rename(spdata,value = DO)
    	}else{
    		spdata <- rename(spdata,value = Temp)
    	}

    	myGeoData <- subset(myGeoData,loggerID %in% spdata$logger)

		# do interpolation
		grid <- createGrid(myGeoData)
		#print(head(grid))
		#print(head(spdata))

		grid$pred <- spatial_interpolation(spdata,grid,method = "IDW")
		grid <- grid[,c("longitude","latitude","pred")] %>% rasterFromXYZ()
		raster::projection(grid)=CRS("+init=epsg:4326")
		return(grid)

		# pal <- colorNumeric(input$colors, values(grid), na.color = "transparent")
		# leafletProxy("mymap") %>% clearControls() %>% clearImages() %>% 
		# addRasterImage(grid, colors = pal, opacity = 0.8) %>% 
		# addLegend(pal = pal, values = values(grid), title = "avg")
	}
  })


  #observe({
    #input$selectAll
    #updateSelectizeInput(session, 'selectedID', choices = unique(geoData()[,"loggerID"]), selected=unique(isolate(geoData())[,"loggerID"]),server = FALSE)
  #})
  
	visData <- reactive({
	  	year <- isolate(input$year)
	  	dataType <- input$dataType
		loggerIndex <- input$selectedID

		if(input$withUpperLogger){
			loggerIndex <- c(loggerIndex,retrivePairLogger(loggerIndex,year)$upper)
		}
		var <- input$var
		groupRange <- input$GroupRange


		if(input$withDO){
			retriveLoggerData_DO_Temp(loggerIndex,year,groupRange,dataType,timeRange=NULL)
		}else{
			retriveLoggerData(loggerIndex,year,var,groupRange,dataType,timeRange=NULL)
		}
		
	})

	geoData <- reactive({
		year <- input$year
   		retriveGeoData(year,"B") %>% lonlat2UTM()
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

	observe({
	   pal <- colorpal()
	   # print(pal)
	   if(input$mapData=="Bathy"){
		   	mygeodata <- geoData()
		    leafletProxy("mymap", data = mygeodata) %>% clearImages() %>% clearShapes() %>% addCircles(layerId=~loggerID,lng=~longitude,lat=~latitude,radius = 3000, weight = 1, color = "#777777",fillColor = ~pal(bathymetry), fillOpacity = 0.8)
	   }
	  else if(input$mapData == "logData"){
	    spdata <- spatialDataAll()
	    
	    if(is.null(spdata)){
	      leafletProxy("mymap", data = spdata) %>% clearImages() %>% clearShapes()
	    }else{
	      names(spdata)[3]="val"
	      leafletProxy("mymap", data = spdata) %>% clearImages() %>% clearShapes() %>% addCircles(layerId=~logger,lng=~longitude,lat=~latitude,radius = 3000, weight = 1, color = "#777777",fillColor = ~pal(val), fillOpacity = 0.8)
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
			addLegend(pal = pal, values = values(grid), title = "avg")

			# change color 
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


	# observe({
	#     if(input$mapData=="Bathy"){
	#       mydata <- geoData()
	#       if(is.null(mydata)){
	#         leafletProxy("mymap",data=mydata) %>% clearControls()
	#       }else{
	#         leafletProxy("mymap",data=mydata) %>% clearControls() %>% addLegend(position = "bottomright", pal = colorpal(), values = ~bathymetry)
	#       }
	#     }
	#     else if(){
	#       spdata <- spatialDataAll()
	#       if(is.null(spdata)){
	#         leafletProxy("mymap",data=spdata) %>% clearControls()
	#       }else{
	#         names(spdata)[3]="avg"
	#         leafletProxy("mymap",data=spdata) %>% clearControls() %>% addLegend(position = "bottomright", pal = colorpal(), values = ~avg)
	#       }
	#     }
 #  	})



	output$timeSeriesPlot <- renderDygraph({
		if(is.null(input$selectedID)){
			leafletProxy("mymap")%>%clearPopups()
			return(dygraph(emptyData) %>% dyRangeSelector())
		}
		selectedGeoData <- subset(geoData(),loggerID %in% as.numeric(input$selectedID))
		leafletProxy("mymap")%>%clearPopups()%>%addPopups(data=selectedGeoData,lng=~longitude,lat=~latitude,paste(selectedGeoData$loggerID,"(",round(selectedGeoData$bathymetry,1),")"),options=popupOptions(maxHeight=20,zoomAnimation=FALSE))
		data <- visData()
		if(input$withDO){
			plot_DO_temp(data)
		}else{
			plot_value(data,isolate(varUnit[[input$var]]),"dygrphs")
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
	#observe({
	#  if(input$year == "2015")
	#    updateDateInput(session,"myDate",value = "2015-08-01")
	#  else
	#    updateDateInput(session,"myDate",value = "2014-08-01")
	#})
	
	
	spatialDataAll <- reactive({

	  	var <- input$var
	  
	  QueryDay <- input$myDate
	  QueryHour <- input$myHour
	  year <- isolate(input$year)
	  
	  # print("hello")
	  if(input$GroupRange=="daily"){
	    sql <- sprintf("Select date(Time) as Time, AVG(%s) as %s, logger from loggerData_%s where date(Time) = '%s' Group by date(Time),logger",var,var,year,QueryDay)
	  }
	  else{
	    sql <- sprintf("Select date(Time) as Time, AVG(%s) as %s, logger from loggerData_%s where date(Time) = '%s' Group by date(Time),logger",var,var,year,QueryDay)
	    
	  }
	  # print(sql)
	  data <- sqlQuery(sql)
	  if(nrow(data)<1){
	    return()
	  }
	  data <- merge(data,geoData(),by.x = "logger", by.y = "loggerID",all.y = FALSE)
	  data$id <- 1:nrow(data)
	  return(data)
	})
	
	
	
	
	spatialData <- reactive({
		# this is for variogram
		# year <- isolate(input$year)
	 #  	dataType <- input$dataType
		# loggerIndex <- input$selectedID
		# var <- input$var
		# groupRange <- input$GroupRange

		# if(input$GroupRange == "daily"){
		# 	QueryTime <- striptime(input$myDate)
		# 	timeRange <- c(QueryTime,QueryDay+3600*24)
		# }else{
		# 	QueryDay <- striptime(paste(input$myDate),input$myHour),format="%y-%m-%d %H")
		# }

		# data <- retriveSnapShot(loggerIndex,year,var,groupRange,dataType,timeRange)
		tmp <- input$selectedID
		tmp <- paste("logger =",tmp)
		tmp <- paste(tmp,collapse=" OR ")
		var <- input$var

		QueryDay <- input$myDate
		QueryHour <- input$myHour

		startDay <- as.character(input$dateRangeVariogram[1])
		endDay <- as.character(input$dateRangeVariogram[2])
		# print(startDay)
		year <- isolate(input$year)
    
		# print("hello")
		if(input$GroupRange=="daily"){
			# sql <- sprintf("Select date(Time) as Time, AVG(%s) as %s, logger from loggerData_%s where (%s) and Time <= '%s' and Time >= '%s' Group by date(Time),logger",var,var,year,tmp,endDay,startDay)
			sql <- sprintf("Select date(Time) as Time, AVG(%s) as %s, logger from loggerData_%s where (%s) and date(Time) = '%s' Group by date(Time),logger",var,var,year,tmp,QueryDay)
		}
		else{
			sql <- sprintf("Select date(Time) as Time, AVG(%s) as %s, logger from loggerData_%s where (%s) and date(Time) = '%s' Group by date(Time),logger",var,var,year,tmp,QueryDay)
		}
		# print(sql)
		data <- sqlQuery(sql)
		if(nrow(data)<1){
		  return()
		}
		data <- merge(data,geoData(),by.x = "logger", by.y = "loggerID",all.y = FALSE)
		data$id <- 1:nrow(data)
		return(data)
	})

	output$Variogram <- renderPlotly({
		if(is.null(input$selectedID)){
			return()
		}
		if(length(input$selectedID)<3){
			return()
		}
		spdata <- spatialData()
		spdata[,3] <- ifelse(spdata[,3]>0.01,spdata[,3],0.01)
		spdata[,3] <- boxcox(spdata[,3],0.2)
		# print(head(spdata))
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
    #print(v)
		p <- plot_ly(v, x = dist, y=gamma, mode="markers",hoverinfo = "text",
          text = paste(v$leftLogger,"(",round(v$leftValue,2),")--",v$rightLogger,"(",round(v$rightValue,2),")",sep=""))
		p

	})

})