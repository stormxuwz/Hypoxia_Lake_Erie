library(shiny)
library(dygraphs)
library(leaflet)
library(RColorBrewer)
library(plotly)

shinyUI(
  fluidPage(
  		fluidRow(
  			column(
  				3,
  				# h4("This is a test to control time and variable"),
  				selectInput("year", 
  						label = h5("Year"),
       					choices = list("2014" = 2014, "2015"=2015),
       					selected=2014)
  			),
        column(
          3,
          selectInput("dataType", 
              label = h5("Variable Type"),
                choices = list("RAW"="Raw","Daily Average"="daily_AVG", "Daily Standard Dev" = "daily_STD","Hourly Average"="hourly_AVG", "Hourly Standard Dev" = "hourly_STD"),
                selected="Raw")
        ),
  			
  			column(
  				3,
  				dateInput("myDate", 
  					label = h5("Date"), 
  					value = "2014-08-01")
  			),
        # column(
        #   3,
        #   selectInput("GroupRange", 
        #       label = h5("Aggregrate Range"),
        #         choices = list("Daily" = "daily", "Hourly"="hourly"),
        #         selected="daily")
        # )
         column(
          2,
          selectInput("colors", h5("Color Scheme"),
            rownames(subset(brewer.pal.info, category %in% c("seq", "div"))))
        )
        
  		),

  		fluidRow(
  			column(
          3,
          selectInput("var", 
              label = h5("Logger Variable"), 
                choices = list("Temperature" = "Temp", "Dissolved Oxygen" = "DO","Temperature & DO" = "All"), 
                selected = "Temp")
        ),

        column(
          3,
          selectInput("mapData", 
              label = h5("Data on the Map"), 
                choices = list("Bathymetry" = "Bathy", "Logger Variable" = "logData","IDW Interpolation" = "Interpolated"), 
                selected = "Bathy")
        ),
  			

  			column(
  				3,
  				sliderInput("myHour", 
  					label = h5("Hour"), 
  					min = 0, max = 23,value=12)
  			)
			  
  		),

  		hr(),
  		
  		fluidRow(
  			column(
  				5,
  				leafletOutput("mymap"),
          selectizeInput("selectedID", label=h5("Selected Loggers"), 
            choices=NULL, 
            selected = NULL, 
            multiple = TRUE,
            options = NULL),
          actionButton("ClearAll","clear all"),
          downloadButton("downloadData", label = "Download", class = NULL)

  				#actionButton("selectAllLogal","select all")

        ),
  			column(
  				7,
          tabsetPanel(
            
            tabPanel("Time Series",
              dygraphOutput('timeSeriesPlot'),
              # checkboxInput("scale", "Scaled?", value = FALSE, width = NULL),
              # checkboxInput("twoy", "Double Y?", value = FALSE, width = NULL),
              # checkboxInput("withOther", "With Both Variables", value = FALSE, width = NULL),
              checkboxInput("withUpperLogger", "With upper logger", value = FALSE, width = NULL),
              checkboxInput("outlier", "Show outlier", value = FALSE, width = NULL)
              ),
            
            tabPanel("Correlation",
              tableOutput('corr')),

            tabPanel("Variogram",
                plotlyOutput("Variogram"),
                # dateRangeInput('dateRangeVariogram',
                      # label = paste("variogram date range"),
                    # separator = " to ", format = "mm/dd/yy",
                    # startview = 'year', weekstart =0
                # ),
                textInput("equation",h5("Variogram Expr"),value = "~longitude+latitude",placeholder="detrend: ~longitude+latitude＋bathymetry")),
            
            tabPanel("Daily Hypoxia Extent",
                dygraphOutput('hypoxiaExtentPlot'),
                # actionButton("calHypoxiaButtom","Calculate Hypoxia Extent"),
                checkboxInput("showArea", "Using hypoxia area (km^2)", value = TRUE, width = NULL),
                downloadButton("downloadHypoxia", label = "Download Hypoxia Area", class = NULL)

              ),
              tabPanel("Instructions",
               label = paste("variogram date range")

              )
           

            # tabPanel("Settings",
            #   selectInput("interpolationMethod", 
            #       label = h5("Method"),
            #       choices = list("IDW" = "IDW"),  
            #       selected ="IDW"),
            #   selectInput("inteprolationPara1", 
            #       label = h5("Para1"),
            #       choices = list("Para" = "Para1"),  
            #       selected ="Para1"),
            #   actionButton("Interpolation","Interpolate") # first is the action name, second is the UI name

              # )

          )
  			)
  		)
	)
)