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
  				selectInput("year", 
  						label = h5("Year"),
       					choices = list("2014" = 2014, "2015"=2015,"2016"=2016),
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
  			),
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

        ),
  			column(
  				7,
          tabsetPanel(
            
            tabPanel("Time Series",
              dygraphOutput('timeSeriesPlot'),
              checkboxInput("withUpperLogger", "With upper logger", value = FALSE, width = NULL),
              checkboxInput("outlier", "Show outlier", value = FALSE, width = NULL)
              ),
            
            tabPanel("Correlation",
              tableOutput('corr')),

            tabPanel("Variogram",
                plotlyOutput("Variogram"),
                textInput("equation",h5("Variogram Expr"),value = "~longitude+latitude",placeholder="detrend: ~longitude+latitude＋bathymetry")),
            
            tabPanel("Daily Hypoxia Extent",
                dygraphOutput('hypoxiaExtentPlot'),
                checkboxInput("showArea", "Using hypoxia area (km^2)", value = TRUE, width = NULL),
               
                downloadButton("downloadHypoxia", label = "Download Hypoxia Area", class = NULL)

              )
            
          )
  			)
  		)
	)
)