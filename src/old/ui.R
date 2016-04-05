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
  				selectInput("var", 
  						label = h5("Variable"), 
        				choices = list("Temperature" = "Temp", "Dissolved Oxygen" = "DO"), 
        				selected = "Temp")
  			),
  			column(
  				3,
  				dateInput("myDate", 
  					label = h5("Date"), 
  					value = "2014-08-01")
  			),
        column(
          3,
          selectInput("mapData", 
              label = h5("Spatial Variable"), 
                choices = list("Bathymetry" = "Bathy", "Logger Data (daily mean)" = "logData"), 
                selected = "Bathy")
        )
  		),

  		fluidRow(
  			column(
  				3,
  				selectInput("dataType", 
  						label = h5("Aggregrate Type"),
       					choices = list("Standard Dev" = "STD", "Average"="AVG","RAW"="Raw"),
       					selected="AVG")
  			),

  			column(
  				3,
  				selectInput("GroupRange", 
  						label = h5("Aggregrate Range"),
       					choices = list("Daily" = "daily", "Hourly"="hourly"),
       					selected="daily")
  			),

  			column(
  				3,
  				sliderInput("myHour", 
  					label = h5("Hour (not implemented yet)"), 
  					min = 0, max = 23,value=12)
  			),
			column(
  				2,
  				selectInput("colors", h5("Color Scheme"),
      			rownames(subset(brewer.pal.info, category %in% c("seq", "div"))))
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
          actionButton("ClearAll","clear all")
  				#actionButton("selectAllLogal","select all")

        ),
  			column(
  				7,
          tabsetPanel(
            
            tabPanel("Time Series",
              dygraphOutput('timeSeriesPlot'),
              checkboxInput("scale", "Scaled?", value = FALSE, width = NULL),
              checkboxInput("twoy", "Double Y?", value = FALSE, width = NULL)
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
                textInput("equation",h5("Variogram Expr"),value = "~longitude+latitude",placeholder="detrend: ~longitude+latitude＋bathymetry"))
          )
  			)
  		)
	)
)