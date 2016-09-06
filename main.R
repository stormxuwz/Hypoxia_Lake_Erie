source("./src/database.R")
source("config.R")
source("./src/plot.R")


# analysis <- function(year){
year <- 2014
bottomLogger <- retriveGeoData(year,"B")
upperLogger <- retriveGeoData(year,"B3")

# loggerData_DO_Temp <- retriveLoggerData_DO_Temp(upperLogger$loggerID,year,"hourly","AVG",c("2014-06-15","2014-10-06"))

# plot_DO_temp(loggerData_DO_Temp)

allPlot <- list()
for(logger in upperLogger$loggerID){
    p <- plot_DO_temp(retriveLoggerData_DO_Temp(c(logger),year,"hourly","raw",NULL))
    allPlot[[paste("logger_",logger,sep="")]] <- p
}


for(logger in bottomLogger$loggerID){
    p <- plot_DO_temp(retriveLoggerData_DO_Temp(c(logger),year,"hourly","raw",NULL))
   allPlot[[paste("logger_",logger,sep="")]] <- p
}

saveRDS(allPlot,"../output/allPlot_dygraphs.rds")

# loggerData_DO <- retriveLoggerData(upperLogger$loggerID,year,"DO","hourly","AVG",c("2014-06-15","2014-10-06"))
# plot_value(loggerData,"dygrphs")



# }

