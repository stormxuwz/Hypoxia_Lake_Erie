require(RMySQL)
require(reshape2)
require(zoo)
source("src/helper.R")

retriveLoggerData_DO_Temp <- function(loggerIndex, year, groupRange, dataType, timeRange = NULL) {
  DOData <- retriveLoggerData(loggerIndex, year, "DO", groupRange, dataType, timeRange = timeRange)
  names(DOData) <- paste("DO", names(DOData), sep = "_")

  tempData <- retriveLoggerData(loggerIndex, year, "Temp", groupRange, dataType, timeRange = timeRange)
  names(tempData) <- paste("Temp", names(tempData), sep = "_")

  timeIndex <- index(DOData)
  data <- cbind(data.frame(DOData), data.frame(tempData))
  data <- zoo(data, order.by = timeIndex)
  return(data)
}

retriveLoggerData <- function(loggerIndex, year, var, groupRange, dataType, timeRange = NULL, transform = TRUE) {
  data <- data.frame()
  sqlRes <- getTimeSeriesSQL(loggerIndex, year, var, groupRange, dataType, timeRange = timeRange)
  sql <- sqlRes$sql
  print(sql)
  timeFormat <- sqlRes$timeFormat
  # print(timeFormat)
  if (!transform) {
    data <- sqlQuery(sql)
    data$Time <- strptime(data$Time, format = timeFormat, tz = "GMT")
  } else {
    data <- sqlQuery(sql) %>% dcast(Time ~ logger, value.var = var)
    data <- zoo(subset(data, select = -Time), order.by = strptime(data$Time, format = timeFormat, tz = "GMT"))
  }

  return(data)
}

retriveGeoData <- function(year, position, loggerIndex = NULL) {
  if (is.null(loggerIndex)) {
    sql <- sprintf("select loggerID,longitude,latitude,bathymetry,site from loggerInfo where available=1 and loggerPosition='%s' and year = %s", position, year)
  } else {
    loggerCondition <- sprintf("(%s)", paste(paste("loggerID =", loggerIndex), collapse = " OR "))

    sql <- sprintf("select loggerID,longitude,latitude,bathymetry,site from loggerInfo where available=1 and year = %s and %s", year, loggerCondition)
  }
  # print(sql)
  res <- sqlQuery(sql) %>% lonlat2UTM()
  res$loggerID <- as.character(res$loggerID)
  return(res)
}


retrivePairLogger <- function(loggerIndex, year) {
  condition <- sprintf("(%s)", paste(loggerIndex, collapse = ","))
  sql <- sprintf("select bottom,upper from loggerBottomUpper where (bottom in %s or upper in %s) and year = %s", condition, condition, year)
  return(sqlQuery(sql))
}


retriveSnapShot <- function(variable, dataType, year, day, hour = NULL, loggerIndex = NULL, timezone = "UTC") {
  loggerCond <- "1=1"
  geoData <- retriveGeoData(year, "B", loggerIndex = loggerIndex)

  if (!is.null(loggerIndex)) {
    # no specific logger info, using all logger
    loggerCond <- paste("logger = ", loggerIndex) %>% paste(collapse = " OR ")
  }


  if (is.null(hour)) {
    # perform daily range
    sql <- sprintf("Select date(Time) as Time, %s(%s) as %s, logger from loggerData_%s where (%s) and date(Time) = '%s' Group by date(Time), logger", dataType, variable, dataType, year, loggerCond, day)
    timeFormat <- "%Y-%m-%d"
  } else {
    # first convert time zone
    time <- as.POSIXct(sprintf("%s %d", day, hour), format = "%Y-%m-%d %H", tz = timezone) %>% strftime(format = "%Y-%m-%d %H", tz = "UTC")
    sql <- sprintf("Select DATE_FORMAT(Time,'%%Y-%%m-%%d %%H') as Time, %s(%s) as %s, logger from loggerData_%s where (%s) and DATE_FORMAT(Time,'%%Y-%%m-%%d %%H') = '%s' Group by DATE_FORMAT(Time,'%%Y-%%m-%%d %%H'),logger", dataType, variable, dataType, year, loggerCond, time)
    timeFormat <- "%Y-%m-%d %H"
  }

  data <- sqlQuery(sql)
  if (nrow(data) < 1) {
    return()
  } else {
    data <- merge(data.frame(data), geoData, by.x = "logger", by.y = "loggerID", all.y = FALSE)
    # print(data$Time)
    data$Time <- as.POSIXct(data$Time, format = timeFormat)
    data$id <- 1:nrow(data)
    return(data)
  }
}

sqlQuery <- function(sql) {
  if (nchar(sql) < 1) {
    stop("wrong sql")
  }
  conn <- dbConnect(MySQL(), dbname = dbConfig$dbname, username = dbConfig$username, password = dbConfig$password, host = dbConfig$host, port = 3306)

  result <- dbGetQuery(conn, sql)
  dbDisconnect(conn)
  return(result)
}

getTimeSeriesSQL <- function(loggerIndex, year, var, groupRange, dataType, timeRange = NULL) {
  loggerCondition <- sprintf("(%s)", paste(paste("logger = '", loggerIndex, "'", sep=""), collapse = " OR "))
  if (is.null(timeRange)) {
    sqlCondition <- loggerCondition
  } else {
    sqlCondition <- sprintf("%s AND Time >= '%s' AND Time < '%s' ", loggerCondition, timeRange[1], timeRange[2])
  }

  if (dataType == "STD") {
    if (groupRange == "daily") {
      sql <- sprintf("Select date(Time) as Time, STD(%s) as %s, logger from loggerData_%s where %s Group by date(Time),logger", var, var, year, sqlCondition)
      timeFormat <- "%Y-%m-%d"
    } else if (groupRange == "hourly") {
      sql <- sprintf(
        "Select DATE_FORMAT(Time,'%%Y-%%m-%%d %%H') as Time, STD(%s) as %s, logger from loggerData_%s where %s Group by DATE_FORMAT(Time,'%%Y-%%m-%%d %%H'), logger", 
        var, var, year, sqlCondition)
      timeFormat <- "%Y-%m-%d %H"
    } else {
      stop("Invalid groupRange type")
    }
  } else if (dataType == "AVG") {
    if (groupRange == "daily") {
      sql <- sprintf("Select date(Time) as Time, AVG(%s) as %s, logger from loggerData_%s where %s Group by date(Time), logger", var, var, year, sqlCondition)
      timeFormat <- "%Y-%m-%d"
    } else if (groupRange == "hourly") {
      sql <- sprintf("Select DATE_FORMAT(Time,'%%Y-%%m-%%d %%H') as Time, AVG(%s) as %s, logger from loggerData_%s where %s Group by DATE_FORMAT(Time,'%%Y-%%m-%%d %%H'),logger", var, var, year, sqlCondition)
      timeFormat <- "%Y-%m-%d %H"
    } else {
      stop("Invalid groupRange type")
    }
  } else {
    sql <- sprintf("Select Time, %s, logger from loggerData_%s where %s", var, year, sqlCondition)
    timeFormat <- "%Y-%m-%d %H:%M:%S"
  }

  return(list(sql = sql, timeFormat = timeFormat))
}



# Following is to send data to database
sendtoSQL_loggerInfo <- function(locationInfo, dbTableName) {
  require(lubridate)
  conn <- dbConnect(MySQL(), dbname = dbConfig$dbname, username = dbConfig$username, password = dbConfig$password, host = dbConfig$host, port = 3306)
  dbWriteTable(conn, dbTableName, locationInfo, append = TRUE, row.names = F)
  dbDisconnect(conn)
}


parsingCSV <- function(fileName) {
  fname <- basename(fileName)
  require(lubridate)
  loggerData <- read.csv(fileName, skip = 1)
  timeZ <- names(loggerData)[2]
  n <- nrow(loggerData)
  loggerData <- loggerData[10:(n - 10), c(2, 3, 4)]
  names(loggerData) <- c("Time", "DO", "Temp")
  loggerData <- subset(loggerData, DO > -100)
  if (max(loggerData$Temp, na.rm = T) > 80) {
    loggerData$Temp <- (loggerData$Temp - 32) * 5 / 9
  }

  if (substr(timeZ, 12, 20) == "GMT.04.00") {
    tzStr <- "etc/GMT+4"
  } else if (substr(timeZ, 12, 20) == "GMT.00.00") {
    tzStr <- "GMT"
  } else {
    print(fname)
    stop("Wrong")
  }
  loggerData$Time <- strptime(loggerData$Time, format = "%m/%d/%y %I:%M:%S %p", tz = tzStr)
  loggerData$Time <- with_tz(loggerData$Time, "GMT")
  loggerData$logger <- as.numeric(substr(fname, 1, 8))
  return(loggerData)
}

sendtoSQL_loggerData <- function(folder, skip = 0, dbTableName = "loggerData_2014") {
  conn <- dbConnect(MySQL(), dbname = dbConfig$dbname, username = dbConfig$username, password = dbConfig$password, host = dbConfig$host, port = 3306)
  fileList <- list.files(folder)
  for (fname in fileList) {
    filePath <- file.path(folder, fname)
    loggerData <- parsingCSV(filePath)
    dbWriteTable(conn, dbTableName, loggerData, append = TRUE, row.names = F)
  }
  dbDisconnect(conn)
}




sendtoSQL_waveData <- function(wave) {
  alldf <- data.frame(Time = c(), uc = c(), vc = c(), wvd = c(), logger = c())
  # wave <- readRDS("2DwaveData.rds")

  for (i in 1:length(wave)) {
    data <- wave[[i]]
    # data$Time <- as.character(data$Time)
    data$logger <- as.numeric(unlist(strsplit(names(wave)[i], "_"))[2])
    names(data)[1] <- "Time"

    # data <- melt(data[,-1],id=c("Time","logger"))
    alldf <- rbind(alldf, data)
  }
  alldf$id <- 1:nrow(alldf)
  dbWriteTable(conn, "waveData", alldf, overwrite = TRUE, row.names = F)
}

sendtoSQL_tempData <- function(tempData) {
  alldf <- data.frame()

  for (i in 1:length(tempData)) {
    data <- tempData[[i]]$myTemp
    # data$Time <- as.character(data$Time)
    data$logger <- as.numeric(unlist(strsplit(names(tempData)[i], "_"))[2])
    names(data)[21] <- "Time"

    # data <- melt(data[,-1],id=c("Time","logger"))
    alldf <- rbind(alldf, data)
  }
  alldf$id <- 1:nrow(alldf)
  dbWriteTable(conn, "3DtempData", alldf, overwrite = TRUE, row.names = F)
}
