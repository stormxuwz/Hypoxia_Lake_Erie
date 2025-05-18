# this script is used to detect outlier
# Including (1) a simple moving average


combindWithOutlier <- function(data, outliers) {
  # data and outlierSeries are zoo class
  t <- index(data)
  data <- as.data.frame(data)

  outliers <- as.data.frame(outliers)
  nameTmp <- names(outliers)
  names(outliers) <- paste("outliers", nameTmp, sep = "_")

  data <- cbind(data, outliers)
  # data$outlier <- outlierSeries
  data <- zoo(data, order.by = t)
  return(data)
}

tsDecomposition <- function(x, season) {
  # stl functions to do time series decomposition
  res <- ts(x, frequency = season) %>% stl(s.window = 7, s.degree = 0, robust = TRUE) # s.window = 7 follows the suggestions
  # print(plot(res))
  return(res$time.series)
  # time.series has three columns: seasonal, trend and remainder
}


remainderOutlier <- function(x, threshold, window = 3 * 24) {
  ut <- function(y) {
    # m  <- median(y)
    # median(y) + threshold * median(abs(y - m))
    stats <- boxplot.stats(y)$stats
    stats[5]
  }

  lt <- function(y) {
    # m = median(y)
    # median(y) - threshold * median(abs(y - m))
    stats <- boxplot.stats(y)$stats
    stats[1]
  }

  z_ut <- rollapply(zoo(x), window, ut, align = "right")
  z_lt <- rollapply(zoo(x), window, lt, align = "right")

  z_ut <- c(rep(z_ut[1], window - 1), z_ut)
  z_lt <- c(rep(z_lt[1], window - 1), z_lt)

  # print(z_lt)

  outliers <- (x > z_ut) | (x < z_lt)

  return(outliers)
}


stlOutlierDetection <- function(x, threshold, testSeasons = c(8, 17, 24)) {
  remainder <- c()
  outlier <- c()
  # remove the outlier
  naIndex <- is.na(x)
  finalOutlier <- rep(NA, length(x))
  x <- x[!naIndex]
  for (season in testSeasons) {
    decomp <- tsDecomposition(x, season)
    remainder <- c(remainder, decomp[, 3])
    outlier <- c(outlier, remainderOutlier(decomp[, 3], threshold))
  }
  remainder <- matrix(remainder, ncol = 3)
  # print(head(outlier))
  outlier <- matrix(outlier, ncol = 3)
  outlier <- rowSums(outlier)
  # print(outlier)
  finalOutlier[!naIndex] <- outlier
  return(finalOutlier)
}


outlierTest <- function() {
  library(tsoutliers)
  year <- 2014
  loggerInfo <- retriveGeoData(year, "B")
  data <- retriveLoggerData(loggerInfo$loggerID, year, "DO", "hourly", "AVG", transform = TRUE) %>% na.omit()
  series <- data[, 1:2]
  for (i in 1:ncol(series)) {
    series2 <- series[, i]
    ol <- stlOutlierDetection(series[, i], 3)

    series2[ol < 3] <- NA # 3 means all three tests are and
    plot(series[, i])
    points(series2, col = "red")
  }

  return(ol)
}

# ol <- outlierTest()
