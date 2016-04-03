library(RMySQL)
sqlQuery <- function (sql,year) {
  #conn <- dbConnect(MySQL(), dbname = paste("DO",year,sep=""), username="root", password="XuWenzhaO", host="127.0.0.1", port=3306)
  conn <- dbConnect(MySQL(), dbname = "DO", username="root", password="XuWenzhaO", host="do.cm1qoaxjisxm.us-west-2.rds.amazonaws.com", port=3306)
  result <- dbGetQuery(conn,sql)
  # cat(sql)
  dbDisconnect(conn)
  # return the dataframe
  return(result)
}