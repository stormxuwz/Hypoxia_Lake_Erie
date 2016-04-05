library(RMySQL)
sqlQuery <- function (sql,year) {
  conn <- dbConnect(MySQL(), dbname = "DO", username="root", password="XuWenzhaO", host="do.cm1qoaxjisxm.us-west-2.rds.amazonaws.com", port=3306)
  result <- dbGetQuery(conn,sql)
  dbDisconnect(conn)
  return(result)
}