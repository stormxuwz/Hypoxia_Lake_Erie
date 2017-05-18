# config.R
dbConfig <- list(dbname = "DO", username="root", password="XuWenzhaO", host="epadodb.cdouolvkxmse.us-west-2.rds.amazonaws.com")
# dbConfig <- list(dbname = "DO", username="root", password="XuWenzhaO", host="127.0.0.1")

varUnit <- list(DO="DO(mg/L)",Temp="Temperature(C)")

# trend <- "1st"
# trend <- as.formula("~ coords+bathymetry+bathymetry^2")
#trend <- as.formula("~x+y+bathymetry+I(bathymetry^2)")
#trend <- ~coords+bathymetry
#trend <- ~coords[,"x"]+coords[,"y"]+I(coords[,"x"]^2)+I(coords[,"y"]^2)+I(coords[,"x"]*coords[,"y"])+bathymetry+I(bathymetry^2)

trend <- ~coords[,"x"]+coords[,"y"]+bathymetry+I(bathymetry^2)
#trend <- ~coords[,"y"]+bathymetry
mapDx <- 0.025
mapDy <- 0.025