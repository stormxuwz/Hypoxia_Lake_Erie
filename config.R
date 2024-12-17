# Define local mySQL database access and other global variables
dbConfig <- list(dbname = "DO", username="root", password="XuWenzhaO", host="127.0.0.1")
varUnit <- list(DO="DO(mg/L)",Temp="Temperature(C)")

trend <- ~coords[,"x"]+coords[,"y"]+bathymetry+I(bathymetry^2)
mapDx <- 0.025
mapDy <- 0.025

erieBathymetryFile <- "/Users/wenzhaoxu/Developer/Hypoxia_Lake_Erie/input/erie_lld/erie_lld.asc"
