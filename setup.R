install.packages(c("cowplot", "RMySQL", "rgdal", "float", "leaflet", "ggmap", "cowplot", 
									 "geoR", "RMySQL", "rgdal", "float", "leaflet", "ggmap", "fields", 
									 "raster", "dismo", "rgdal", "sp", "float", "doParallel"))

source("src/helper.R")

# create google basemap, requires google API
# resourceFolder <- "./resources"
# createFolder(resourceFolder)
# for(year in c(2014, 2015, 2016)) {
# 	createGogleMapFiles(year, resourceFolder)
# }
