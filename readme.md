Code structure:
- ui.R, server.R: RShiny application ui and server code
- config.R: define 
- main.R: main analysis code to calculate results
- bayeSensitivityAnalysis.R: analysis code with different bayesian priors
- src: source code. The code is based on R S3 classes

To run:

(1) run the following cmd to install packages

```
install.packages(c("cowplot", "RMySQL", "rgdal", "float", "leaflet", "ggmap", "geoR", "fields", "raster", "dismo", "sp", "doParallel", "dplyr", "reshape2", "zoo" ,"gstat", "ggplot"))
```

(2) Download and create some resources. 

a. Run the following code to get google map on Lake Erie for further plot and save into `./resources` folder (a google API is needed, so the generated files (e.g. `erieGoogleMap_2014_new.rds`) are provided for convenience)

```
source("src/helper.R")
resourceFolder <- "./resources"
createFolder(resourceFolder)
for(year in c(2014, 2015, 2016)) {
	createGogleMapFiles(year, resourceFolder)
}
```

b. Download Lake Erie bathymetry file from https://www.ngdc.noaa.gov/mgg/greatlakes/erie.html (ARC ASCII version), and modify `config.R` to specify the file location


(3) Create a SQL database to store the data

Download the data from https://www.dropbox.com/sh/m1x90felp9bqp18/AAAkzGFITIIHG8Ld4MN7kA5Oa?dl=0 and use `LakeErieDO.sql` to create a mySQL database. The raw sensor data and sensor locations files are also provided. Modify config.R to give access to the database

The database contains 5 tables, they are
loggerBottomUpper: listing locations relationship of bottom and upper loggers.
loggerData_2014 (loggerData_2015 or loggerData_2016): List all logger DO and temperature data in 2014 (2015 or 2016). Columns are `id`, `Time`, `DO`, `Temp` and `logger`, where the timezone is GMT
loggerInfo: List the location of loggers in each year. Columns are `loggerID`, `latitude`, `longitude`, `loggerPosition`, `available`, `bathymetry`, `year`, `site`, `notedDepth`, where bathymetry is extracted from Lake Erie bathymetry file based on latitude and longitude for 2014 and 2015. For 2016, the offshore logger depth is from operator's notes and nearshore loggers are extracted based on latitude and longtitude. 


(4) To run the Rshiny app, in R concole, type:
```
shiny::runApp('src')
```

(5) To run algorithm without Rshiny app, check `main` function in `main.R`. 

