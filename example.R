rm(list= ls())
## uncomment the following line to set working directory
# setwd("/Users/wenzhaoxu/Developer/Hypoxia_Lake_Erie")
source("src/database.R")
source("src/classPrediction.R")
source("src/classDef.R")
source("src/classReConstruct.R")
source("src/classSummary.R")
source("src/helper.R")
source("src/postAnalysis.R")
source("src/basisDecomposition.R")
source("config.R")

# create output folder
outputBaseName <- "./output/"
resultFolder = paste0(outputBaseName, "results")
if (!dir.exists(outputBaseName)) {dir.create(outputBaseName)}
if (!dir.exists(resultFolder)) {dir.create(resultFolder)}

randomSeed <- 42 # random seed for reproducibility

main <- function(year = 2014, aggType = "daily"){
	trend <- ~coords[,"x"]+ coords[,"y"] + bathymetry + I(bathymetry^2)

	# create a as.lakeDO object that handles data and interpolation 
	erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()
	
	# create interpolation area with grid size of mapDx * mapDy. 
    # The interpolation area is created by CONVEX HULL formed by logger locations
    # grid size: 0.025, 0.025 long,lat, from the config.R
	grid <- createGrid(erieDO$loggerInfo, mapDx, mapDy)
	
	# using IDW to interpolate. nmax is the number of nearest points to be considered
	print("interpolationg using IDW")
	hypoxia_idw <- predict(obj = erieDO, 
				grid = grid, 
				method = "idw", 
				predictType = "extent",
				metaFolder = sprintf("%s%d_%s_idw/",outputBaseName, year, aggType), 
				nmax = 5)

	saveRDS(hypoxia_idw, sprintf("%s%d_%s_idw/extent.rds",outputBaseName, year, aggType))

    for(r in rList){
		
        if("Reml" %in% interpolation_methods) {
            print(paste("r:",r))
            print("interpolating using Reml")
            
            # Use basis decomposition, then do kriging interpolation on basis coefficient by MLE estimation
            # finally estimate hypoxia extent based on conditional simulation
            hypoxia_reml <- predict(obj = erieDO,
                        grid = grid, 
                        method = "Reml", 
                        predictType = "extent",
                        trend = trend, 
                        r = r, 
                        totalSim = 1000,
                        nmax = 5,
                        metaFolder = sprintf("%s%d_%s_Reml_%d/",outputBaseName, year, aggType, r))
            saveRDS(hypoxia_reml, sprintf("%s%d_%s_Reml_%d/extent.rds",outputBaseName, year, aggType, r))
        }

        if("Baye" %in% interpolation_methods) {
        
            print("interpolating using Baye")
            # Use basis decomposition, then do interpolation on basis coefficient by bayesian kriging
            # finally estimate hypoxia extent based on conditional simulation
            hypoxia_baye <- predict(obj = erieDO, 
                        grid = grid, 
                        method ="Baye", 
                        predictType = "extent", 
                        trend = trend, 
                        r = r, 
                        totalSim = 1000,
                        nmax = 5,
                        metaFolder = sprintf("%s%d_%s_Baye_%d/",outputBaseName, year, aggType, r))

            saveRDS(hypoxia_baye, sprintf("%s%d_%s_Baye_%d/extent.rds",outputBaseName, year, aggType, r))
        }
	}
}

# set configurations
# rList <- c(2,5,10,15) # the number of basis function to try, use 10 for now
rList <- c(10)
interpolation_methods <- c("Reml", "Baye")


# run for 2014 daily DO data
# the intermediate results are saved in the output folder
# the final results are saved in the output folder
main(2014, aggType="daily") # you can change the year and aggregation type here. possible aggType: "daily", "hourly"

# plot the DO interpolations for each time
plot_gif(2014, aggType="daily", method="Reml", r=10)


# to retrieve hypoxia extent for a specific year and aggregation type
hypoxia <- readRDS(sprintf("%s%d_%s_Baye_%d/extent.rds",outputBaseName, 2014, "daily", 10))
hypoxia

# to plot time series of hypoxia extent
getHypoxiaExtent(2014, aggType="daily", method="Baye", r=10) # result saved to output/results