rm(list = ls())
## uncomment the following line to set working directory
setwd("/Users/wenzhaoxu/Developer/Hypoxia_Lake_Erie")
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
resultFolder <- paste0(outputBaseName, "results")
if (!dir.exists(outputBaseName)) {
  dir.create(outputBaseName)
}
if (!dir.exists(resultFolder)) {
  dir.create(resultFolder)
}

randomSeed <- 42 # random seed for reproducibility

main <- function(year = 2014, aggType = "daily") {
  trend <- ~ coords[, "x"] + coords[, "y"] + bathymetry + I(bathymetry^2)

  # create a as.lakeDO object that handles data and interpolation
  erieDO <- getLakeDO(year, "B", aggType) %>% na.omit()

  # create interpolation area with grid size of mapDx * mapDy.
  # The interpolation area is created by CONVEX HULL formed by logger locations
  # grid size: 0.025, 0.025 long,lat, from the config.R
  grid <- createGrid(erieDO$loggerInfo, mapDx, mapDy)

  # using IDW to interpolate. nmax is the number of nearest points to be considered
  print("interpolationg using IDW")
  hypoxia_idw <- predict(
    obj = erieDO,
    grid = grid,
    method = "idw",
    predictType = "extent",
    metaFolder = sprintf("%s%d_%s_idw/", outputBaseName, year, aggType),
    nmax = 5
  )

  saveRDS(hypoxia_idw, sprintf("%s%d_%s_idw/extent.rds", outputBaseName, year, aggType))

  for (r in rList) {
    if ("Reml" %in% interpolation_methods) {
      print(paste("r:", r))
      print("interpolating using Reml")

      # Use basis decomposition, then do kriging interpolation on basis coefficient by MLE estimation
      # finally estimate hypoxia extent based on conditional simulation
      hypoxia_reml <- predict(
        obj = erieDO,
        grid = grid,
        method = "Reml",
        predictType = "extent",
        trend = trend,
        r = r,
        totalSim = 1000,
        nmax = 5,
        metaFolder = sprintf("%s%d_%s_Reml_%d/", outputBaseName, year, aggType, r)
      )
      saveRDS(hypoxia_reml, sprintf("%s%d_%s_Reml_%d/extent.rds", outputBaseName, year, aggType, r))
    }

    if ("Baye" %in% interpolation_methods) {
      print("interpolating using Baye")
      # Use basis decomposition, then do interpolation on basis coefficient by bayesian kriging
      # finally estimate hypoxia extent based on conditional simulation
      hypoxia_baye <- predict(
        obj = erieDO,
        grid = grid,
        method = "Baye",
        predictType = "extent",
        trend = trend,
        r = r,
        totalSim = 1000,
        nmax = 5,
        metaFolder = sprintf("%s%d_%s_Baye_%d/", outputBaseName, year, aggType, r)
      )

      saveRDS(hypoxia_baye, sprintf("%s%d_%s_Baye_%d/extent.rds", outputBaseName, year, aggType, r))
    }
  }
}

# erieDO <- getLakeDO(year, "B", "daily") %>% na.omit()
# erieDO$samplingData

###
# Run interpolations
###
# rList <- c(2,5,10,15) # the number of basis function to try, use 10 for now
rList <- c(10)
interpolation_methods <- c("Baye")

# the intermediate results are saved in the output folder
# the final results are saved in the output folder
main(2021, aggType = "hourly")
main(2022, aggType = "hourly")
main(2023, aggType = "hourly")


#### 
# Start analysis
####

year <- 2023
# target_method = "Reml"
target_method <- "Baye"
target_time_agg <- "daily"

# to plot time series of hypoxia extent
getHypoxiaExtent(year, aggType = target_time_agg, method = target_method, r = 10) # result saved to output/results

# plot the DO interpolations for each time
plot_gif(year, aggType = target_time_agg, method = target_method, r = 10)



# To get detail data
intermediate_result_folder = sprintf("%s%d_%s_%s_%d/", outputBaseName, year, target_time_agg, target_method, 10)
basisModel <- readRDS(sprintf("%s/basisModelRes.rds", intermediate_result_folder))
trendModel <- readRDS(sprintf("%s/trendModel.rds", intermediate_result_folder))
hypoxia <- readRDS(sprintf("%s/extent.rds", intermediate_result_folder))
grid_size <- calculate_average_grid_tile_area(trendModel$grid)

# hypoxia$less0 contains the number of grid tiles with DO < 0.01, convert to km2 by multiplying grid_size
timeSeries <- index(basisModel$residuals$samplingData)
hypoxia$less0 <- hypoxia$less0 * grid_size
hypoxia$less2 <- hypoxia$less2 * grid_size
hypoxia$less4 <- hypoxia$less4 * grid_size

# average hypoxia
hypoxia$average_DO

# plot logger locations
erieDO <- getLakeDO(2022, "B", "daily") %>% na.omit()
loggerInfo <- erieDO$loggerInfo

baseMap <- readRDS("./resources/larger_google_map.rds") + labs(x = "Longitude", y = "Latitude") 
p <- baseMap + geom_point(
  aes(longitude, latitude),
  size = I(2),
  color = "black",
  shape = 21,
  data = loggerInfo
)
print(p)
