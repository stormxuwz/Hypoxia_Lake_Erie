library(ggmap)

getGoogleMap <- function(year) {
	# function to fix the map rds saved by previous code
	if(year == 2014) {
		longRange <- c(-82.45828, -80.70047)
		latRange <- c(41.34530, 42.65157)
	} else if (year == 2015) {
		longRange <- c(-82.40251, -80.64470)
		latRange <- c(41.38663, 42.69205)
	} else if (year == 2016) {
		longRange <- c(-82.38872, -80.63091)
		latRange <- c(41.40675, 42.71177)
	}
	
	myMap <- get_googlemap(center = c(lon = mean(lonRange), lat = mean(latRange)), crop=TRUE, maptype = "terrain", scale = 2, zoom = 9)
	ggmap(myMap) %>% saveRDS(sprintf("./resources/erieGoogleMap_%d_new.rds",year))
}