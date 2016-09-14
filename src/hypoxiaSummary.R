## script to summarize the hypoxia area

hypoxiaExtent <- function(DOGrid,threshold = 2){
	# DOGrid is a matrix
	# each row is the interpolation of one time
	# column is the different points
	hypoxiaExtent <- c()
	for(i in 1:nrow(DOGrid)){
		hypoxiaExtent <- c(hypoxiaExtent,sum(DOGrid[i,]<threshold,na.rm = TRUE))
	}
	return(hypoxiaExtent)
}

