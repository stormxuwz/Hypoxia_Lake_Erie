rm(list=ls())
require("ggplot2")
require("gridExtra")
source("settings.R")

meta <- readRDS(paste(data_source,"meta.rds",sep=""))
DO_raw_list <- readRDS(paste(data_source,"rawFile.rds",sep=""))
meta=meta[-c(12,18,19,27),]

stationList <- unique(meta$station)

dynamic_time_warping <- function(DOdata,ref_station,query_station,stationName="tmp"){
    require(dtw)
    
    ref_series=DOdata[,ref_station]
    query_series=DOdata[,query_station]
    
    alignment <- dtw(query_series,ref_series,keep=TRUE,step=rabinerJuangStepPattern(6,"c",T))
    
    # alignment<-dtw(query_series,S[1002:5002],keep=TRUE,step=rabinerJuangStepPattern(6,"c"))
    plot(alignment,type="twoway",xlab=query_station,ylab=ref_station,offset=5)
    
    #plot(dtw(N[1000:5000],S[1002:5002],keep=TRUE,step=rabinerJuangStepPattern(6,"c")),type="twoway",offset=-2);
    return(alignment)
}

dtw_test <- dynamic_time_warping(data_sync_DO,"logger_10523437","logger_10523445")


require("combinat")
combination_list=combn(nrow(meta_B),2)

for(i in 1:ncol(combination_list)){
    
    dynamic_time_warping(combination_list[1,i],combination_list[2,i],stationName)
}
