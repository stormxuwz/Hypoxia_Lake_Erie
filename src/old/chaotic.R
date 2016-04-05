findMinimumPeak <- function(series){
    # Automatically find the minimum value of a series
    backward <- diff(series)
    n <- length(backward)
    
    switches <- backward[1:n-1]*backward[2:n]
    minimumPoints <- which(switches<0)
    
    for(i in minimumPoints){
        if(backward[i]<0 && backward[i+1]>0)
        {
            return(i+1)
        }
    }
    
    return(NULL)
}

transformData <- function(DOseries,delay,dims){
    # unfold the embeded space
    n=length(DOseries)
    
    m=n-(dims-1)*delay
    data=matrix(0,m,dims)
    
    for(d in 1:dims){
        # print(d)
        data[,d]=DOseries[seq(from=1+delay*(d-1),to=n,by=1)][1:m]
    }
    
    return(data)
}


chaoticAnalysis <- function(DOseries){
    require(rgl)
    mutualInfo <- mutual(DOseries,lag.max=200)
    minimumIndexMutual <- findMinimumPeak(mutualInfo)
    theilerWindow <- 3*minimumIndexMutual
    fn <- as.numeric(false.nearest(DOseries,m=20,d=minimumIndexMutual,t=theilerWindow))
    embedDim <- findMinimumPeak(fn)
    lp=lyap_k(DOseries,m=embedDim,d=minimumIndexMutual,k=40,t=theilerWindow,s=4000,eps=30,ref=10000)
    
    embededData <- transformData(DOseries, minimumIndexMutual,embedDim)
    plot3d(embededData[,1],embededData[,2],embededData[,3])
}

