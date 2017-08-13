for(year in 2014:2015){
  NMF_analysis(year, "hourly", explore = FALSE, new = TRUE, method = "snmf/r") # sparse coefficients
  NMF_analysis(year, "hourly", explore = FALSE, new = TRUE, method = "snmf/l") # sparse basis
}

