scale2Unit <- function(s){
	return((s-min(s,na.rm=T))/(max(s,na.rm=T)-min(s,na.rm=T)))
}     

ft2meter <- function(ft){
	return(ft/3.2808399)
}