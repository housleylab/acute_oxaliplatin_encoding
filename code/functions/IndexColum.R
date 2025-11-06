IndexColum<-function(x) {     ###creates a column for lapply to easily index off of
  cbind(x, Index_times=x$Time)
}