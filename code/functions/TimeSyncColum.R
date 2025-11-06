TimeSyncColum<-function(x) {     ###creates a column for lapply to easily index off of
  cbind(x, Time_sync=seq.int(nrow(x))/10000)
}
