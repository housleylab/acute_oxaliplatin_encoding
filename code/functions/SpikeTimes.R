SpikeTimes<-function(x) {     
  x%>%filter(X8.nw.1==1)%>%select(Time)
}