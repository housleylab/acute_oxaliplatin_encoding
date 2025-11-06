SplitFilterBinding<-function(x) { 
  ############## vibs   ############## 
  # idx<-c(0, cumsum(diff(x$Trigger_Vib) >0.9))
  # x_split<-split(x, idx)
  # x_split<-x_split[c(2)] 
  # res <- lapply(x_split, function(x){ 
  #   x <- top_n(x, -13000)  
  # })
  # res<-bind_rows(res, .id = "column_label")
  # Vibration_time_corrected<-res
  
  ##############  ramps   ############## 
  idx<-c(0, cumsum(diff(x$Trigger_Ramp) >0.9))
  x_split<-split(x, idx)
  # x_split<-x_split[c(2,3,4,5,  6,7,8,9, 10,11,12,13)]
  x_split<-x_split[c(-1)] ##just remove 1st for consistency

  rampSamples=19000
  res <- lapply(x_split, function(x){
    x <- top_n(x, -rampSamples)
  })
  res<-bind_rows(res, .id = "column_label")
  Ramp_time_corrected<-res[(0:(rampSamples*12)),]

  ##############  triangles   ############## 
  idx<-c(0, cumsum(diff(x$Trigger_Tri) >0.9))
  x_split<-split(x, idx)
  x_split<-x_split[c(2:3)]
  triSamples=53000
  
  res <- lapply(x_split, function(x){
    x <- top_n(x, -triSamples)
  })
  res<-bind_rows(res, .id = "column_label")
  Triangles_time_corrected<-res[(0:(triSamples*2)),]
  
  ############## all bind ############## 
  Train_10_g<-rbind(#Vibration_time_corrected,
    Triangles_time_corrected[((triSamples*0):(triSamples*1)),],
    Ramp_time_corrected[((rampSamples*0):(rampSamples*4)),],
    Ramp_time_corrected[((rampSamples*0):(rampSamples*11)),],
    
    ### test/validation ###
    Triangles_time_corrected[((triSamples*1)+1):(triSamples*2),],
    Ramp_time_corrected[(((rampSamples*11)+1):(rampSamples*12)),]
    
                  #   Triangles_time_corrected[((triSamples*0):(triSamples*1)),],
                  #   Ramp_time_corrected[((rampSamples*4):(rampSamples*8)),], 
                  #   Ramp_time_corrected[((rampSamples*1):(rampSamples*4)),],  
                  #   Ramp_time_corrected[((rampSamples*4):(rampSamples*8)),], 
                  #   Ramp_time_corrected[((rampSamples*8):(rampSamples*11)),], 
                  #   #Triangles_time_corrected[((triSamples*0):(triSamples*1)),],
                  #   
                  #   ### test/validation ###
                  #   #Ramp_time_corrected[(((rampSamples*3)+1):(rampSamples*4)-1),],
                  #   #Ramp_time_corrected[(((rampSamples*7)+1):(rampSamples*8)-1),],
                  #   #Ramp_time_corrected[(((rampSamples*11)+1):(rampSamples*12)),],
                  #   
                  # Triangles_time_corrected[((triSamples*1)+1):(triSamples*2),],
                  # Ramp_time_corrected[(((rampSamples*11)+1):(rampSamples*12)),]
                  # 
                  )
  
  
}
