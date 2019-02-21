
all.reachable.streams<-function(PCA,mobility){
  
  stream.mobility.lst<-list()
  
  for(n in 1:length(PCA)){             # Particular PCA2
    SegNo<-PCA[n]                                      # n
    
    # Upstream
    upstream.mobility<-Ups_Mobility(SegNo = SegNo,max.mobility = mobility)
    
    # Downstream
    downstream<-alldownstream(hierarchy = hierarchy,catchname = SegNo)
    downstream.mobility.1<-downstream[SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(downstream,SEQ.Clip$SegmentNo)]<=mobility]
    downstream.mobility.1<-downstream.mobility.1[!is.na(downstream.mobility.1)]  # remove NA from the vector
    
    mobility.down<-c()
    for(i in 2:length(downstream.mobility.1)){
      
      other.branch<-setdiff(hierarchy$site[downstream.mobility.1[i]==hierarchy$nextds],downstream.mobility.1)   # begin with 2
      max.mobility.temp<-mobility-(SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(downstream.mobility.1[i],SEQ.Clip$SegmentNo)])
      
      if(length(other.branch)==1){
        mobility.up.down<-Ups_Mobility(SegNo = other.branch,max.mobility = max.mobility.temp)
        mobility.down<-append(mobility.down,mobility.up.down)
      }
      
      if(length(other.branch)>1){
        for(j in 1:length(other.branch)){
          mobility.up.down<-Ups_Mobility(SegNo = other.branch[j],max.mobility = max.mobility.temp)
          mobility.down<-append(mobility.down,mobility.up.down)
        }
      }
    }
    
    stream.mobility<-c(upstream.mobility,mobility.down,downstream.mobility.1)
    stream.mobility.lst[[n]]<-stream.mobility
  }
  return(stream.mobility.lst)
}

