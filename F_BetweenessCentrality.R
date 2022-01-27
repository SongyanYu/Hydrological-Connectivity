# here functions to calculate the metric of betweenness centrality
# one with one core, the other is a paralell calculation version.

betweenness.centrality<-function(SegNo,hierarcy){
  
  pathway.full<-list()
  for(n in 1:3){
    pathway<-lapply(1:length(SegNo),FUN = function(i) pairReaches(hierarchy = hierarchy,site1 = SegNo[n], site2 = SegNo[i]))
    pathway.full<-append(pathway.full,pathway)
    cat(n," out of ", length(SegNo),"\n")
  }
  
  unlist.pathway.full<-unlist(pathway.full)
  bwtweenness.centrality<-sapply(1:length(SegNo),FUN = function(i) sum(unlist.pathway.full==SegNo[i],na.rm=TRUE))
  BC.df<-data.frame(SegNo=SegNo,BC=bwtweenness.centrality)
  return(BC.df)
}


betweenness.centrality.cl<-function(SegNo,hierarcy){
  
  pathway.full<-list()
  
  pathway.full<-parLapply(cl,1:3,fun = function(n){
    pathway<-lapply(1:length(SegNo),FUN = function(i) pairReaches(hierarchy = hierarchy,site1 = SegNo[n], site2 = SegNo[i]))
    pathway.full<-append(pathway.full,pathway)
  })
  
  unlist.pathway.full<-unlist(pathway.full)
  bwtweenness.centrality<-sapply(1:length(SegNo),FUN = function(i) sum(unlist.pathway.full==SegNo[i],na.rm=TRUE))
  BC.df<-data.frame(SegNo=SegNo,BC=bwtweenness.centrality)
  return(BC.df)
}
