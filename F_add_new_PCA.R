
add_new_PCA<-function(SegNo,reachable.stream,mobility){
  SegNo.missed<-SegNo[!(SegNo %in% reachable.stream)]
  reachable.missed<-all_reachable_PCAs(PCA = SegNo.missed,mobility = mobility)
  
  routine.M.lst<-list()
  for(i in 1:99){
    remaining.seg<-SegNo.missed
    routine.M<-c()
    
    while(length(remaining.seg)>0){
      n.lst<-match(remaining.seg[sample(1:length(remaining.seg),1)],SegNo.missed)
      remaining.seg.next<-setdiff(remaining.seg,reachable.missed[[n.lst]])
      
      if(length(remaining.seg)>length(remaining.seg.next)){
        routine.M<-append(routine.M,SegNo.missed[n.lst])
        remaining.seg<-remaining.seg.next
      }
    }
    routine.M.lst[[i]]<-routine.M
  }
  n.min<-match(min(lengths(routine.M.lst)),lengths(routine.M.lst))
  extra.routine<-routine.M.lst[[n.min]]
  
  return(extra.routine)
}
