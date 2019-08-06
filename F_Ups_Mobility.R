
Ups_Mobility<-function(SegNo,max.mobility){
  upstream<-allupstream(hierarchy = hierarchy,catchname = SegNo)
  upstream.mobility<-upstream[SEQ.Clip$D2OUTLET[match(upstream,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]<=max.mobility]
  
  if(sum(dam.segment %in% upstream.mobility)){
    
    if(sum(dam.segment %in% upstream.mobility)==1){
      block.segment<-dam.segment[dam.segment %in% upstream.mobility]
      upstream.block<-allupstream(hierarchy = hierarchy,catchname = block.segment)
      upstream.mobility<-setdiff(upstream.mobility,upstream.block)
    }
    
    if(sum(dam.segment %in% upstream.mobility)>1){
      block.segments<-dam.segment[dam.segment %in% upstream.mobility]
      upstream.block<-unique(unlist(list_all_upstream(hierarchy = hierarchy,catchnames = block.segments)))
      upstream.mobility<-setdiff(upstream.mobility,upstream.block)
    }
  }
  
  return(upstream.mobility)
}
