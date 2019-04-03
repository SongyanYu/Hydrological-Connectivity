library(catchstats)
library(maptools)

#shapfile<-readShapePoly("H:/My Drive/PhD at GU/Part 2 Daily stream flow simulation/Gauges in Tamar/spatial/AFGH_Catch_Tamar_clean.shp")
shapfile<-SEQ.Clip
pour.seg<-SegNo

extract.ups(shapfile = shapfile,pour.seg =SegNo)

extract.ups<-function(shapfile,pour.seg){
  hierarchy_1<-data.frame(site=shapfile$SegmentNo,nextds=shapfile$DWNID1)
  lst_1<-allupstream(hierarchy_1,pour.seg)
  for(i in 1:length(lst_1)){
    cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",pour.seg,".txt"),append = TRUE)
  }
}

lst_1<-PCA3
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}
