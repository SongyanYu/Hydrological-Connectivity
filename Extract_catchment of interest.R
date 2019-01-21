library(catchstats)
library(maptools)

shapfile<-readShapePoly("H:/My Drive/PhD at GU/Part 2 Daily stream flow simulation/Gauges in Tamar/spatial/AFGH_Catch_Tamar_clean.shp")
pour.seg<-c("5799")

extract.ups<-function(shapfile,pour.seg){
  hierarchy_1<-data.frame(site=shapfile$SegmentNo,nextds=shapfile$DWNID1)
  lst_1<-allupstream(hierarchy_1,pour.seg)
  for(i in 1:length(lst_1)){
    cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",pour.seg,".txt"),append = TRUE)
  }
}