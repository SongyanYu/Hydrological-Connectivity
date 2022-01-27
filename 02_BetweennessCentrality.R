#---------------------------
# This script calculate 'betweenness centrality' for aquatic refuge prioritisation.
# author: Songyan Yu
# date updated: 27/01/2022
#---------------------------

# read in shapefiles
SEQ.clip<-maptools::readShapePoly("../../Data/Shapfile/SEQ_Clip.shp")
names(SEQ.clip)
hierarchy<-data.frame(site=SEQ.clip$SegmentNo,nextds=SEQ.clip$DWNID1)

Brisbane.clip<-readShapePoly("../../Data/Shapfile/BrisbaneRiver_Clip.shp")
plot(Brisbane.clip)
Logan.clip<-readShapePoly("../../Data/Shapfile/LoganAlbertRiver_Clip.shp")
Maroochy.clip<-readShapePoly("../../Data/Shapfile/MaroochyRiver_Clip.shp")
Pine.clip<-readShapePoly("../../Data/Shapfile/PineRiver_Clip.shp")
Southcoast.clip<-readShapePoly("../../Data/Shapfile/SouthCoast_Clip.shp")


# calculate betweenness centrality for each individual river network
library(catchstats)
#source("F_catchstatsFunc.R")
source("F_BetweenessCentrality.R")

#SegNo<-Maroochy.clip$SegmentNo
#BC.mar<-betweenness.centrality(SegNo = SegNo,hierarcy = hierarchy)
#saveRDS(BC.mar,file = "BC.maroochy")
BC.mar<-readRDS("../../Data/R_data/BC.maroochy")

#SegNo<-Southcoast.clip$SegmentNo
#BC.sth<-betweenness.centrality(SegNo,hierarchy)
#saveRDS(BC.sth,file="BC.southcoast")
BC.sth<-readRDS("../../Data/R_data/BC.southcoast")

#SegNo<-Pine.clip$SegmentNo
#BC.pin<-betweenness.centrality(SegNo,hierarchy)
#saveRDS(BC.pin,file = "BC.pine")
BC.pin<-readRDS("../../Data/R_data/BC.pine")

#SegNo<-Logan.clip$SegmentNo
#BC.log<-betweenness.centrality(SegNo,hierarchy)
#saveRDS(BC.log,file = "BC.logan")
BC.log<-readRDS("../../Data/R_data/BC.logan")

#SegNo<-Brisbane.clip$SegmentNo
#BC.bne<-betweenness.centrality(SegNo,hierarchy)
BC.bne<-readRDS("../../Data/R_data/Betweenness Centrality_BNE")

BC.SEQ.raw<-rbind(BC.mar,BC.sth,BC.pin,BC.log,BC.bne)
summary(BC.SEQ.raw$BC)

# To make BC for stream segments from different river networks comparable,
# BC values were normalised within each river network before they were compared.
BC.mar$BC.nor<-(BC.mar$BC-min(BC.mar$BC))/(max(BC.mar$BC)-min(BC.mar$BC))
BC.sth$BC.nor<-(BC.sth$BC-min(BC.sth$BC))/(max(BC.sth$BC)-min(BC.sth$BC))
BC.pin$BC.nor<-(BC.pin$BC-min(BC.pin$BC))/(max(BC.pin$BC)-min(BC.pin$BC))
BC.log$BC.nor<-(BC.log$BC-min(BC.log$BC))/(max(BC.log$BC)-min(BC.log$BC))
BC.bne$BC.nor<-(BC.bne$BC-min(BC.bne$BC))/(max(BC.bne$BC)-min(BC.bne$BC))

BC.SEQ<-rbind(BC.mar,BC.sth,BC.pin,BC.log,BC.bne)

# integrate to river networks
SEQ.networks <- maptools::readShapeLines("../../Data/Shapfile/SEQ_networks_strahler02.shp")
names(SEQ.networks)

library(dplyr)
SEQ.networks@data <-
  SEQ.networks@data %>%
  left_join(., BC.SEQ, by = c("SegmentNo" = "SegNo"))

maptools::writeLinesShape(SEQ.networks, fn = "../../Data/Shapfile/BetweennessCentrality/SEQ_BC")


