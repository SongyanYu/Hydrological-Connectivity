#--------------------
# Produce "candidate.df" for "04_PCA3_Mobility.R" to use.
# Objective function = stream segment penalty + feature penalty + position penalty
# Water-only refuges as candidate refuges
# RDI as the cost of stream segments(the lower RDI, the less cost)
# BC as the cost associated with positional importance
# Feature penalty for not achieving conservation target for all species
# Use species mobility to separate priority refuges, maximising dispersal potential after flow resumes.
# Date: 04/07/2019
#--------------------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

#-------------
# water-only refuges
#-------------
library(maptools)
temp<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only.shp")
names(temp)
summary(temp$Freq)
candidate.refuges<-temp$SegmentNo[temp$Freq>=1]  # water-only refuges
candidate.df<-data.frame(SegNo=candidate.refuges)

#------------
# RDI
#------------
SEQ.clip<-readShapePoly("Data/Shapfile/SEQ_Clip.shp")
library(dplyr)
RDI.df<-data.frame(SegNo=SEQ.clip$SegmentNo,RDI=SEQ.clip$RDI)
candidate.df<-left_join(candidate.df,RDI.df,by=c("SegNo"))
sum(is.na(candidate.df$RDI))  # should be "0"

#--------------
# Betweenness centrality (BC)
#--------------
BC.mar<-readRDS("Data/R data/BC.maroochy")
BC.sth<-readRDS("Data/R data/BC.southcoast")
BC.pin<-readRDS("Data/R data/BC.pine")
BC.log<-readRDS("Data/R data/BC.logan")
BC.bne<-readRDS("Data/R data/Betweenness Centrality_BNE")

# To make BC for stream segments from different river networks comparable,
# BC values were normalised within each river network before they were compared.
BC.mar$BC.nor<-(BC.mar$BC-min(BC.mar$BC))/(max(BC.mar$BC)-min(BC.mar$BC))
BC.sth$BC.nor<-(BC.sth$BC-min(BC.sth$BC))/(max(BC.sth$BC)-min(BC.sth$BC))
BC.pin$BC.nor<-(BC.pin$BC-min(BC.pin$BC))/(max(BC.pin$BC)-min(BC.pin$BC))
BC.log$BC.nor<-(BC.log$BC-min(BC.log$BC))/(max(BC.log$BC)-min(BC.log$BC))
BC.bne$BC.nor<-(BC.bne$BC-min(BC.bne$BC))/(max(BC.bne$BC)-min(BC.bne$BC))

BC.SEQ<-rbind(BC.mar,BC.sth,BC.pin,BC.log,BC.bne)

candidate.df<-left_join(candidate.df,BC.SEQ[,-2],by="SegNo")

#----------------
# Plotting BC
#---------------_
names(temp)

BC.SEQ$class[BC.SEQ$BC.nor>=0.71]<-1
BC.SEQ$class[BC.SEQ$BC.nor>=0.38&BC.SEQ$BC.nor<0.71]<-2
BC.SEQ$class[BC.SEQ$BC.nor>=0.12&BC.SEQ$BC.nor<0.38]<-3
BC.SEQ$class[BC.SEQ$BC.nor<0.12]<-4

top5<-quantile(BC.SEQ$BC.nor,probs = 0.95)
BC.SEQ$top5[BC.SEQ$BC.nor>=top5]<-1
BC.SEQ$top5[BC.SEQ$BC.nor<top5]<-0

library(dplyr)
temp@data<-left_join(temp@data,BC.SEQ,by=c("SegmentNo"="SegNo"))
writeLinesShape(temp,fn = "Data/Shapfile/Betweenness centrality/SEQ_BC")





