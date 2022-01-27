setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

library(maptools)
RDI<-readShapeLines("Data/Shapfile/Species distribution model/Ensemble_forecasts.shp")
RDI$RDI

RDI$RDI.class<-1
RDI$RDI.class[RDI$RDI>0.15&RDI$RDI<=0.29]<-2
RDI$RDI.class[RDI$RDI>0.29&RDI$RDI<=0.50]<-3
RDI$RDI.class[RDI$RDI>0.50]<-4

writeLinesShape(RDI,fn="Data/Shapfile/Species distribution model/Ensemble_forecasts.shp")
