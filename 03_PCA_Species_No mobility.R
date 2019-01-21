#----------------------------
# This script aims to identify PCA by taking species distribution with naive mobility.
# Author: Songyan Yu
# Date create: 21/12/2018
#----------------------------

# Required variables
# 1) PCA.water.only (running "02_PCA_Water only.R")

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

# 1.1 Species distribution modelled by Ross et al. 2016
library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/Ensemble_forecasts.shp")

# 1.2 Annual PDA within SEQ
PCA<-readShapeLines(fn="Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05.shp")

sum(sdm$SEGMENTNO %in% PCA$SegmentNo)==nrow(sdm)  # should be "TRUE", meaning all segments in sdm is a subset of PCA.
plot(sdm)
plot(PCA,add=TRUE)

# 2. identify annual PDA for each fish species (total 25)
names(sdm)
names(sdm)[45]  #"AmbAga"
names(sdm)[73]  #"TanTan"

sdm.species<-data.frame(SegNo=sdm$SEGMENTNO,sdm@data[,c(45:73)])
sdm.threshold<-500   # Equivalent to probability threshold of 0.5

PCA.species.lst<-list()
for(i in 2:ncol(sdm.species)){
  habitat.segment<-sdm.species$SegNo[sdm.species[,i]>sdm.threshold]
  
  PCA.species<-c()
  for(j in 1:107){
    sw.segment<-PCA.water.only$SegNo[PCA.water.only[,j]==1]
    PCA.species<-append(PCA.species,intersect(habitat.segment,sw.segment))
  }
  
  PCA.species.lst[[i-1]]<-table(PCA.species)
  cat(i, " out of ",ncol(sdm.species),"\n")
}

library(dplyr)

PCA.species.df<-data.frame(SegNo=sdm$SEGMENTNO)     # how many times a stream segment had been a habitat for a species.
for(i in 1:length(PCA.species.lst)){
  temp<-data.frame(PCA.species.lst[[i]])
  #is.factor(temp$PCA.species)
  temp$PCA.species<-as.numeric(as.character(temp$PCA.species))
  PCA.species.df<-left_join(PCA.species.df,temp,by=c("SegNo"="PCA.species"))
}

species.name<-colnames(sdm.species)[-1]
colnames(PCA.species.df)[-1]<-species.name
PCA.species.df[is.na(PCA.species.df)]<-0

PCA.species.disc<-PCA.species.df
PCA.species.disc[,-1][PCA.species.disc[,-1]<25]<-4
PCA.species.disc[,-1][PCA.species.disc[,-1]>=25&PCA.species.disc[,-1]<50]<-3
PCA.species.disc[,-1][PCA.species.disc[,-1]>=50&PCA.species.disc[,-1]<75]<-2
PCA.species.disc[,-1][PCA.species.disc[,-1]>=75]<-1

# 4. Combine Ensemble models (river networks) and "PCA.speces.df" together for visualisation.
plot(sdm)
names(sdm)
head(sdm@data)

nrow(sdm)==nrow(PCA.species.disc)
sdm@data<-cbind(sdm@data,PCA.species.disc)

#writeLinesShape(sdm,fn="Data/Shapfile/Species distribution model/PCA_Naive_Species")

# 5. Overall spatial distribution of PCA taking species distribution with naive mobility.
PCA.naive<-unique(names(unlist(PCA.species.lst)))

PCA.species.overall<-data.frame(SegNo=sdm$SEGMENTNO)
PCA.species.overall[match(PCA.naive,PCA.species.overall$SegNo),2]<-1
colnames(PCA.species.overall)[2]<-"PCA_Naive"

sum(sdm$SEGMENTNO==PCA.species.overall$SegNo)  # "3538"

names(sdm)
sdm@data<-cbind(sdm@data,PCA.species.overall)
writeLinesShape(sdm,fn="Data/Shapfile/Species distribution model/PCA_Naive_Species")
