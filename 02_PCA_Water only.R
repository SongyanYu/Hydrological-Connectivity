#-----------------------------
# This script aims to 
# 1) identify water-only refuges based on surface water persistence and hydrological connectivity, and
# 2) intersect water-only refuges with species distribution to evaluate representation.
# Author: Songyan Yu
# Date Create: 19/12/2018
# Date update: 08/08/2019
#-----------------------------

# Pre-required variables
# 1) daily.sw.1911.2017.df (surface water persistence, derived from running "01_Surface water persistence.R")
# 2) flow.pulse.1911.2017 (annual number of flow pulses, derived from running "01_Flow pulses.R")

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

#---
# surface water persistence
#---
sw.threshold<-365  # threshold of surface water persistence 
daily.sw.1911.2017.df[daily.sw.1911.2017.df<sw.threshold]<-0
daily.sw.1911.2017.df[daily.sw.1911.2017.df>=sw.threshold]<-1
daily.sw.1911.2017.df$SegNo<-rownames(daily.sw.1911.2017.df)

#---
# hydrological connectivity (# flow pulses)
#---
fp.threshold<-5  # threshold of flow pulses
colnames(flow.pulse.1911.2017)
flow.pulse.1911.2017[,c(1:107)][flow.pulse.1911.2017[,c(1:107)]<fp.threshold]<-0
flow.pulse.1911.2017[,c(1:107)][flow.pulse.1911.2017[,c(1:107)]>=fp.threshold]<-1
head(flow.pulse.1911.2017)

#---
# combine sw and fp
#---
row.order<-match(daily.sw.1911.2017.df$SegNo,flow.pulse.1911.2017$SegNo)
flow.pulse.1911.2017.sorted<-flow.pulse.1911.2017[row.order,]

PCA.water.only<-flow.pulse.1911.2017.sorted[,c(1:107)]*daily.sw.1911.2017.df[,c(1:107)]
head(PCA.water.only)
PCA.water.only$SegNo<-daily.sw.1911.2017.df$SegNo

#---
# Combine water-only refuges with SEQ river network
#---
library(maptools)
SEQ.networks<-readShapeLines("Data/Shapfile/SEQ_networks.shp")
plot(SEQ.networks)
names(SEQ.networks)
PCA.water.only$SegNo<-as.numeric(PCA.water.only$SegNo)

library(dplyr)
SEQ.networks@data<-left_join(SEQ.networks@data,PCA.water.only,by=c("SegmentNo"="SegNo"))
#writeLinesShape(SEQ.networks,fn="Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05")

#---
# prioritise water-only refuges (top 5%, 15% and 25%)
# intersect with species distribution to evaluate representation
#---

# prioritise water-only refuges (top 5%, 15% and 25% of number of years meeting sw and fp criteria)
PCA.freq<-data.frame(SegNo=PCA.water.only$SegNo,Freq=rowSums(PCA.water.only[,c(1:107)]))
summary(PCA.freq)
top<-c(0.05,0.15,0.25)
freq.threshold<-quantile(PCA.freq$Freq,probs = (1-top))
PCA.freq$Freq[PCA.freq$Freq<freq.threshold]<-0
PCA.freq$Freq[PCA.freq$Freq>=freq.threshold]<-1

SEQ.networks@data<-left_join(SEQ.networks@data,PCA.freq,by=c("SegmentNo"="SegNo"))
names(SEQ.networks)
#writeLinesShape(SEQ.networks,fn="Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05")
PCA.water.only<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05")  # top 25%
names(PCA.water.only)
PCA1_SegNo<-PCA.water.only$SegmentNo[PCA.water.only$Freq==1]

# read in species distribution
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
#names(sdm)
species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,18,21,26)]  # delete 4 non-selecte species

# combine water only refuges with species distribution
sp.PCA1<-species.distribution.df[na.omit(match(PCA1_SegNo,species.distribution.df$SegNo)),]

# calcualte sp representation
mean(colSums(sp.PCA1[,c(2:26)])/colSums(species.distribution.df[,c(2:26)]))






