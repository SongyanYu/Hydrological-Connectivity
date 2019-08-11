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
#setwd("D:/My Drive/PhD at GU/Part 4 Hydrologic connectivity")

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
PCA.water.only$SegNo<-as.numeric(PCA.water.only$SegNo)
PCA.water.only$Freq<-rowSums(PCA.water.only[,c(1:107)])

#---
# Combine water-only refuges with SEQ river network
#---
library(maptools)
SEQ.networks<-readShapeLines("Data/Shapfile/SEQ_networks_strahler02.shp")
names(SEQ.networks)
nrow(SEQ.networks)

library(dplyr)
SEQ.networks@data<-left_join(SEQ.networks@data,PCA.water.only,by=c("SegmentNo"="SegNo"))
names(SEQ.networks)
#writeLinesShape(SEQ.networks,fn="Data/Shapfile/Threshold of quant 0.5/PCA_Water only")

#---
# prioritise water-only refuges (top 5%, 15% and 25%)
# intersect with species distribution to evaluate representation
#---

# inundated SegNo
innudt.shp<-readShapeLines("Data/Shapfile/Innudt_SegNo.shp")
inundt.SegNo<-innudt.shp$SegmentNo
delete.seg<-c(859398,859529,856156)
inundt.SegNo<-inundt.SegNo[-match(delete.seg,inundt.SegNo)]

PCA1_SegNo<-setdiff(SEQ.networks$SegmentNo,inundt.SegNo)

numenator<-sum(SEQ.networks$Freq[match(PCA1_SegNo,SEQ.networks$SegmentNo)]>=101)
denominator<-sum(SEQ.networks$Freq[match(PCA1_SegNo,SEQ.networks$SegmentNo)]>0)
numenator/denominator   #0.172, 0.272 and 0.491

# prioritise water-only refuges (top 15%, 25% and 50% of number of years meeting sw and fp criteria)
freq.threshold<-c(106,105,101) # top 15%, 25% and 50%, need to check if the threshold is still true.

SEQ.networks$Freq_15[SEQ.networks$Freq<freq.threshold[1]]<-0
SEQ.networks$Freq_15[SEQ.networks$Freq>=freq.threshold[1]]<-1
SEQ.networks$Freq_25[SEQ.networks$Freq<freq.threshold[2]]<-0
SEQ.networks$Freq_25[SEQ.networks$Freq>=freq.threshold[2]]<-1
SEQ.networks$Freq_50[SEQ.networks$Freq<freq.threshold[3]]<-0
SEQ.networks$Freq_50[SEQ.networks$Freq>=freq.threshold[3]]<-1

names(SEQ.networks)
#writeLinesShape(SEQ.networks,fn="Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
SEQ.networks<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
names(SEQ.networks)

# read in species distribution
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
#names(sdm)
species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,18,21,26)]  # delete 4 non-selecte species

# intersect water-only refuges (top 15%, 25% and 50%) with species distribution 
PCA1_SegNo<-SEQ.networks$SegmentNo[SEQ.networks$Freq_15==1]  # top 15%
PCA1_SegNo<-SEQ.networks$SegmentNo[SEQ.networks$Freq_25==1]  # top 25%
PCA1_SegNo<-SEQ.networks$SegmentNo[SEQ.networks$Freq_50==1]  # top 50%

sp.PCA1<-species.distribution.df[na.omit(match(PCA1_SegNo,species.distribution.df$SegNo)),]

# calcualte sp representation
mean(colSums(sp.PCA1[,c(2:26)])/colSums(species.distribution.df[,c(2:26)]))
# sp.representation: 28%, 41% and 70%.
