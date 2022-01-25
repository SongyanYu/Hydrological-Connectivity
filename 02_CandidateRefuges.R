#-----------------------------
# This script identify candidate refuges.
# Author: Songyan Yu
# Date Create: 19/12/2018
# Date update: 25/01/2022
#-----------------------------

setwd("../../")

# Pre-required variables
# 1) daily.sw.1911.2017.df (surface water persistence)
source("Scripts/Hydrological-Connectivity/01_Surface water persistence.R")
# 2) flow.pulse.1911.2017 (annual number of flow pulses)
source("Scripts/Hydrological-Connectivity/01_Flow pulses.R")

#---
# surface water persistence
#---
sw.threshold<-365  # threshold of surface water persistence 
daily.sw.1911.2017.df[daily.sw.1911.2017.df<sw.threshold]<-0
daily.sw.1911.2017.df[daily.sw.1911.2017.df>=sw.threshold]<-1
daily.sw.1911.2017.df$SegNo<-as.numeric(rownames(daily.sw.1911.2017.df))

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

both.selec<-flow.pulse.1911.2017.sorted[,c(1:107)]*daily.sw.1911.2017.df[,c(1:107)]
head(both.selec)
both.selec$SegNo<-as.numeric(daily.sw.1911.2017.df$SegNo)
both.selec$Freq<-rowSums(both.selec[,c(1:107)])

library(maptools)
SEQ.networks<-readShapeLines("Data/Shapfile/SEQ_networks_strahler02.shp")
names(SEQ.networks)
nrow(SEQ.networks)

library(dplyr)
SEQ.networks@data<-left_join(SEQ.networks@data,both.selec,by=c("SegmentNo"="SegNo"))
names(SEQ.networks)

#---
# exclude inundated segments
#---
# inundated SegNo
innudt.shp<-readShapeLines("Data/Shapfile/Innudt_SegNo.shp")
inundt.SegNo<-innudt.shp$SegmentNo
delete.seg<-c(859398,859529,856156)
inundt.SegNo<-inundt.SegNo[-match(delete.seg,inundt.SegNo)]

candidate.seg<-SEQ.networks$SegmentNo[SEQ.networks$Freq>0]
candidate.seg<-setdiff(candidate.seg,inundt.SegNo)

# form a data frame with frequency
candidate.df <- data.frame(SegNo=candidate.seg)
library(dplyr)
candidate.df<-left_join(candidate.df,both.selec[,c(108,109)],by="SegNo")
sum(is.na(candidate.freq$Freq))  # "0"
