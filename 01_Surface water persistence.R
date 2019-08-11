#------------------------
# This script aims to quantify annual surface water persistence
# Author: Songyan Yu
# Creat date: 18/12/2018
#------------------------

setwd("D:/New folder/Google Drive/PhD at GU")
#setwd("D:/My Drive/PhD at GU")

#daily.sw.1911.2017<-readRDS("Part 3 Surface water availability/Data/Model extrapolation/wetCont-days_1911-2017")
daily.sw.1911.2017<-readRDS("Part 3 Surface water availability/Data/Model extrapolation/wet50-days_1911-2017")

daily.sw.1911.2017.df<-data.frame(daily.sw.1911.2017)
colnames(daily.sw.1911.2017.df)<-c(1911:2017)
sum(is.na(daily.sw.1911.2017.df))  # should be "0"

# Combine the sw results with SEQ river network for visulisation
library(maptools)
SEQ.networks<-readShapeLines("Part 4 Hydrologic connectivity/Data/Shapfile/SEQ_networks.shp")
plot(SEQ.networks)
names(SEQ.networks)

sw.threshold<-365  # threshold of surface water persistence 

daily.sw.1911.2017.df[daily.sw.1911.2017.df<sw.threshold]<-0
daily.sw.1911.2017.df[daily.sw.1911.2017.df>=sw.threshold]<-1
daily.sw.1911.2017.df$SegNo<-as.numeric(rownames(daily.sw.1911.2017.df))

library(dplyr)
SEQ.networks@data<-left_join(SEQ.networks@data,daily.sw.1911.2017.df,by=c("SegmentNo"="SegNo"))

writeLinesShape(SEQ.networks,fn="Part 4 Hydrologic connectivity/Data/Shapfile/Daily SW 1911_2017_v2")
