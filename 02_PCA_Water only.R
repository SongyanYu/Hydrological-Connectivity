#-----------------------------
# This script aims to identify prioritised conservation area (PCA) for the water only scenario
# Author: Songyan Yu
# Date Create: 19/12/2018
#-----------------------------

# The variables you need to have before running this script
# 1) daily.sw.1911.2017.df (surface water persistence, derived from running "01_Surface water persistence.R")
# 2) flow.pulse.1911.2017 (annual number of flow pulses, derived from running "01_Flow pulses.R")

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

sw.threshold<-365  # threshold of surface water persistence 

daily.sw.1911.2017.df[daily.sw.1911.2017.df<sw.threshold]<-0
daily.sw.1911.2017.df[daily.sw.1911.2017.df>=sw.threshold]<-1
daily.sw.1911.2017.df$SegNo<-rownames(daily.sw.1911.2017.df)

fp.threshold<-5
colnames(flow.pulse.1911.2017)
flow.pulse.1911.2017[,c(1:107)][flow.pulse.1911.2017[,c(1:107)]<fp.threshold]<-0
flow.pulse.1911.2017[,c(1:107)][flow.pulse.1911.2017[,c(1:107)]>=fp.threshold]<-1
head(flow.pulse.1911.2017)

row.order<-match(daily.sw.1911.2017.df$SegNo,flow.pulse.1911.2017$SegNo)
flow.pulse.1911.2017.sorted<-flow.pulse.1911.2017[row.order,]

PCA.water.only<-flow.pulse.1911.2017.sorted[,c(1:107)]*daily.sw.1911.2017.df[,c(1:107)]
head(PCA.water.only)
PCA.water.only$SegNo<-daily.sw.1911.2017.df$SegNo

# 1. Combine the PCA results with SEQ river network for visulisation
library(maptools)
SEQ.networks<-readShapeLines("Data/Shapfile/SEQ_networks.shp")
plot(SEQ.networks)
names(SEQ.networks)
PCA.water.only$SegNo<-as.numeric(PCA.water.only$SegNo)

library(dplyr)
SEQ.networks@data<-left_join(SEQ.networks@data,PCA.water.only,by=c("SegmentNo"="SegNo"))

#writeLinesShape(SEQ.networks,fn="Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05")

# 2. Spatial pattern of PCA (top 25%)
PCA.freq<-data.frame(SegNo=PCA.water.only$SegNo,Freq=rowSums(PCA.water.only[,c(1:107)]))
summary(PCA.freq)

freq.threshold.perc<-0.75   # top 25%
freq.threshold<-quantile(PCA.freq$Freq,probs = freq.threshold.perc)

PCA.freq$Freq[PCA.freq$Freq<freq.threshold]<-0
PCA.freq$Freq[PCA.freq$Freq>=freq.threshold]<-1

# 2.1 Combine the PCA results with SEQ river network for visulisation
SEQ.networks@data<-left_join(SEQ.networks@data,PCA.freq,by=c("SegmentNo"="SegNo"))
names(SEQ.networks)

#writeLinesShape(SEQ.networks,fn="Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05")

# 3. Temporal dynamic of PCA
total.length<-sum(SEQ.networks$RCHLEN)

prop.PCA<-c()
for(i in 128:234){
  prop.PCA[i-127]<-sum(SEQ.networks$RCHLEN[SEQ.networks@data[,i]==1])/total.length
}

prop.PCA.df<-data.frame(Year=1911:2017,Prop_PCA=prop.PCA)
write.csv(prop.PCA.df,file = "Anal/Proportion PCA 1911_2017.csv")

#which(prop.PCA==max(prop.PCA))+1910

