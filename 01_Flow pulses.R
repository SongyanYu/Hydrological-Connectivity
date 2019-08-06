#---------------------------
# This script aims to calculate the annual number of flow pulses for each stream segments across SEQ.
# Author: Songyan Yu
# Crate date: 17/12/2018
#---------------------------

#setwd("H:/My Drive/PhD at GU")
setwd("D:/New folder/Google Drive/PhD at GU")


# Directly read in the annual flow pulse from 1911-2017
#flow.pulse.1911.2017<-readRDS("Part 4 Hydrologic connectivity/Data/R data/Annual flow pulse 1911_2017")
flow.pulse.1911.2017<-readRDS("Part 4 Hydrologic connectivity/Data/R data/Annual flow pulse 1911_2017_quant0.5")

# Or 
# You can calculate it by yourself
#------------------------ Calculation code--------------------
# chopped daily discharge simulations from AWRA-L without river routing
flow.files<-list.files("Part 3 Surface water availability/ENV_Variables/CatRunoff_1911_2017/",
                       pattern = "chopped_discharge_",full.names = T)

library(hydrostats)
library(plyr)
library(reshape)

flow.pulse.list<-list()
for(i in 1:length(flow.files)){
  flow.daily.annual<-read.csv(flow.files[i],row.names = 1)
  flow.daily.annual[flow.daily.annual[]<0]<-0
  
  temp<-data.frame(t(flow.daily.annual))
  temp$Date=gsub("X","",rownames(temp))
  temp.melt<-melt(temp,id="Date")
  
  temp.ts<-ts.format(temp.melt[,c(1,3)],format = "%Y.%m.%d")
  temp.ts$SegNo<-gsub("X","",temp.melt$variable)
  
  #temp.ddply<-ddply(temp.ts,.(SegNo),function(x) high.spells(x,threshold = 0.01,plot = FALSE,ann.stats = FALSE))
  temp.ddply<-ddply(temp.ts,.(SegNo),function(x) high.spells(x,quant = 0.5,plot = FALSE,ann.stats = FALSE))
  
  flow.pulse.list[[i]]<-temp.ddply$n.events
  
  cat(i," out of ",length(flow.files),"\n")
}

flow.pulse.1911.2017<-data.frame(flow.pulse.list)
colnames(flow.pulse.1911.2017)<-c(1911:2017)
flow.pulse.1911.2017$SegNo=temp.ddply$SegNo

saveRDS(flow.pulse.1911.2017,file = "Part 4 Hydrologic connectivity/Data/R data/Annual flow pulse 1911_2017_quant0.5")
#--------------------- Calculation code end -----------------------

# Combine the fp results with SEQ river network for visulisation

# mean annual number of flow pulses
head(flow.pulse.1911.2017)
mean.fp.df<-data.frame(SegNo=flow.pulse.1911.2017$SegNo,mean.fp=rowMeans(flow.pulse.1911.2017[,-108]))

library(maptools)
SEQ.networks<-readShapeLines("Part 4 Hydrologic connectivity/Data/Shapfile/SEQ_networks.shp")
plot(SEQ.networks)

mean.fp.df$SegNo<-as.numeric(flow.pulse.1911.2017$SegNo)
mean.fp.df$fp.class[mean.fp.df$mean.fp>10]<-1
mean.fp.df$fp.class[mean.fp.df$mean.fp<=10&mean.fp.df$mean.fp>9]<-2
mean.fp.df$fp.class[mean.fp.df$mean.fp<=9&mean.fp.df$mean.fp>8]<-3
mean.fp.df$fp.class[mean.fp.df$mean.fp<=8&mean.fp.df$mean.fp>5]<-4
mean.fp.df$fp.class[mean.fp.df$mean.fp<=5]<-5


library(dplyr)
SEQ.networks@data<-left_join(SEQ.networks@data,mean.fp.df,by=c("SegmentNo"="SegNo"))
names(SEQ.networks)

writeLinesShape(SEQ.networks,fn="Part 4 Hydrologic connectivity/Data/Shapfile/Annual mean FP 1911_2017")

# use flow pulse to indicate potential hydro connectivity
fp.threshold<-5

colnames(flow.pulse.1911.2017)
flow.pulse.1911.2017[,c(1:107)][flow.pulse.1911.2017[,c(1:107)]<fp.threshold]<-0
flow.pulse.1911.2017[,c(1:107)][flow.pulse.1911.2017[,c(1:107)]>=fp.threshold]<-1

library(maptools)
SEQ.networks<-readShapeLines("Part 4 Hydrologic connectivity/Data/Shapfile/SEQ_networks.shp")
plot(SEQ.networks)
names(SEQ.networks)

flow.pulse.1911.2017$SegNo<-as.numeric(flow.pulse.1911.2017$SegNo)

library(dplyr)
SEQ.networks@data<-left_join(SEQ.networks@data,flow.pulse.1911.2017,by=c("SegmentNo"="SegNo"))

writeLinesShape(SEQ.networks,fn="Part 4 Hydrologic connectivity/Data/Shapfile/Annual FP 1911_2017")

