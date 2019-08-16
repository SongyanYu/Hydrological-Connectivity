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

both.selec<-flow.pulse.1911.2017.sorted[,c(1:107)]*daily.sw.1911.2017.df[,c(1:107)]
head(both.selec)
both.selec$SegNo<-daily.sw.1911.2017.df$SegNo
both.selec$SegNo<-as.numeric(both.selec$SegNo)
both.selec$Freq<-rowSums(both.selec[,c(1:107)])

library(maptools)
SEQ.networks<-readShapeLines("Data/Shapfile/SEQ_networks_strahler02.shp")
names(SEQ.networks)
nrow(SEQ.networks)

library(dplyr)
SEQ.networks@data<-left_join(SEQ.networks@data,both.selec,by=c("SegmentNo"="SegNo"))
names(SEQ.networks)

#---
# candidate refuges
#---
# inundated SegNo
innudt.shp<-readShapeLines("Data/Shapfile/Innudt_SegNo.shp")
inundt.SegNo<-innudt.shp$SegmentNo
delete.seg<-c(859398,859529,856156)
inundt.SegNo<-inundt.SegNo[-match(delete.seg,inundt.SegNo)]

candidate.seg<-SEQ.networks$SegmentNo[SEQ.networks$Freq>0]
candidate.seg<-setdiff(candidate.seg,inundt.SegNo)
candidate.freq<-data.frame(SegNo=candidate.seg)
library(dplyr)
candidate.freq<-left_join(candidate.freq,both.selec[,c(108,109)],by="SegNo")
sum(is.na(candidate.freq$Freq))  # "0"

#---
# plot freq fdc
#---
library(hydroTSM)
width=17.47
ratio=0.6
ppi=30
png(filename = "Figures/01_Water Only/HighRefugia356/Candidate refuges_FDC.png",width = width*ppi,height = width*ratio*ppi)
fdc(candidate.freq$Freq,xlab = "Exceedance probability",ylab = "# years meeting both selection criteria",
    lQ.thr = NA,hQ.thr = NA,thr.shw=FALSE,
    main=NULL,lwd = 3,pch = NA,cex.lab = 1.5)
dev.off()

#---
# prioritise water-only refuges
#---
top<-c(0.85,0.75,0.65)
freq.15<-quantile(candidate.freq$Freq,probs = top[1])  # actual top 17.2%
freq.25<-quantile(candidate.freq$Freq,probs = top[2])  # actual top 27.2%
freq.35<-quantile(candidate.freq$Freq,probs = top[3])  # actual top 40.4%
# forcibly change freq.35 to be "104".
freq.35<-104  # actual top 33.8%

freq.15.seg<-candidate.freq$SegNo[candidate.freq$Freq>=freq.15]
length(freq.15.seg)/nrow(candidate.freq)
freq.25.seg<-candidate.freq$SegNo[candidate.freq$Freq>=freq.25&candidate.freq$Freq<freq.15]
freq.35.seg<-candidate.freq$SegNo[candidate.freq$Freq>=freq.35&candidate.freq$Freq<freq.25]

SEQ.networks$Freq_class<-4
SEQ.networks$Freq_class[match(freq.15.seg,SEQ.networks$SegmentNo)]<-1
SEQ.networks$Freq_class[match(freq.25.seg,SEQ.networks$SegmentNo)]<-2
SEQ.networks$Freq_class[match(freq.35.seg,SEQ.networks$SegmentNo)]<-3
names(SEQ.networks)
#writeLinesShape(SEQ.networks,fn="Data/Shapfile/Threshold of quant 0.5/PCA_Water only")

#-----
# when re-visit, start from here to calculate species representation.
#----
library(maptools)
SEQ.networks<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
names(SEQ.networks)

#---
# species representation of water-only refuges
#---
# read in species distribution
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
#names(sdm)
species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,18,21,26)]  # delete 4 non-selecte species
n.sp<-colSums(species.distribution.df[,c(2:26)])

# calcualte sp representation: 25.5% (n=476), 38.0% (n=752) and 47.5% (n=935)
freq.15.seg<-SEQ.networks$SegmentNo[SEQ.networks$Freq_class==1]  # top 15%
freq.25.seg<-SEQ.networks$SegmentNo[SEQ.networks$Freq_class==2]  # top 25%
freq.35.seg<-SEQ.networks$SegmentNo[SEQ.networks$Freq_class==3]  # top 35%

sp.15<-species.distribution.df[na.omit(match(freq.15.seg,species.distribution.df$SegNo)),]
rep.15<-colSums(sp.15[,c(2:26)])/n.sp
rep.mean.15<-mean(colSums(sp.15[,c(2:26)])/n.sp)

sp.25<-species.distribution.df[na.omit(match(c(freq.15.seg,freq.25.seg),species.distribution.df$SegNo)),]
rep.25<-colSums(sp.25[,c(2:26)])/n.sp
rep.mean.25<-mean(colSums(sp.25[,c(2:26)])/n.sp)

sp.35<-species.distribution.df[na.omit(match(c(freq.15.seg,freq.25.seg,freq.35.seg),species.distribution.df$SegNo)),]
rep.35<-colSums(sp.35[,c(2:26)])/n.sp
rep.mean.35<-mean(colSums(sp.35[,c(2:26)])/n.sp)

#---
# ggplot sp.rep
#---
bar.data<-data.frame(rest=1-rep.35,top35=rep.35-rep.25,top25=rep.25-rep.15,top15=rep.15)
bar.data$sp<-substr(rownames(bar.data),1,6)

library(reshape)
bar.melt<-melt(bar.data,id="sp")
library(ggplot2)
ggplot()+geom_bar(data=bar.melt,aes(x=sp,y=value,fill=variable),stat = "identity")

