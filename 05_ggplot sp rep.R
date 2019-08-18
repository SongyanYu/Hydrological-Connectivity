#-------------
# ggplot species representaiton across all three methods
# water-only refuges, positional refuges and systematic refuges
# Author: Songyan Yu
# Date create: 16/08/2019
#-------------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

#---
# 1. detailed species representaiton
#---
bar.data.water.only<-readRDS("Data/R data/Bar data_water_only")
bar.data.BC<-readRDS("Data/R data/Bar data_BC")
bar.data.PCA3<-readRDS("Data/R data/Bar data_PCA3")

bar.data.water.only$method<-"Water-only refuges"
bar.data.BC$method<-"Positional refuges"
bar.data.PCA3$method<-"Systematic refuges"

bar.all<-rbind(bar.data.water.only,bar.data.BC,bar.data.PCA3)

library(reshape)
bar.melt<-melt(bar.all,id=c("sp","method"))
is.factor(bar.melt$method)
bar.melt$method<-factor(bar.melt$method,levels=c("Water-only refuges","Positional refuges","Systematic refuges"))
bar.melt$value<-bar.melt$value*100

library(ggplot2)
ggplot()+geom_bar(data=bar.melt,aes(x=sp,y=value,fill=variable),stat="identity")+
  facet_grid(method~ .)+theme_classic()+
  scale_fill_grey(start = 1,end = 0.3,labels=c("","Top35% / 35%","Top25% / 25%","Top15% / 15%"))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust=1))+
  xlab("Fish species")+ylab("% of species total distribution")+labs(fill=c("Threshold / Target"))+
  ggsave(filename = "Figures/03_Objective function/Cumulative representation of each method.png")

#---
# 2. overall species representation (mean)
#---
# read in species distribution
library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
#names(sdm)
species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,18,21,26)]  # delete 4 non-selecte species
n.sp<-colSums(species.distribution.df[,c(2:26)])

#---
# 2.1 sp rep of water-only refuges: 25.5% (n=476), 38.0% (n=752) and 47.5% (n=935)
#---
SEQ.networks<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
names(SEQ.networks)

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

refuge.size.water.only<-data.frame(top15=nrow(sp.15),top25=nrow(sp.25),top35=nrow(sp.35),method="Water-only refuges")
sp.rep.water.only<-data.frame(top15=rep.mean.15,sd15=sd(rep.15),top25=rep.mean.25,sd25=sd(rep.25),top35=rep.mean.35,sd35=sd(rep.35),method="Water-only refuges")

#---
# 2.2 sp rep of positional refuges: 27.6% (n=415), 39.4% (n=693) and 48.9% (n=969)
#---
# read in positional refuges
SEQ.bc<-readShapeLines("Data/Shapfile/Betweenness centrality/SEQ_BC")
names(SEQ.bc)

# the number of positional refuges should be the same as that of water-only refuges 
BC.15.seg<-SEQ.bc$SegmentNo[SEQ.bc$BC_class==1]
BC.25.seg<-SEQ.bc$SegmentNo[SEQ.bc$BC_class==2]
BC.35.seg<-SEQ.bc$SegmentNo[SEQ.bc$BC_class==3]

sp.15<-species.distribution.df[na.omit(match(BC.15.seg,species.distribution.df$SegNo)),]
rep.15<-colSums(sp.15[,c(2:26)])/n.sp
rep.mean.15<-mean(colSums(sp.15[,c(2:26)])/n.sp)

sp.25<-species.distribution.df[na.omit(match(c(BC.25.seg,BC.15.seg),species.distribution.df$SegNo)),]
rep.25<-colSums(sp.25[,c(2:26)])/n.sp
rep.mean.25<-mean(colSums(sp.25[,c(2:26)])/n.sp)

sp.35<-species.distribution.df[na.omit(match(c(BC.35.seg,BC.25.seg,BC.15.seg),species.distribution.df$SegNo)),]
rep.35<-colSums(sp.35[,c(2:26)])/n.sp
rep.mean.35<-mean(colSums(sp.35[,c(2:26)])/n.sp)

refuge.size.positional<-data.frame(top15=nrow(sp.15),top25=nrow(sp.25),top35=nrow(sp.35),method="Positional refuges")
sp.rep.positional<-data.frame(top15=rep.mean.15,sd15=sd(rep.15),top25=rep.mean.25,sd25=sd(rep.25),top35=rep.mean.35,sd35=sd(rep.35),method="Positional refuges")

#---
# 2.3 sp rep of systematic refuges: 15.4% (n=321),26.2% (n=560) and 35.9% (n=801).
#---
SEQ.best<-readShapeLines("Data/Shapfile/PCA3 non Mob/Best solution")
names(SEQ.best)

best.solution.top15<-SEQ.best$SegmentNo[SEQ.best$best_top15==1]
best.solution.top25<-SEQ.best$SegmentNo[SEQ.best$best_top25==1]
best.solution.top35<-SEQ.best$SegmentNo[SEQ.best$best_top35==1]

# sp representation of the best solution: 
sp.15<-species.distribution.df[match(best.solution.top15,species.distribution.df$SegNo),]
sp.25<-species.distribution.df[match(best.solution.top25,species.distribution.df$SegNo),]
sp.35<-species.distribution.df[match(best.solution.top35,species.distribution.df$SegNo),]

rep.15<-colSums(sp.15[,c(2:26)])/n.sp
rep.25<-colSums(sp.25[,c(2:26)])/n.sp
rep.35<-colSums(sp.35[,c(2:26)])/n.sp

rep.mean.15<-mean(colSums(sp.15[,c(2:26)])/n.sp)
rep.mean.25<-mean(colSums(sp.25[,c(2:26)])/n.sp)
rep.mean.35<-mean(colSums(sp.35[,c(2:26)])/n.sp)

refuge.size.systematic<-data.frame(top15=nrow(sp.15),top25=nrow(sp.25),top35=nrow(sp.35),method="Systematic refuges")
sp.rep.systematic<-data.frame(top15=rep.mean.15,sd15=sd(rep.15),top25=rep.mean.25,sd25=sd(rep.25),top35=rep.mean.35,sd35=sd(rep.35),method="Systematic refuges")

#---
# 2.4 plot mean sp rep for all three methods
#---
refuge.size<-rbind(refuge.size.water.only,refuge.size.positional,refuge.size.systematic)
sp.rep<-rbind(sp.rep.water.only,sp.rep.positional,sp.rep.systematic)

library(reshape)
refuge.size.melt<-melt(refuge.size,id="method")
sp.rep.mean.melt<-melt(sp.rep[,c(1,3,5,7)],id=c("method"))
sp.rep.sd.melt<-melt(sp.rep[,c(2,4,6,7)],id=c("method"))
sp.rep.mean.melt$sd<-sp.rep.sd.melt$value*100
sp.rep.mean.melt$value<-sp.rep.mean.melt$value*100
names(sp.rep.mean.melt)[1]<-"Method"

library(ggplot2)
ggplot(data = refuge.size.melt,aes(x=variable,y=value,group=method,color=method))+
  geom_line()+
  geom_point(aes(shape=method,size=4))+theme_classic()+
  xlab("Threshold / Target")+ylab("Number of stream segments")+
  scale_x_discrete(labels=c("top15% / 15%", "top25% / 25%","top35% / 35%"))+
  guides(size=FALSE)+
  labs(title=c("Priority refuge network size"))+
  theme(legend.position = "none")+
  ggsave(filename = "Figures/03_Objective function/Overal representation of each method_1.png",width = 4,height = 4)

ggplot(data = sp.rep.mean.melt,aes(x=variable,y=value,group=Method,color=Method))+
  geom_errorbar(aes(ymin=value-sd,ymax=value+sd),width=.1,position = position_dodge(0.05))+
  geom_line()+
  geom_point(aes(shape=Method,size=4))+theme_classic()+
  xlab("Threshold / Target")+ylab("Mean % of species total distribution")+
  labs(title=c("Mean species representation"))+guides(size=FALSE)+
  scale_x_discrete(labels=c("top15% / 15%", "top25% / 25%","top35% / 35%"))+
  ggsave(filename = "Figures/03_Objective function/Overal representation of each method_2.png",width = 6,height = 4)
