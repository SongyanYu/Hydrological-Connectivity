#-------------
# ggplot species representaiton across all three methods
# water-only refuges, positional refuges and systematic refuges
# Author: Songyan Yu
# Date create: 16/08/2019
#-------------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

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
  scale_fill_grey(start = 1,end = 0.3,labels=c("","Top15% / 15%","Top25% / 25%","Top35% / 35%"))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust=1))+
  xlab("Fish species")+ylab("% of Total Species Distribution")+labs(fill=c("Threshold / Target"))+
  ggsave(filename = "Figures/03_Objective function/Cumulative representation of each method.png")
