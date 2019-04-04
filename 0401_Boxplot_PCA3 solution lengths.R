#-------------------------
# This script aims to boxplot the lengths of PCA3 solutions for the 5 incremental conservation targets
# Author: Songyan Yu
# Date create: 05/04/2019
#-------------------------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

lgt.PCA3.mob<-readRDS("Data/R data/Lengths_PCA3_Solutions Mob")
lgt.PCA3.nonmob<-readRDS("Data/R data/Lengths_PCA3_Solutions non Mob")

lgt.mob.df<-data.frame(lgt.PCA3.mob)
lgt.nonmob.df<-data.frame(lgt.PCA3.nonmob)

colnames(lgt.mob.df)<-c("5%","10%","15%","20%","25%")
colnames(lgt.nonmob.df)<-c("5%","10%","15%","20%","25%")

lgt.mob.df$group<-"Mobility"
lgt.nonmob.df$group<-"Naive mobility"

library(reshape)
mob.melt<-melt(lgt.mob.df,id.vars = "group")
nonmon.melt<-melt(lgt.nonmob.df,id.vars = "group")
lgt.PCA3<-rbind(mob.melt,nonmon.melt)

library(ggplot2)
p<-ggplot(lgt.PCA3,aes(x=variable,y=value,fill=group))+geom_boxplot()
p<-p+labs(x="Conservation target",y="# priority conservation areas \n identified systematiclly")
p<-p+theme_classic()+theme(legend.title = element_blank(),legend.position = c(0.15,0.85))
p+ggsave(filename = "Figures/03_Mobility/Boxplot_lgt_PCA3 solutions.png",width = 5.5,height = 5.5*0.6)

axis.title = element_text(size=10)