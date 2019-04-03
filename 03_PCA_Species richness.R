#----------------------------
# Scenario 2
# This script aims to identify the species richness for each of PCA1
# Author: Songyan Yu
# Date create: 03/04/2019
#----------------------------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

# PCA.water only
PCA.water.only<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05")
names(PCA.water.only)
PCA.water.only$Freq

PCA1_SegNo<-PCA.water.only$SegmentNo[PCA.water.only$Freq==1]

# Species distribution
library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)
species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,15,18,20,21,26)]  # delete 6 non-selecte species

# delete those segNo that have no speices dwelling.
species.distribution.df<-species.distribution.df[rowSums(species.distribution.df[,-1])>0,]

# intersect PCA1 and SP
PCA1_SP_SegNo<-intersect(PCA1_SegNo,species.distribution.df$SegNo)

# calculate species richness for each PCA1
PCA1_SP_distribution<-species.distribution.df[match(PCA1_SP_SegNo,species.distribution.df$SegNo),]
SP.richness.df<-data.frame(SegNo=PCA1_SP_distribution$SegNo,SP.richness=rowSums(PCA1_SP_distribution[,-1]))
summary(SP.richness.df$SP.richness)
SP.richness.df$class<-5
SP.richness.df$class[SP.richness.df$SP.richness<5&SP.richness.df$SP.richness>=1]<-4
SP.richness.df$class[SP.richness.df$SP.richness<8&SP.richness.df$SP.richness>=5]<-3
SP.richness.df$class[SP.richness.df$SP.richness<11&SP.richness.df$SP.richness>=8]<-2
SP.richness.df$class[SP.richness.df$SP.richness<17&SP.richness.df$SP.richness>=11]<-1

# calculate # of refugia (PCA1) each species has
n.PCA1.SP<-colSums(PCA1_SP_distribution[,-1])
write.csv(n.PCA1.SP,file = "Anal/number of PCA1 per species.csv")

# write the species richness down to river network for visulisation
sdm@data$SP_richness<-0
sdm@data$SP_class<-5
sdm$SP_richness[match(SP.richness.df$SegNo,sdm$SegNo)]<-SP.richness.df$SP.richness
sdm$SP_class[match(SP.richness.df$SegNo,sdm$SegNo)]<-SP.richness.df$class
head(sdm)
writeLinesShape(sdm,fn="Data/Shapfile/SP richness/SP richness")

