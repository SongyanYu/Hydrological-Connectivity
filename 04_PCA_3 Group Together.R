#-----------------------
# This script aims to use systematic methods to identify priority conservation areas
# Three mobility group speices are considered at the same time.
# Author: Songyan Yu
# Date create: 01/04/2019
#-----------------------

# require "PCA.water.only"

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)
plot(sdm)

SEQ.Clip<-readShapePoly(fn="Data/Shapfile/SEQ_Clip.shp")
plot(SEQ.Clip)
names(SEQ.Clip)
hierarchy<-data.frame(site=SEQ.Clip$SegmentNo,nextds=SEQ.Clip$DWNID1)

species.distribution.df<-sdm@data[,c(185:214)]
species.distribution.df<-species.distribution.df[,-c(14,15,18,20,21,26)]  # delete 6 non-selecte species

sp.mobility<-read.csv("Data/R data/Species mobility.csv")

PCA.water.only<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05")
names(PCA.water.only)
PCA.water.only$Freq

# randomly select a PCA1
# check whether the conservation target is met. if yes, break; if no, do next step.
# check the dwelling species that has the most limited mobility
# use that mobility to search for potentially recolonisable streams.
# remove all streams that are within the recolonisable extent
# then randomly selecte another PCA1
# repeat step 2~5 unilt the conservation target (5 planning unit for each species) is met.

low.mobility<-10
medium.mobility<-50
high.mobility<-100

dam.segment<-c(859803,853302,856476,874709)  # consider the blocking effect of dams

PCA1_SegNo<-PCA.water.only$SegmentNo[PCA.water.only$Freq==1]
reachable.PCA.L<-all_reachable_PCAs(PCA=PCA1_SegNo,mobility = low.mobility)
reachable.PCA.M<-all_reachable_PCAs(PCA=PCA1_SegNo,mobility = medium.mobility)
reachable.PCA.H<-all_reachable_PCAs(PCA=PCA1_SegNo,mobility = high.mobility)

# delete those segNo that have no speices dwelling.
species.distribution.df<-species.distribution.df[rowSums(species.distribution.df[,-1])>0,]
PCA1_SP_SegNo<-intersect(PCA1_SegNo,species.distribution.df$SegNo)

cons.target<-5

prot.species.df<-data.frame()
remaining.seg<-PCA1_SP_SegNo

while(sum(colSums(prot.species.df[,-1])>4)<23){
  
  # randomly select a PCA1_SP_SegNo
  n<-sample(1:length(remaining.seg),1)
  
  # check whether conservation target is met
  prot.species<-species.distribution.df[match(remaining.seg[n],species.distribution.df$SegNo),]
  prot.species.df<-rbind(prot.species.df,prot.species)
  
  cat(sum(colSums(prot.species.df[,-1])>4),"\n")
  
  # identify which species protection is still needed
  deficit.species<-colnames(prot.species.df[,-1])[colSums(prot.species.df[,-1])<5]
  col.deficit<-match(deficit.species,colnames(species.distribution.df[,-1]))
  
  # update the selection pool by deleting the segNo that does not contribute species feature.
  if(length(col.deficit)>1){
    update.SP.SegNo<-species.distribution.df$SegNo[rowSums(species.distribution.df[,(col.deficit+1)]==1)>0]
    reduce.SP.SegNo<-setdiff(species.distribution.df$SegNo,update.SP.SegNo)
  }
  if(length(col.deficit)==1){
    update.SP.SegNo<-species.distribution.df$SegNo[species.distribution.df[,(col.deficit+1)]==1]
    reduce.SP.SegNo<-setdiff(species.distribution.df$SegNo,update.SP.SegNo)
  }

  # update the selection pool by deleting the reachable PCA1
  species.selected<-colnames(species.distribution.df[,-1])[species.distribution.df[match(remaining.seg[n],species.distribution.df$SegNo),-1]==1]
  least.mobility<-min(sp.mobility$Mobility[match(substr(species.selected,1,6),sp.mobility$Abbrev)])
  if(least.mobility==low.mobility){
    reduce.PCA1.SegNo<-reachable.PCA.L[match(remaining.seg[n],PCA1_SegNo)]
  }
  if(least.mobility==medium.mobility){
    reduce.PCA1.SegNo<-reachable.PCA.M[match(remaining.seg[n],PCA1_SegNo)]
  }
  if(least.mobility==high.mobility){
    reduce.PCA1.SegNo<-reachable.PCA.H[match(remaining.seg[n],PCA1_SegNo)]
  }
  
  # update the selection pool formally
  remaining.seg<-setdiff(remaining.seg,union(reduce.SP.SegNo,reduce.PCA1.SegNo))
}







