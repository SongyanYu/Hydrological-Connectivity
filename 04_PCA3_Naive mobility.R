#----------------------
# Scenario 3
# This script aims to identify PCA3 without mobility (3 mobility group considered at the same time)
# Author: Songyan Yu
# Date create: 04/04/2019
#----------------------

# require "n.PCA1.SP" from "03_PCA_Species richness.R"

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)

SEQ.Clip<-readShapePoly(fn="Data/Shapfile/SEQ_Clip.shp")
names(SEQ.Clip)
hierarchy<-data.frame(site=SEQ.Clip$SegmentNo,nextds=SEQ.Clip$DWNID1)

species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,15,18,20,21,26)]  # delete 6 non-selecte species

sp.mobility<-read.csv("Data/R data/Species mobility.csv")

PCA.water.only<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05")
names(PCA.water.only)
PCA.water.only$Freq

dam.segment<-c(859803,853302,856476,874709)  # consider the blocking effect of dams

# inundated SegNo
innudt.shp<-readShapeLines("Data/Shapfile/Innudt_SegNo.shp")
inundt.SegNo<-innudt.shp$SegmentNo
delete.seg<-c(859398,859529,856156)
inundt.SegNo<-inundt.SegNo[-match(delete.seg,inundt.SegNo)]

PCA1_SegNo<-PCA.water.only$SegmentNo[PCA.water.only$Freq==1]

# delete those segNo that have no speices dwelling.
species.distribution.df<-species.distribution.df[rowSums(species.distribution.df[,-1])>0,]
PCA1_SP_SegNo<-intersect(PCA1_SegNo,species.distribution.df$SegNo)
PCA1_SP_SegNo<-setdiff(PCA1_SP_SegNo,inundt.SegNo)  # also exclude inundated SegNo

scaling.factor<-c(0.05,0.1,0.15,0.2,0.25)

frequency.PCA3.df<-data.frame(No=c(1:1108))  # 1108 is the number of PCA1 after excluding the inundated and non-dwelling segments.
best.PCA3.lst<-list()
lgt.PCA3.lst<-list()

for(m in 1:length(scaling.factor)){
  cons.target<-floor(n.PCA1.SP*scaling.factor[m])
  PCA3.lst<-list()
  
  for(i in 1:1000){
    prot.species.df<-data.frame()
    remaining.seg<-PCA1_SP_SegNo
    
    while(sum(colSums(prot.species.df[,-1])>=cons.target)<23){
      
      # randomly select a PCA1_SP_SegNo
      n<-sample(1:length(remaining.seg),1)
      
      # check whether conservation target is met
      prot.species<-species.distribution.df[match(remaining.seg[n],species.distribution.df$SegNo),]
      prot.species.df<-rbind(prot.species.df,prot.species)
      
      #cat(sum(colSums(prot.species.df[,-1])>(cons.target-1)),"\n")
      
      # identify which species protection is still needed
      deficit.species<-colnames(prot.species.df[,-1])[colSums(prot.species.df[,-1])<cons.target]
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
      
      # update the selection pool finally
      remaining.seg<-setdiff(remaining.seg,reduce.SP.SegNo)
      
      if(length(remaining.seg)==0){
        break
      }
    }
    
    PCA3.lst[[i]]<-prot.species.df$SegNo
    cat(i,"\n")
    
  }
  
  lgt.PCA3.lst[[m]]<-lengths(PCA3.lst)
  
  # irreplacibility (selection frequency)
  frequency.PCA3<-data.frame(table(unlist(PCA3.lst)))
  colnames(frequency.PCA3)[1]<-"SegNo"
  frequency.PCA3$Per<-frequency.PCA3$Freq/length(PCA3.lst)
  frequency.PCA3$class[frequency.PCA3$Per>=0.75]<-1
  frequency.PCA3$class[frequency.PCA3$Per>=0.5&frequency.PCA3$Per<0.75]<-2
  frequency.PCA3$class[frequency.PCA3$Per>=0.25&frequency.PCA3$Per<0.5]<-3
  frequency.PCA3$class[frequency.PCA3$Per<0.25]<-4
  
  frequency.PCA3.df<-cbind(frequency.PCA3.df,frequency.PCA3)
  
  # write the PCA3 best solution and irreplability into the shapfile
  #summary(lengths(PCA3.lst))
  n.best.PCA3<-which(lengths(PCA3.lst)==min(lengths(PCA3.lst)))
  
  if(length(n.best.PCA3)==1){
    best.PCA3<-PCA3.lst[[n.best.PCA3]]
  }
  # use RDI (a measure of human impact) to choose the best PCA3
  if(length(n.best.PCA3)>1){
    n.low.RDI<-c()
    for(r in 1:length(n.best.PCA3)){
      n.low.RDI[r]<-sum(SEQ.Clip$RDI[match(PCA3.lst[[n.best.PCA3[r]]],SEQ.Clip$SegmentNo)]<0.5)
    }
    best.PCA3<-PCA3.lst[[n.best.PCA3[match(max(n.low.RDI),n.low.RDI)]]]
  }
  
  best.PCA3.lst[[m]]<-best.PCA3
  #temp<-frequency.PCA3$Per[match(best.PCA3,frequency.PCA3$SegNo)]
  
}

saveRDS(best.PCA3.lst,file="Data/R data/Best_PCA3 non_Mob")
saveRDS(frequency.PCA3.df,file = "Data/R data/Frequency_PCA3 non_Mob")
saveRDS(lgt.PCA3.lst,file = "Data/R data/Lengths_PCA3_Solutions non Mob")

frequency.PCA3.df<-frequency.PCA3.df[,c(2,5,9,13,17,21)]
colnames(frequency.PCA3.df)[2:6]<-c("Target_5%","Target_10%","Target_15%","Target_20%","Target_25%")

#names(sdm)
#sdm@data<-sdm@data[,-c(216:227)]

library(dplyr)
frequency.PCA3.df$SegNo<-as.numeric(as.character(frequency.PCA3.df$SegNo))
sdm@data<-left_join(sdm@data,frequency.PCA3.df,by="SegNo")
sdm@data[,216:220][is.na(sdm@data[,216:220])]<-5

writeLinesShape(sdm,fn="Data/Shapfile/PCA3 Mob/PCA3_Mob")


#----------------FUNCTIONS---------------------
allupstream <- function(hierarchy, catchname) {
  
  names(hierarchy) <- c("site", "nextds")
  
  if (length(which(hierarchy$site == catchname)) > 0) {
    catchname <- as.vector(catchname)
    allsc <- as.vector(hierarchy$site[hierarchy$nextds == catchname])
    allsc <- allsc[!is.na(allsc)]
    # subcatchments immediately upstream
    nbrnch <- end <- length(allsc)
    # number of branches immediately upstream
    start <- 1
    while (nbrnch > 0) {
      for (i in start:end) {
        allsc <- c(allsc, as.vector(hierarchy$site[hierarchy$nextds == allsc[i]]))
        allsc <- allsc[!is.na(allsc)]
      }
      start <- end + 1
      end <- length(allsc)
      nbrnch <- end - (start - 1)
    }
    allsc <- c(catchname, allsc)
    allsc
  } else cat(paste(catchname, "is not a site listed in the hierarchy table", "\n"))
}

alldownstream <- function(hierarchy, catchname) {
  
  names(hierarchy) <- c("site", "nextds")
  
  if (length(which(hierarchy$site == catchname)) > 0) {
    catchname <- as.vector(catchname)
    allsc <- as.vector(hierarchy$nextds[hierarchy$site == catchname])
    allsc <- allsc[!is.na(allsc)]
    # subcatchments immediately upstream
    nbrnch <- end <- length(allsc)
    # number of branches immediately upstream
    start <- 1
    while (nbrnch > 0 & !-1 %in% allsc) {
      for (j in start:end) {
        allsc <- c(allsc, as.vector(hierarchy$nextds[hierarchy$site == allsc[j]]))
        allsc <- allsc[!is.na(allsc)]
      }
      start <- end + 1
      end <- length(allsc)
      nbrnch <- end - (start - 1)
    }
    allsc <- c(catchname, allsc)
    allsc <- allsc[allsc != -1]
    allsc
  } else cat(paste(catchname, "is not a site listed in the hierarchy table", "\n"))
}

list_all_upstream <- function(hierarchy, catchnames) {
  
  all.us.sites <- vector("list", length(catchnames))
  
  for (i in 1:length(catchnames)) {
    us.sites <- allupstream(hierarchy, catchnames[i])
    all.us.sites[[i]] <- us.sites
  }
  return(all.us.sites)
}






