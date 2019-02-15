#---------------------------
# This script aims to calculate potential recolonisable length based on species mobility capacity.
# Author: Songyan Yu
# Date Create: 24/12/2018
#---------------------------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

# 1.1 Species distribution modelled by Ross et al. 2016
library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)
head(sdm@data)
sdm_old<-readShapeLines(fn="Data/Shapfile/Species distribution model/Ensemble_forecasts.shp")

PCA<-readShapeLines(fn="Data/Shapfile/Threshold of quant 0.5/PCA_Water only_sw08_fp05.shp")
names(PCA)

SEQ.Clip<-readShapePoly(fn="Data/Shapfile/SEQ_Clip.shp")
plot(SEQ.Clip)
names(SEQ.Clip)
hierarchy<-data.frame(site=SEQ.Clip$SegmentNo,nextds=SEQ.Clip$DWNID1)

#max.mobility<-10
#max.mobility<-50
max.mobility<-100

mobility.segment<-c()   # Potentially re-colonisable length (PRL) for each stream segment

# Calculate PRL
# Consider the dam effect on  fish passage (only downstream, never upstream)
dam.segment<-c(859803,853302,856476,874709)

for(m in 1:nrow(sdm)){
  SegNo=sdm$SEGMENTNO[m]
  
  mobility.up<-Mobility_Up(SegNo = SegNo,max.mobility = max.mobility)
  
  downstream<-alldownstream(hierarchy = hierarchy,catchname = SegNo)
  # downstream in the main stem
  downstream.mobility.1<-downstream[SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(downstream,SEQ.Clip$SegmentNo)]<=max.mobility]
  downstream.mobility.1<-downstream.mobility.1[!is.na(downstream.mobility.1)]  # remove NA from the vector
  mobility.down.1<-sum(PCA$RCHLEN[match(downstream.mobility.1,PCA$SegmentNo)],na.rm = TRUE)
  
  mobility.down<-c()
  for(i in 2:length(downstream.mobility.1)){
    
    other.branch<-setdiff(hierarchy$site[downstream.mobility.1[i]==hierarchy$nextds],downstream.mobility.1)   # begin with 2
    max.mobility.temp<-max.mobility-(SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(downstream.mobility.1[i],SEQ.Clip$SegmentNo)])
    
    #cat(length(other.branch),"\n")
    
    if(length(other.branch)==1){
      mobility.up.down<-Mobility_Up(SegNo = other.branch,max.mobility = max.mobility.temp)
      mobility.down<-append(mobility.down,mobility.up.down)
    }
    
    if(length(other.branch)>1){
      for(j in 1:length(other.branch)){
        mobility.up.down<-Mobility_Up(SegNo = other.branch[j],max.mobility = max.mobility.temp)
        mobility.down<-append(mobility.down,mobility.up.down)
      }
    }
  }
  
  mobility.segment[m]<-sum(mobility.up,mobility.down,mobility.down.1)
  
  cat(m," out of ",nrow(sdm),"\n")
}

sum(sdm$SEGMENTNO==sdm_old$SEGMENTNO)==nrow(sdm)

sdm$HydCon_10<-mobility.segment
sdm$HydCon_50<-mobility.segment
sdm$HydCon_100<-mobility.segment

names(sdm@data)

#writeLinesShape(sdm,fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")


#--------------------------- Custermised functions-------------------------
Mobility_Up<-function(SegNo,max.mobility){
  upstream<-allupstream(hierarchy = hierarchy,catchname = SegNo)
  upstream.mobility<-upstream[SEQ.Clip$D2OUTLET[match(upstream,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]<=max.mobility]
  
  if(sum(dam.segment %in% upstream.mobility)){
    
    if(sum(dam.segment %in% upstream.mobility)==1){
      block.segment<-dam.segment[dam.segment %in% upstream.mobility]
      upstream.block<-allupstream(hierarchy = hierarchy,catchname = block.segment)
      upstream.mobility<-setdiff(upstream.mobility,upstream.block)
    }
    
    if(sum(dam.segment %in% upstream.mobility)>1){
      block.segments<-dam.segment[dam.segment %in% upstream.mobility]
      upstream.block<-unique(unlist(list_all_upstream(hierarchy = hierarchy,catchnames = block.segments)))
      upstream.mobility<-setdiff(upstream.mobility,upstream.block)
    }
  }
  
  mobility.up<-sum(PCA$RCHLEN[match(upstream.mobility,PCA$SegmentNo)],na.rm=TRUE)
  return(mobility.up)
}

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

list_all_upstream <- function(hierarchy, catchnames) {
  
  all.us.sites <- vector("list", length(catchnames))
  
  for (i in 1:length(catchnames)) {
    us.sites <- allupstream(hierarchy, catchnames[i])
    all.us.sites[[i]] <- us.sites
  }
  return(all.us.sites)
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
