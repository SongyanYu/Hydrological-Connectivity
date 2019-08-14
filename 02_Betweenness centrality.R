#---------------------------
# This script calcualte betweenness centrality for aquatic refugia prioritisation.
#---------------------------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity")

library(maptools)
SEQ.clip<-readShapePoly("Data/Shapfile/SEQ_Clip.shp")
names(SEQ.clip)

hierarchy<-data.frame(site=SEQ.clip$SegmentNo,nextds=SEQ.clip$DWNID1)

Brisbane.clip<-readShapePoly("Data/Shapfile/BrisbaneRiver_Clip.shp")
plot(Brisbane.clip)
Logan.clip<-readShapePoly("Data/Shapfile/LoganAlbertRiver_Clip.shp")
Maroochy.clip<-readShapePoly("Data/Shapfile/MaroochyRiver_Clip.shp")
Pine.clip<-readShapePoly("Data/Shapfile/PineRiver_Clip.shp")
Southcoast.clip<-readShapePoly("Data/Shapfile/SouthCoast_Clip.shp")

#-------------------
# load functions
#-------------------

#library(catchstats)
pairReaches <- function(hierarchy, site1, site2) {
  
  sites <- c(site1, site2)
  all.ds.sites <- list_all_downstream(hierarchy, sites)
  
  # check if the two sites can be linked within stream network.
  if(length(intersect(all.ds.sites[[1]],all.ds.sites[[2]]))>1){
    ds.sites <- unique(c(setdiff(all.ds.sites[[1]], all.ds.sites[[2]]), setdiff(all.ds.sites[[2]], all.ds.sites[[1]])))
    ds.sites <- unique(c(ds.sites, site1, site2))
    ds.sites <- setdiff(ds.sites,sites)
    return(ds.sites)
  }
  else{
    #cat("the two sites can be linked within stream network.\n")
    return(NA)
  }

}

list_all_downstream <- function(hierarchy, catchnames) {
  
  all.ds.sites <- vector("list", length(catchnames))
  
  for (i in 1:length(catchnames)) {
    ds.sites <- alldownstream(hierarchy, catchnames[i])
    all.ds.sites[[i]] <- ds.sites
  }
  return(all.ds.sites)
}

alldownstream <- function(hierarchy, catchname) {
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

allupstream <- function(hierarchy, catchname) {
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


betweenness.centrality<-function(SegNo,hierarcy){
  
  pathway.full<-list()
  for(n in 1:3){
    pathway<-lapply(1:length(SegNo),FUN = function(i) pairReaches(hierarchy = hierarchy,site1 = SegNo[n], site2 = SegNo[i]))
    pathway.full<-append(pathway.full,pathway)
    cat(n," out of ", length(SegNo),"\n")
  }
  
  unlist.pathway.full<-unlist(pathway.full)
  bwtweenness.centrality<-sapply(1:length(SegNo),FUN = function(i) sum(unlist.pathway.full==SegNo[i],na.rm=TRUE))
  BC.df<-data.frame(SegNo=SegNo,BC=bwtweenness.centrality)
  return(BC.df)
}


betweenness.centrality.cl<-function(SegNo,hierarcy){
  
  pathway.full<-list()
  
  pathway.full<-parLapply(cl,1:3,fun = function(n){
    pathway<-lapply(1:length(SegNo),FUN = function(i) pairReaches(hierarchy = hierarchy,site1 = SegNo[n], site2 = SegNo[i]))
    pathway.full<-append(pathway.full,pathway)
  })
  
  unlist.pathway.full<-unlist(pathway.full)
  bwtweenness.centrality<-sapply(1:length(SegNo),FUN = function(i) sum(unlist.pathway.full==SegNo[i],na.rm=TRUE))
  BC.df<-data.frame(SegNo=SegNo,BC=bwtweenness.centrality)
  return(BC.df)
}


#---------------------
# calcualting betweenness centrality for each individual river network has the following benefit
# significantly reduce the computional time
#---------------------

#SegNo<-Maroochy.clip$SegmentNo
#BC.mar<-betweenness.centrality(SegNo = SegNo,hierarcy = hierarchy)
#saveRDS(BC.mar,file = "BC.maroochy")
BC.mar<-readRDS("Data/R data/BC.maroochy")

#SegNo<-Southcoast.clip$SegmentNo
#BC.sth<-betweenness.centrality(SegNo,hierarchy)
#saveRDS(BC.sth,file="BC.southcoast")
BC.sth<-readRDS("Data/R data/BC.southcoast")

#SegNo<-Pine.clip$SegmentNo
#BC.pin<-betweenness.centrality(SegNo,hierarchy)
#saveRDS(BC.pin,file = "BC.pine")
BC.pin<-readRDS("Data/R data/BC.pine")

#SegNo<-Logan.clip$SegmentNo
#BC.log<-betweenness.centrality(SegNo,hierarchy)
#saveRDS(BC.log,file = "BC.logan")
BC.log<-readRDS("Data/R data/BC.logan")

#SegNo<-Brisbane.clip$SegmentNo
#BC.bne<-betweenness.centrality(SegNo,hierarchy)
BC.bne<-readRDS("Data/R data/Betweenness Centrality_BNE")

# To make BC for stream segments from different river networks comparable,
# BC values were normalised within each river network before they were compared.
BC.mar$BC.nor<-(BC.mar$BC-min(BC.mar$BC))/(max(BC.mar$BC)-min(BC.mar$BC))
BC.sth$BC.nor<-(BC.sth$BC-min(BC.sth$BC))/(max(BC.sth$BC)-min(BC.sth$BC))
BC.pin$BC.nor<-(BC.pin$BC-min(BC.pin$BC))/(max(BC.pin$BC)-min(BC.pin$BC))
BC.log$BC.nor<-(BC.log$BC-min(BC.log$BC))/(max(BC.log$BC)-min(BC.log$BC))
BC.bne$BC.nor<-(BC.bne$BC-min(BC.bne$BC))/(max(BC.bne$BC)-min(BC.bne$BC))

BC.SEQ<-rbind(BC.mar,BC.sth,BC.pin,BC.log,BC.bne)

#-------------
# plotting BC fdc
#-------------
library(maptools)
SEQ.networks<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
names(SEQ.networks)
candidate.refuges<-SEQ.networks$SegmentNo[SEQ.networks$Freq>0]
candidate.BC<-data.frame(SegNo=candidate.refuges)

library(dplyr)
candidate.BC<-left_join(candidate.BC,BC.SEQ,by="SegNo")
sum(is.na(candidate.BC))  # should be "0"

library(hydroTSM)
width=17.47
ratio=0.8
ppi=30
png(filename = "Figures/02_Betweenness/BC_fdc.png",width = width*ppi,height = width*ratio*ppi)
fdc(candidate.BC$BC.nor,xlab = "Exceedance probability",ylab = "Normalised BC",
    lQ.thr = NA,hQ.thr = NA,thr.shw=FALSE,
    main=NULL,lwd = 3,pch = NA,cex.lab = 1.5)
dev.off()

#----
# prioritise BC
#----
top<-c(0.85,0.75,0.65)  # top 15%, 25% and 35%
BC.15<-quantile(candidate.BC$BC.nor,probs = top[1])
BC.25<-quantile(candidate.BC$BC.nor,probs = top[2])
BC.35<-quantile(candidate.BC$BC.nor,probs = top[3])

BC.15.seg<-candidate.BC$SegNo[candidate.BC$BC.nor>=BC.15]
BC.25.seg<-candidate.BC$SegNo[candidate.BC$BC.nor>=BC.25&candidate.BC$BC.nor<BC.15]
BC.35.seg<-candidate.BC$SegNo[candidate.BC$BC.nor>=BC.35&candidate.BC$BC.nor<BC.25]
BC.rest<-candidate.BC$SegNo[candidate.BC$BC.nor<BC.35]

SEQ.networks$BC_class<-4
SEQ.networks$BC_class[match(BC.15.seg,SEQ.networks$SegmentNo)]<-1
SEQ.networks$BC_class[match(BC.25.seg,SEQ.networks$SegmentNo)]<-2
SEQ.networks$BC_class[match(BC.35.seg,SEQ.networks$SegmentNo)]<-3
sum(SEQ.networks$BC_class==3)

names(SEQ.networks)
writeLinesShape(SEQ.networks,fn="Data/Shapfile/Betweenness centrality/SEQ_BC")


