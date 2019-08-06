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

SegNo<-Maroochy.clip$SegmentNo
#BC.mar<-betweenness.centrality(SegNo = SegNo,hierarcy = hierarchy)
#saveRDS(BC.mar,file = "BC.maroochy")
BC.mar<-readRDS("Data/R data/BC.maroochy")

SegNo<-Southcoast.clip$SegmentNo
#BC.sth<-betweenness.centrality(SegNo,hierarchy)
#saveRDS(BC.sth,file="BC.southcoast")
BC.sth<-readRDS("Data/R data/BC.southcoast")

SegNo<-Pine.clip$SegmentNo
#BC.pin<-betweenness.centrality(SegNo,hierarchy)
#saveRDS(BC.pin,file = "BC.pine")
BC.pin<-readRDS("Data/R data/BC.pine")

SegNo<-Logan.clip$SegmentNo
#BC.log<-betweenness.centrality(SegNo,hierarchy)
#saveRDS(BC.log,file = "BC.logan")
BC.log<-readRDS("Data/R data/BC.logan")

SegNo<-Brisbane.clip$SegmentNo
#BC.bne<-betweenness.centrality(SegNo,hierarchy)
BC.bne<-readRDS("Data/R data/Betweenness Centrality_BNE")

top.percent<-0.95  # 5%

priority.mar<-BC.mar$SegNo[BC.mar$BC>=quantile(BC.mar$BC,probs = top.percent)]
priority.sth<-BC.sth$SegNo[BC.sth$BC>=quantile(BC.sth$BC,probs = top.percent)]
priority.pin<-BC.pin$SegNo[BC.pin$BC>=quantile(BC.pin$BC,probs = top.percent)]
priority.log<-BC.log$SegNo[BC.log$BC>=quantile(BC.log$BC,probs = top.percent)]
priority.bne<-BC.bne$SegNo[BC.bne$BC>=quantile(BC.bne$BC,probs = top.percent)]

