#---------------------------
# This script calcualte betweenness centrality for aquatic refugia prioritisation.
# This is a Parallel version for HPC
# Only for Brisbane River catchment
#---------------------------

library(maptools)
SEQ.clip<-readShapePoly("/scratch/s2974665/Shapefile/SEQ_Clip.shp")
hierarchy<-data.frame(site=SEQ.clip$SegmentNo,nextds=SEQ.clip$DWNID1)

Brisbane.clip<-readShapePoly("/scratch/s2974665/Shapefile/BrisbaneRiver_Clip.shp")
SegNo<-Brisbane.clip$SegmentNo

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

betweenness.centrality.cl<-function(SegNo,hierarcy){
  
  pathway.full<-list()
  
  pathway.full<-parLapply(cl,1:length(SegNo),fun = function(n){
    pathway<-lapply(1:length(SegNo),FUN = function(i) pairReaches(hierarchy = hierarchy,site1 = SegNo[n], site2 = SegNo[i]))
    pathway.full<-append(pathway.full,pathway)
  })
  
  unlist.pathway.full<-unlist(pathway.full)
  bwtweenness.centrality<-sapply(1:length(SegNo),FUN = function(i) sum(unlist.pathway.full==SegNo[i],na.rm=TRUE))
  BC.df<-data.frame(SegNo=SegNo,BC=bwtweenness.centrality)
  return(BC.df)
}


#---------------------
# calcualting betweenness centrality for the Brisbane River catchment only
# Scripts for calculating BC for the other four river catchments were in "02_Betweenness centrality.R"
#---------------------

library(parallel)
no_cores<-detectCores()
cl<-makeCluster(no_cores,type = "FORK")

#clusterExport(cl,varlist = c("SegNo","hierarchy","pairReaches","list_all_downstream","alldownstream"))  # for WINDOWS OS

BC.bne.cl<-betweenness.centrality.cl(SegNo = SegNo,hierarcy = hierarchy)
stopCluster(cl)

saveRDS(BC.bne.cl,file="Betweenness Centrality_BNE")



