#---------------------------
# This script calcualte betweenness centrality for aquatic refugia prioritisation.
#---------------------------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity")

library(maptools)
network<-readShapeLines("Data/Shapfile/SEQ_networks.shp")
plot(network)
names(network)

SEQ.clip<-readShapePoly("Data/Shapfile/SEQ_Clip.shp")
names(SEQ.clip)

Brisbane.clip<-readShapePoly("Data/Shapfile/BrisbaneRiver_Clip.shp")
plot(Brisbane.clip)
Logan.clip<-readShapePoly("Data/Shapfile/LoganAlbertRiver_Clip.shp")
Maroochy.clip<-readShapePoly("Data/Shapfile/MaroochyRiver_Clip.shp")
Pine.clip<-readShapePoly("Data/Shapfile/PineRiver_Clip.shp")
Southcoast.clip<-readShapePoly("Data/Shapfile/SouthCoast_Clip.shp")



#library(catchstats)
pairReaches <- function(hierarchy, site1, site2) {
  
  sites <- c(site1, site2)
  all.ds.sites <- list_all_downstream(hierarchy, sites)
  
  # check if the two sites can be linked within stream network.
  if(length(intersect(all.ds.sites[[1]],all.ds.sites[[2]]))>1){
    ds.sites <- unique(c(setdiff(all.ds.sites[[1]], all.ds.sites[[2]]), setdiff(all.ds.sites[[2]], all.ds.sites[[1]])))
    ds.sites <- unique(c(ds.sites, site1, site2))
    return(ds.sites)
  }
  else{
    cat("the two sites can be linked within stream network.\n")
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

SegNo<-Brisbane.clip$SegmentNo
SegNo<-Maroochy.clip$SegmentNo

pairReaches(hierarchy = hierarchy,site1 = 104838, site2 = 104931)

BC.bne<-lapply(1:length(SegNo),FUN = function(i) pairReaches(hierarchy = hierarchy.bne,site1 = SegNo[1], site2 = SegNo[i]))
sum(is.na(BC.bne))


BC.mar<-lapply(1:length(SegNo),FUN = function(i) pairReaches(hierarchy = hierarchy,site1 = SegNo[1], site2 = SegNo[i]))
sum(is.na(BC.mar))

