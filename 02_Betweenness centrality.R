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

hierarchy<-data.frame(site=SEQ.clip$SegmentNo,nextds=SEQ.clip$DWNID1)

#library(catchstats)
pairReaches <- function(hierarchy, site1, site2) {
  sites <- c(site1, site2)
  all.ds.sites <- list_all_downstream(hierarchy, sites)
  ds.sites <- unique(c(setdiff(all.ds.sites[[1]], all.ds.sites[[2]]), setdiff(all.ds.sites[[2]], all.ds.sites[[1]])))
  ds.sites <- unique(c(ds.sites, site1, site2))
  return(ds.sites)
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

pairReaches(hierarchy = hierarchy,site1 = SEQ.clip$SegmentNo[250], site2 = SEQ.clip$SegmentNo[260])

