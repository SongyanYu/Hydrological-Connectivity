#-----------------------------
# This script aims to identify PCAs by taking species distribution with mobility (efficient vs. non-efficient).
# Essentially, it reduces the number of PCAs identified by naive mobility in the condition that all PCAs can be recolonised
# by species with mobility when flow resumes.
# Author: Songyan Yu
# Date Create: 27/12/2018
#-----------------------------

#---------------------
# Just read in
#---------------------
setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")
PCA.species.mobility.lst<-readRDS(file = "Data/R data/PCA_Varying_Mobility_Nonefficient")  # Non-efficient
PCA.species.mobility.efficient<-readRDS(file="Data/R data/PCA_Varing_Mobility_Efficient")  # Efficient

# Already read in...
# Next, combine it with river network for visulisation
sp.mobility<-read.csv("Data/R data/Species mobility.csv")

library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)

# For each species
PCA.species.mobility.df<-data.frame(SegNo=sdm$SEGMENTNO)

# Non-efficient PCA3s
overall.PCAs<-unique(unlist(PCA.species.mobility.lst))
PCA.species.mobility.df[match(overall.PCAs,PCA.species.mobility.df$SegNo),2]<-1   # Check whether '2' is still correct.
colnames(PCA.species.mobility.df)[2]<-"PCA_NoEFI"
PCA.species.mobility.df$PCA_NoEFI[is.na(PCA.species.mobility.df$PCA_NoEFI)]<-0

# Efficient PCA3s
PCA.species.mobility.df[match(PCA.species.mobility.efficient,PCA.species.mobility.df$SegNo),3]<-1
colnames(PCA.species.mobility.df)[3]<-"PCA_EFI"
PCA.species.mobility.df$PCA_EFI[is.na(PCA.species.mobility.df$PCA_EFI)]<-0

# Write the results down
#sdm@data<-sdm@data[,-c(220:243)]
sdm@data<-cbind(sdm@data,PCA.species.mobility.df[,c(2,3)])
names(sdm)
#writeLinesShape(sdm,"Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")

# Identify reachable PCAs for visulisation
# Take the 1st fish species for example. "AmbAga"
SEQ.Clip<-readShapePoly(fn="Data/Shapfile/SEQ_Clip.shp")
plot(SEQ.Clip)
names(SEQ.Clip)
hierarchy<-data.frame(site=SEQ.Clip$SegmentNo,nextds=SEQ.Clip$DWNID1)

SegNo.species<-PCA.species.mobility.lst[[1]]
max.mobility<-50

stream.mobility.lst<-list()
reachable.stream.lst<-list()
for(n in 1:length(SegNo.species)){             # Particular PCA
  SegNo<-SegNo.species[n]                                      # n
  
  # Upstream
  upstream.mobility<-Ups_Mobility(SegNo = SegNo,max.mobility = max.mobility)
  
  # Downstream
  downstream<-alldownstream(hierarchy = hierarchy,catchname = SegNo)
  downstream.mobility.1<-downstream[SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(downstream,SEQ.Clip$SegmentNo)]<=max.mobility]
  downstream.mobility.1<-downstream.mobility.1[!is.na(downstream.mobility.1)]  # remove NA from the vector
  
  mobility.down<-c()
  for(i in 2:length(downstream.mobility.1)){
    
    other.branch<-setdiff(hierarchy$site[downstream.mobility.1[i]==hierarchy$nextds],downstream.mobility.1)   # begin with 2
    max.mobility.temp<-max.mobility-(SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(downstream.mobility.1[i],SEQ.Clip$SegmentNo)])
    
    if(length(other.branch)==1){
      mobility.up.down<-Ups_Mobility(SegNo = other.branch,max.mobility = max.mobility.temp)
      mobility.down<-append(mobility.down,mobility.up.down)
    }
    
    if(length(other.branch)>1){
      for(j in 1:length(other.branch)){
        mobility.up.down<-Ups_Mobility(SegNo = other.branch[j],max.mobility = max.mobility.temp)
        mobility.down<-append(mobility.down,mobility.up.down)
      }
    }
  }
  stream.mobility<-c(upstream.mobility,mobility.down,downstream.mobility.1)
  stream.mobility.lst[[n]]<-stream.mobility
  
  
  cat(n," out of ",length(SegNo.species),"\n")
}

for(i in 1:length(stream.mobility.lst)){
  lst_1<-stream.mobility.lst[[i]]
  
  for(x in 1:length(lst_1)){
    cat("\"SegmentNo\"","=",lst_1[x],"OR",
        file = paste0("Ups-",lst_1[1],".txt"),
        append = TRUE)
  }
}


#---------------------
# OR, calculate it... then do the combination mentioned above
#---------------------

#------------------
# Efficient method to calcualte PCAs for varying mobility
#------------------

# prerequisite variables 
# 1) PCA.species.disc (running "03_PCA_Species_No mobility.R")
PCA.species.disc<-PCA.species.disc[,c(1:13,16:17,19,22:25,27:30)]  # Only 23 species are selected.
PCA.species.disc<-PCA.species.disc[,c(1:9,11,10,12:24)]   # Make sure the col order is consistent with other datasets.

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)

SEQ.Clip<-readShapePoly(fn="Data/Shapfile/SEQ_Clip.shp")
plot(SEQ.Clip)
names(SEQ.Clip)
hierarchy<-data.frame(site=SEQ.Clip$SegmentNo,nextds=SEQ.Clip$DWNID1)

sp.mobility<-read.csv("Data/R data/Species mobility.csv")

dam.segment<-c(859803,853302,856476,874709)  # consider the blocking effect of dams

species.L<-sp.mobility$Abbrev[sp.mobility$Mobility==10]   # species in Low mobility group
n.L<-match(species.L,colnames(PCA.species.disc))  
species.M<-sp.mobility$Abbrev[sp.mobility$Mobility==50]   # Midium mobility group
n.M<-match(species.M,colnames(PCA.species.disc))
species.H<-sp.mobility$Abbrev[sp.mobility$Mobility==100]  # High mobility group
n.H<-match(species.H,colnames(PCA.species.disc))


SegNo.L<-PCA.species.disc$SegNo[rowSums(PCA.species.disc[,n.L])>0]  # PCA2 for species in L group
SegNo.M<-PCA.species.disc$SegNo[rowSums(PCA.species.disc[,n.M])>0]  # PCA2 for species in M group
SegNo.H<-PCA.species.disc$SegNo[rowSums(PCA.species.disc[,n.H])>0]  # PCA2 for species in H group

low.mobility<-10  # starting from the Low mobility group
medium.mobility<-50
high.mobility<-100

reachable.PCA.L<-all_reachable_PCAs(PCA=SegNo.L,mobility = low.mobility)

num<-0
m<-1
n.PCAs<-c()
routine.lst<-list()
for(x in 1:2000){
  # PCA.species.L
  remaining.seg<-SegNo.L
  routine<-c()
  
  while(length(remaining.seg)>0){
    n.lst<-match(remaining.seg[sample(1:length(remaining.seg),1)],SegNo.L)
    remaining.seg.next<-setdiff(remaining.seg,reachable.PCA.L[[n.lst]])
    
    if(length(remaining.seg)>length(remaining.seg.next)){
      routine<-append(routine,SegNo.L[n.lst])
      remaining.seg<-remaining.seg.next
    }
  }
  
  # whether the conservation target of 5 PCA for each species is met.
  target.meet<-c()
  for(i in 2:ncol(PCA.species.disc)){
    target.meet[i-1]<-sum(PCA.species.disc$SegNo[PCA.species.disc[,i]>0] %in% routine)
  }
  
  if(sum(target.meet>4)==23){
  # whether PCA.species.L can make all PCA2s.M connected
  stream.mobility.M<-all.reachable.streams(PCA = routine, mobility = medium.mobility)
  reachable.stream.M<-unique(unlist(stream.mobility.M))  # all stream segments that can be reached from PCA3.L when mobility is changed to be 50km.
  
  if(sum(SegNo.M %in% reachable.stream.M)==length(SegNo.M)){
    
    if(sum(SegNo.H %in% reachable.stream.M)==length(SegNo.M)){
      solution<-routine
      num<-num+1
    }
    else{
      stream.mobility.H<-all.reachable.streams(PCA = routine, mobility = high.mobility)
      reachable.stream.H<-unique(unlist(stream.mobility.H))  # all stream segments that can be reached from PCA3.L when mobility is changed to be 50km.
      
      # whether PCA.species.L can make all PCA2s.H connected
      if(sum(SegNo.H %in% reachable.stream.H)==length(SegNo.H)){
        solution<-routine
        num<-num+1
      }
      else{
        extra.routine<-add_new_PCA(SegNo = SegNo.H,reachable.stream = reachable.stream.H,mobility = high.mobility)
        solution<-c(routine,extra.routine)
      }
    }
  }
  
  else{
    extra.routine<-add_new_PCA(SegNo = SegNo.M,reachable.stream = reachable.stream.M, mobility = medium.mobility)
    solution<-c(routine,extra.routine)
    
    stream.mobility.H<-all.reachable.streams(PCA = solution, mobility = high.mobility)
    reachable.stream.H<-unique(unlist(stream.mobility.H))  # all stream segments that can be reached from PCA3.L when mobility is changed to be 50km.
    
    if(sum(SegNo.H %in% reachable.stream.H)==length(SegNo.H)){
      solution<-solution
    }
    else{
      extra.routine<-add_new_PCA(SegNo = SegNo.H,reachable.stream = reachable.stream.H,mobility = high.mobility)
      solution<-c(solution,extra.routine)
    }
  }

  routine.lst[[m]]<-solution
  n.PCAs[m]<-length(solution)
  m<-m+1
  }
  
  #plot(table(n.PCAs))
  cat(x," out of ",900,"\n")
}

saveRDS(object = routine.lst,file = "Data/R data/2000_solutions")
extra.solutions<-readRDS(file = "Data/R data/99_solutions")

routine.lst<-append(routine.lst,extra.solutions)

# conservation target of 5 PCAs for each species
#selected.routine.lst<-list()
#m<-1
#for(j in 1:length(routine.lst)){
#  target.meet<-c()
#  for(i in 2:ncol(PCA.species.disc)){
#    target.meet[i-1]<-sum(PCA.species.disc$SegNo[PCA.species.disc[,i]>0] %in% routine.lst[[j]])
#  }
#  if(sum(target.meet>4)==23){
#    selected.routine.lst[[m]]<-routine.lst[[j]]
#    m<-m+1
#  }
#  cat(j, "out of 9999",'\n')
#}

#n.PCAs<-lengths(selected.routine.lst)
#plot(n.PCAs)
#saveRDS(object = selected.routine.lst,file = "Data/R data/PCA3_")

# Selection frequency
frequency.PCA2<-data.frame(table(unlist(routine.lst)))
colnames(frequency.PCA2)[1]<-"SegNo"
frequency.PCA2$Per<-frequency.PCA2$Freq/length(routine.lst)
frequency.PCA2$class[frequency.PCA2$Per>=0.75]<-1
frequency.PCA2$class[frequency.PCA2$Per>=0.5&frequency.PCA2$Per<0.75]<-2
frequency.PCA2$class[frequency.PCA2$Per>=0.25&frequency.PCA2$Per<0.5]<-3
frequency.PCA2$class[frequency.PCA2$Per<0.25]<-4

summary(frequency.PCA2$Per)
plot(frequency.PCA2$Per)

# mobility = 10km
total.PRL<-sapply(routine.lst[which(n.PCAs==min(n.PCAs))],FUN = function(x) sum(sdm$HydCon_10[match(x,sdm$SEGMENTNO)]))  

# optimal PCA3s for all species 
PCA3.best<-unlist(routine.lst[which(n.PCAs==min(n.PCAs))[match(max(total.PRL),total.PRL)]])
PCA3.best.df<-data.frame(SegNo=sdm$SegNo)
PCA3.best.df$PCA3_best[PCA3.best.df$SegNo %in% PCA3.best]<-1
PCA3.best.df$PCA3_best[!(PCA3.best.df$SegNo %in% PCA3.best)]<-0

# combinw with sdm
frequency.PCA2$SegNo<-as.numeric(as.character(frequency.PCA2$SegNo))
sdm@data<-left_join(sdm@data,frequency.PCA2,by=("SegNo"))
sdm@data<-cbind(sdm@data,PCA3.best.df$PCA3_best)
names(sdm@data)[224]<-"PCA3_best"
names(sdm@data)

writeLinesShape(x=sdm,fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")

#-------------------
# non-efficient method (for comparison)
#-------------------
# prerequisite variables 
# 1) PCA.species.disc (running "02_PCA_Species_No mobility.R")
#PCA.species.disc<-PCA.species.disc[,c(1:13,16:17,19,22:25,27:30)]  # Only 23 species are selected.
#PCA.species.disc<-PCA.species.disc[,c(1:9,11,10,12:24)]   # Make sure the col order is consistent with other datasets.

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)

SEQ.Clip<-readShapePoly(fn="Data/Shapfile/SEQ_Clip.shp")
plot(SEQ.Clip)
names(SEQ.Clip)
hierarchy<-data.frame(site=SEQ.Clip$SegmentNo,nextds=SEQ.Clip$DWNID1)

sp.mobility<-read.csv("Data/R data/Species mobility.csv")

dam.segment<-c(859803,853302,856476,874709)  # consider the blocking effect of dams
PCA.species.mobility.lst<-list()

for(m in 2:ncol(PCA.species.disc)){                # particular species
  SegNo.species<-PCA.species.disc$SegNo[PCA.species.disc[,m]>0]      # m
  max.mobility<-sp.mobility[(m-1),2]                             # m
  
  stream.mobility.lst<-list()
  reachable.stream.lst<-list()
  for(n in 1:length(SegNo.species)){             # Particular PCA2
    SegNo<-SegNo.species[n]                                      # n
    
    # Upstream
    upstream.mobility<-Ups_Mobility(SegNo = SegNo,max.mobility = max.mobility)
    
    # Downstream
    downstream<-alldownstream(hierarchy = hierarchy,catchname = SegNo)
    downstream.mobility.1<-downstream[SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(downstream,SEQ.Clip$SegmentNo)]<=max.mobility]
    downstream.mobility.1<-downstream.mobility.1[!is.na(downstream.mobility.1)]  # remove NA from the vector
    
    mobility.down<-c()
    for(i in 2:length(downstream.mobility.1)){
      
      other.branch<-setdiff(hierarchy$site[downstream.mobility.1[i]==hierarchy$nextds],downstream.mobility.1)   # begin with 2
      max.mobility.temp<-max.mobility-(SEQ.Clip$D2OUTLET[match(SegNo,SEQ.Clip$SegmentNo)]-SEQ.Clip$D2OUTLET[match(downstream.mobility.1[i],SEQ.Clip$SegmentNo)])
      
      if(length(other.branch)==1){
        mobility.up.down<-Ups_Mobility(SegNo = other.branch,max.mobility = max.mobility.temp)
        mobility.down<-append(mobility.down,mobility.up.down)
      }
      
      if(length(other.branch)>1){
        for(j in 1:length(other.branch)){
          mobility.up.down<-Ups_Mobility(SegNo = other.branch[j],max.mobility = max.mobility.temp)
          mobility.down<-append(mobility.down,mobility.up.down)
        }
      }
    }
    
    stream.mobility<-c(upstream.mobility,mobility.down,downstream.mobility.1)
    stream.mobility.lst[[n]]<-stream.mobility
    reachable.streams<-base::intersect(stream.mobility,SegNo.species)
    reachable.stream.lst[[n]]<-reachable.streams
    
    cat(n," out of ",length(SegNo.species),"\n")
  }
  
  n.PCAs<-c()
  routine.lst<-list()
  for(x in 1:9999){
    remaining.seg<-SegNo.species
    routine<-c()
    
    while(length(remaining.seg)>0){
      n.lst<-match(remaining.seg[sample(1:length(remaining.seg),1)],SegNo.species)
      remaining.seg.next<-setdiff(remaining.seg,reachable.stream.lst[[n.lst]])
      
      if(length(remaining.seg)>length(remaining.seg.next)){
        routine<-append(routine,SegNo.species[n.lst])
        remaining.seg<-remaining.seg.next
      }
    }
    routine.lst[[x]]<-routine
    n.PCAs[x]<-length(routine)
    plot(table(n.PCAs))
    cat(x," out of ",9999,"\n")
  }
  
  if(max.mobility==50){
    total.PRL<-sapply(routine.lst[which(n.PCAs==min(n.PCAs))],FUN = function(x) sum(sdm$HydCon_50[match(x,sdm$SEGMENTNO)]))
  }
  if(max.mobility==10){
    total.PRL<-sapply(routine.lst[which(n.PCAs==min(n.PCAs))],FUN = function(x) sum(sdm$HydCon_10[match(x,sdm$SEGMENTNO)]))
  }
  if(max.mobility==100){
    total.PRL<-sapply(routine.lst[which(n.PCAs==min(n.PCAs))],FUN = function(x) sum(sdm$HydCon_100[match(x,sdm$SEGMENTNO)]))
  }
  
  PCA.species.mobility.lst[[m-1]]<-unlist(routine.lst[which(n.PCAs==min(n.PCAs))[match(max(total.PRL),total.PRL)]])
  
}

# save
#saveRDS(PCA.species.mobility.lst,file = "Data/R data/PCA_Varying_Mobility_Nonefficient")

#------------------ Check out the catchment of interest in ArcGIS------------------------
SegNo.species[348]

lst_1<-stream.mobility.lst[[348]]
lst_1<-reachable.stream.lst[[348]]
lst_1<-routine

for(x in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[x],"OR",
      file = paste0("Ups-",lst_1[1],".txt"),
      append = TRUE)
}


#--------------------------- Custermised functions-------------------------
Ups_Mobility<-function(SegNo,max.mobility){
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
  
  return(upstream.mobility)
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


