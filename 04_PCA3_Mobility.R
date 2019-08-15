#-----------------------
# This script aims to use systematic methods to identify priority conservation areas
# Three mobility group speices are considered at the same time.
# Author: Songyan Yu
# Date create: 01/04/2019
#-----------------------

#-----
# required variable: 
# 1) candidate.df, obtained from running "04_Candidate df.R"
#------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity/")

#----
# when re-visit, skip the optimisation process and directly read in previous outcome.
#---
# read in species distribution data
library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)
species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,18,21,26)]  # delete 4 non-selecte species
n.sp<-colSums(species.distribution.df[,-1])

# read in 100 potential solutions for each conservation target (15%, 25% and 35%).
solution.lst<-readRDS("Data/R data/PCA3_solution_top15_non_Mob")
solution.lst<-readRDS("Data/R data/PCA3_solution_top25_non_Mob")
solution.lst<-readRDS("Data/R data/PCA3_solution_top35_non_Mob")

solution.lst<-readRDS("Data/R data/PCA3_solution_top15_Mob")
solution.lst<-readRDS("Data/R data/PCA3_solution_top25_Mob")
solution.lst<-readRDS("Data/R data/PCA3_solution_top35_Mob")

# size of priority refuge network
n.seg<-lengths(solution.lst)
summary(n.seg)

# species representation
scaling.factor<-c(15,25,35)
cons.target<-floor(n.sp*(scaling.factor[3]/100))  # change the index of the scaling.factor.
rep.sp<-sapply(solution.lst,FUN = function(x){prot.species.df<-species.distribution.df[match(x,species.distribution.df$SegNo),]; rep.sp<-sum(colSums(prot.species.df[,-1])>=cons.target)})
summary(rep.sp)

#---
# selection frequency
#---
frequency.seg<-data.frame(table(unlist(solution.lst)))
summary(frequency.seg$Freq/100)

seq.network<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
names(seq.network)
seq.network$"selec_freq"<-0
seq.network$"selec_freq"[match(frequency.seg$Var1,seq.network$SegmentNo)]<-frequency.seg$Freq
names(seq.network)

writeLinesShape(seq.network,fn="Data/Shapfile/PCA3 non Mob/Selec freq_top35")

#---
# best solution with the minimum objective function value
#---

para.a<-0.4   # position penalty weight
para.b<-0.5  # feature penalty weight

obj.func<-sapply(solution.lst,FUN = function(x){
  prot.species.df<-species.distribution.df[match(x,species.distribution.df$SegNo),]
  RDI<-candidate.df$RDI[match(prot.species.df$SegNo,candidate.df$SegNo)]
  BC<-candidate.df$BC.nor[match(prot.species.df$SegNo,candidate.df$SegNo)]
  feature<-cons.target-colSums(prot.species.df[,-1])
  feature<-sapply(feature,FUN = function(x) {ifelse(x<0,0,x)})
  objective.func<-sum(RDI)+sum(para.a*(1-BC))+sum(para.b*feature)
})
summary(obj.func)

best.solution<-solution.lst[[which(obj.func==min(obj.func))]]
length(best.solution)




library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)
species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,18,21,26)]  # delete 4 non-selecte species
n.sp<-colSums(species.distribution.df[,-1])

#SEQ.Clip<-readShapePoly(fn="Data/Shapfile/SEQ_Clip.shp")
#hierarchy<-data.frame(site=SEQ.Clip$SegmentNo,nextds=SEQ.Clip$DWNID1)
#sp.mobility<-read.csv("Data/R data/Species mobility.csv")

#low.mobility<-5
#medium.mobility<-10
#high.mobility<-15

candidate<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
names(candidate)

#dam.segment<-c(859803,853302,856476,874709)  # consider the blocking effect of dams

# inundated SegNo
innudt.shp<-readShapeLines("Data/Shapfile/Innudt_SegNo.shp")
inundt.SegNo<-innudt.shp$SegmentNo
delete.seg<-c(859398,859529,856156)
inundt.SegNo<-inundt.SegNo[-match(delete.seg,inundt.SegNo)]

candidate.seg<-candidate$SegmentNo[candidate$Freq>0]
candidate.seg<-setdiff(candidate.seg,inundt.SegNo)

#reachable.PCA.L<-all_reachable_PCAs(PCA=candidate.seg,mobility = low.mobility)
#reachable.PCA.M<-all_reachable_PCAs(PCA=candidate.seg,mobility = medium.mobility)
#reachable.PCA.H<-all_reachable_PCAs(PCA=candidate.seg,mobility = high.mobility)

# delete those segNo that have no speices dwelling.
#species.distribution.df<-species.distribution.df[rowSums(species.distribution.df[,-1])>0,]
#PCA1_SP_SegNo<-intersect(PCA1_SegNo,species.distribution.df$SegNo)
#PCA1_SP_SegNo<-setdiff(PCA1_SP_SegNo,inundt.SegNo)  # also exclude inundated SegNo

# calculate # of refugia (PCA1) each species has (i.e. species distribution)
#PCA1_SP_distribution<-species.distribution.df[match(PCA1_SP_SegNo,species.distribution.df$SegNo),]
#n.PCA1.SP<-colSums(PCA1_SP_distribution[,-1])

scaling.factor<-c(15,25,35)  # used to set conservation target (% of species distribution).

for(m in 1:length(scaling.factor)){
  
  #cons.target<-floor(n.PCA1.SP*scaling.factor)
  cons.target<-floor(n.sp*(scaling.factor[m]/100))
  
  para.a<-0.4   # position penalty weight
  para.b<-0.5  # feature penalty weight
  rep.sp<-c()
  n.seg<-c()
  solution.lst<-list()
  obj.func<-c()
  zero.remaining.pool<-0
  for(i in 1:100){
    
    remaining.seg<-candidate.seg
    #------
    # 1. start an initial configuration by randomly selecting one stream segment.
    #------
    n<-sample(1:length(remaining.seg),1)
    
    #------
    # 2. calculate the value of objective function
    #------
    prot.species.df<-species.distribution.df[match(remaining.seg[n],species.distribution.df$SegNo),]
    
    RDI<-candidate.df$RDI[match(prot.species.df$SegNo,candidate.df$SegNo)]
    BC<-candidate.df$BC.nor[match(prot.species.df$SegNo,candidate.df$SegNo)]
    feature<-cons.target-colSums(prot.species.df[,-1])
    feature<-sapply(feature,FUN = function(x) {ifelse(x<0,0,x)})
    
    objective.func<-sum(RDI)+sum(para.a*(1-BC))+sum(para.b*feature)
    
    #------
    # 3. remove candidate segments from remaining pool
    # 1) reachable segments from the retained (NA)
    # 2) segments not contributing to species feature
    #------
    
    # 1) reachable segments from the retained
    #species.selected<-colnames(species.distribution.df[,-1])[species.distribution.df[match(remaining.seg[n],species.distribution.df$SegNo),-1]==1]
    #least.mobility<-min(sp.mobility$Mobility[match(substr(species.selected,1,6),sp.mobility$Abbrev)])
    #if(least.mobility==low.mobility){
    #  reduce.PCA1.SegNo<-reachable.PCA.L[[match(remaining.seg[n],candidate.seg)]]
    #}
    #if(least.mobility==medium.mobility){
    #  reduce.PCA1.SegNo<-reachable.PCA.M[[match(remaining.seg[n],candidate.seg)]]
    #}
    #if(least.mobility==high.mobility){
    #  reduce.PCA1.SegNo<-reachable.PCA.H[[match(remaining.seg[n],candidate.seg)]]
    #}
    
    # 2) segments not contributing to species feature
    # identify which species protection is still needed
    #deficit.species<-colnames(prot.species.df[,-1])[colSums(prot.species.df[,-1])<cons.target]
    #col.deficit<-match(deficit.species,colnames(species.distribution.df[,-1]))
    #if(length(col.deficit)>1){
    #  update.SP.SegNo<-species.distribution.df$SegNo[rowSums(species.distribution.df[,(col.deficit+1)]==1)>0]
    #  reduce.SP.SegNo<-setdiff(species.distribution.df$SegNo,update.SP.SegNo)
    #}
    #if(length(col.deficit)==1){
    #  update.SP.SegNo<-species.distribution.df$SegNo[species.distribution.df[,(col.deficit+1)]==1]
    #  reduce.SP.SegNo<-setdiff(species.distribution.df$SegNo,update.SP.SegNo)
    #}
    
    # 3) update the selection pool formally
    #remaining.seg<-setdiff(remaining.seg,base::union(reduce.SP.SegNo,prot.species.df$SegNo))
    #remaining.seg<-setdiff(remaining.seg,prot.species.df$SegNo)
    #remaining.seg<-setdiff(remaining.seg,reduce.PCA1.SegNo)
    
    #------
    # 4. iteratively try a new segment to add to the current configuration
    # 1) randomly select a new segment
    # 2) check whether the new segment can improve/worsen the value of objective function
    # 3) if improving (reducing), keep it in; if worsening (increasing), remove it.
    #------
    iteration.n<-0
    #while(iteration.n<=1000&sum(colSums(prot.species.df[,-1])>=cons.target)<25){
    while(iteration.n<=1000){
      
      # 1) randomly select a new segment
      n<-sample(1:length(remaining.seg),1)
      
      # check whether the selected segment is in the current configuration.
      if(remaining.seg[n] %in% prot.species.df$SegNo){
        prot.species.df.2<-prot.species.df[-match(remaining.seg[n],prot.species.df$SegNo),]
        
        RDI<-candidate.df$RDI[match(prot.species.df.2$SegNo,candidate.df$SegNo)]
        BC<-candidate.df$BC.nor[match(prot.species.df.2$SegNo,candidate.df$SegNo)]
        feature<-cons.target-colSums(prot.species.df.2[,-1])
        feature<-sapply(feature,FUN = function(x) {ifelse(x<0,0,x)})
        
        objective.func.2<-sum(RDI)+sum(para.a*(1-BC))+sum(para.b*feature)  # need to adjust para.b
        
        if(objective.func>=objective.func.2){
          prot.species.df<-prot.species.df.2   # update current configuration
          objective.func<-objective.func.2    # update the value of objective function
          
          iteration.n<-0
          # update remaining seg
          #species.selected<-colnames(species.distribution.df[,-1])[species.distribution.df[match(remaining.seg[n],species.distribution.df$SegNo),-1]==1]
          #least.mobility<-min(sp.mobility$Mobility[match(substr(species.selected,1,6),sp.mobility$Abbrev)])
          #if(least.mobility==low.mobility){
          #  reduce.PCA1.SegNo<-reachable.PCA.L[[match(remaining.seg[n],candidate.seg)]]
          #}
          #if(least.mobility==medium.mobility){
          #  reduce.PCA1.SegNo<-reachable.PCA.M[[match(remaining.seg[n],candidate.seg)]]
          #}
          #if(least.mobility==high.mobility){
          #  reduce.PCA1.SegNo<-reachable.PCA.H[[match(remaining.seg[n],candidate.seg)]]
          #}    
          
          #remaining.seg<-c(remaining.seg,reduce.PCA1.SegNo)
        }
        else{
          #cat(iteration.n,"\n")
          iteration.n<-iteration.n+1
        }
      }
      else{
        prot.species<-species.distribution.df[match(remaining.seg[n],species.distribution.df$SegNo),]
        prot.species.df.2<-rbind(prot.species.df,prot.species)
        
        # 2) check whether the new segment can improve/wosen the value of objective function
        # calculate the value of objective function
        RDI<-candidate.df$RDI[match(prot.species.df.2$SegNo,candidate.df$SegNo)]
        BC<-candidate.df$BC.nor[match(prot.species.df.2$SegNo,candidate.df$SegNo)]
        feature<-cons.target-colSums(prot.species.df.2[,-1])
        feature<-sapply(feature,FUN = function(x) {ifelse(x<0,0,x)})
        
        objective.func.2<-sum(RDI)+sum(para.a*(1-BC))+sum(para.b*feature)  # need to adjust para.b
        
        if(objective.func>=objective.func.2){
          prot.species.df<-prot.species.df.2   # update current configuration
          objective.func<-objective.func.2    # update the value of objective function
          
          # update remaining pool
          # 1) reachable segments from the retained
          #species.selected<-colnames(species.distribution.df[,-1])[species.distribution.df[match(remaining.seg[n],species.distribution.df$SegNo),-1]==1]
          #least.mobility<-min(sp.mobility$Mobility[match(substr(species.selected,1,6),sp.mobility$Abbrev)])
          #if(least.mobility==low.mobility){
          #  reduce.PCA1.SegNo<-reachable.PCA.L[[match(remaining.seg[n],candidate.seg)]]
          #}
          #if(least.mobility==medium.mobility){
          #  reduce.PCA1.SegNo<-reachable.PCA.M[[match(remaining.seg[n],candidate.seg)]]
          #}
          #if(least.mobility==high.mobility){
          #  reduce.PCA1.SegNo<-reachable.PCA.H[[match(remaining.seg[n],candidate.seg)]]
          #}
          
          # 2) segments not contributing to species feature
          # identify which species protection is still needed
          #deficit.species<-colnames(prot.species.df[,-1])[colSums(prot.species.df[,-1])<cons.target]
          #col.deficit<-match(deficit.species,colnames(species.distribution.df[,-1]))
          #if(length(col.deficit)>1){
          #  update.SP.SegNo<-species.distribution.df$SegNo[rowSums(species.distribution.df[,(col.deficit+1)]==1)>0]
          #  reduce.SP.SegNo<-setdiff(species.distribution.df$SegNo,update.SP.SegNo)
          #}
          #if(length(col.deficit)==1){
          #  update.SP.SegNo<-species.distribution.df$SegNo[species.distribution.df[,(col.deficit+1)]==1]
          #  reduce.SP.SegNo<-setdiff(species.distribution.df$SegNo,update.SP.SegNo)
          #}
          
          # 3) update the selection pool formally
          #remaining.seg<-setdiff(remaining.seg,union(reduce.SP.SegNo,reduce.PCA1.SegNo))
          #remaining.seg<-setdiff(remaining.seg,base::union(reduce.SP.SegNo,prot.species.df$SegNo))
          #remaining.seg<-setdiff(remaining.seg,prot.species.df$SegNo)
          #remaining.seg<-setdiff(remaining.seg,reduce.PCA1.SegNo)
          
          
          if(length(remaining.seg)==0){
            zero.remaining.pool<-zero.remaining.pool+1
            break
          }
          iteration.n<-0
        }
        else{
          #cat(iteration.n,"\n")
          iteration.n<-iteration.n+1
        }
      }
      #cat(nrow(prot.species.df),"\n",file = "dynamics of selection_b015.txt",append = TRUE)
    }
    n.seg[i]<-nrow(prot.species.df)
    solution.lst[[i]]<-prot.species.df$SegNo
    rep.sp[i]<-sum(colSums(prot.species.df[,-1])>=cons.target)
    obj.func[i]<-objective.func
    
    cat("new loop ",i,"\n")
  }
  saveRDS(solution.lst,file = "Data/R data/PCA3_solution_top15_non_Mob_v3")
  saveRDS(n.seg,file = "Data/R data/PCA3_size_top15_non_Mob_v3")
  saveRDS(rep.sp,file = "Data/R data/PCA3_repSp_top15_non_Mob_v3")
  saveRDS(obj.func,file = "Data/R data/PCA3_objFunc_top15_non_Mob_v3")
}
#saveRDS(solution.lst,file = paste0("Data/R data/PCA3_solution_top",scaling.factor[m],"_non_Mob_v2"))

zero.remaining.pool
frequency.seg<-data.frame(table(unlist(solution.lst)))
summary(frequency.seg$Freq/100)
summary(n.seg)
summary(rep.sp)




# sp representation of the best solution
plot(colSums(prot.species.df[,c(2:26)])/n.sp)
mean(colSums(prot.species.df[,c(2:26)])/colSums(species.distribution.df[,c(2:26)]))

# sp rep for randomly selected segments in the same number as that systematically selected
random.seg<-PCA1_SegNo[runif(357,min = 1,max=length(PCA1_SegNo))]
plot(colSums(species.distribution.df[match(random.seg,species.distribution.df$SegNo),-1])/n.sp)



# near-optimal solutions
summary(objective.func)
up.quant<-which(objective.func<=quantile(objective.func,probs = 0.25))
frequency.PCA3<-data.frame(table(unlist(PCA3.lst[up.quant])))
is.factor(frequency.PCA3$Var1)
frequency.PCA3$Var1=as.numeric(as.character(frequency.PCA3$Var1))
colnames(frequency.PCA3)<-c("SegNo","SlcFrq")

frequency.PCA3$class[frequency.PCA3$SlcFrq<=15]<-4
frequency.PCA3$class[frequency.PCA3$SlcFrq>15&frequency.PCA3$SlcFrq<=58]<-3
frequency.PCA3$class[frequency.PCA3$SlcFrq>58&frequency.PCA3$SlcFrq<=150]<-2
frequency.PCA3$class[frequency.PCA3$SlcFrq>150]<-1

sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)
sdm@data<-sdm@data[,-c(216:227)]
names(sdm)
sdm@data<-left_join(sdm@data,frequency.PCA3,by=c("SegNo"))

# the best solution
best.slt<-unlist(PCA3.lst[objective.func==min(objective.func)])
sdm$best.slt<-0
sdm$best.slt[match(best.slt,sdm$SegNo)]<-1
length(sdm$SegNo)

writeLinesShape(sdm,fn="Data/Shapfile/PCA3_objective func")

# identify reachable streams from best solution
low.reachable.streams.lst<-all.reachable.streams(best.slt,mobility = 10)
medium.reachable.streams.lst<-all.reachable.streams(best.slt,mobility = 50)
high.reachable.streams.lst<-all.reachable.streams(best.slt,mobility = 100)

low.reachable.streams<-unique(unlist(low.reachable.streams.lst))
medium.reachable.streams<-unique(unlist(medium.reachable.streams.lst))
high.reachable.streams<-unique(unlist(high.reachable.streams.lst))

sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)
sdm@data<-sdm@data[,-c(216:227)]

sdm@data$low.reachable[match(low.reachable.streams,sdm$SegNo)]<-1
sdm@data$low.reachable[is.na(sdm$low.reachable)]<-0

sdm@data$med.reachable[match(medium.reachable.streams,sdm$SegNo)]<-1
sdm$med.reachable[is.na(sdm$med.reachable)]<-0

sdm@data$high.reachable[match(high.reachable.streams,sdm$SegNo)]<-1
sdm$high.reachable[is.na(sdm$high.reachable)]<-0

writeLinesShape(sdm,fn="Data/Shapfile/PCA3_reachable")




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



