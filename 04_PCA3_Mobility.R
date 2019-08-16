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
# when re-visit, skip the optimisation process and directly read in optimisation outcome.
#---
# read in species distribution data
library(maptools)
sdm<-readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)
species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,18,21,26)]  # delete 4 non-selecte species
n.sp<-colSums(species.distribution.df[,-1])

# read in 100 potential solutions for each conservation target (15%, 25% and 35%).
solution.top15<-readRDS("Data/R data/PCA3_solution_top15_non_Mob")
solution.top25<-readRDS("Data/R data/PCA3_solution_top25_non_Mob")
solution.top35<-readRDS("Data/R data/PCA3_solution_top35_non_Mob")

#solution.lst<-readRDS("Data/R data/PCA3_solution_top15_Mob")
#solution.lst<-readRDS("Data/R data/PCA3_solution_top25_Mob")
#solution.lst<-readRDS("Data/R data/PCA3_solution_top35_Mob")

# size of priority refuge network
n.seg.top15<-lengths(solution.top15)
n.seg.top25<-lengths(solution.top25)
n.seg.top35<-lengths(solution.top35)

summary(n.seg.top15)
summary(n.seg.top25)
summary(n.seg.top35)

# species representation
scaling.factor<-c(15,25,35)

cons.target.top15<-floor(n.sp*(scaling.factor[1]/100))  # change the index of the scaling.factor.
cons.target.top25<-floor(n.sp*(scaling.factor[2]/100))
cons.target.top35<-floor(n.sp*(scaling.factor[3]/100))

rep.sp.top15<-sapply(solution.top15,FUN = function(x){prot.species.df<-species.distribution.df[match(x,species.distribution.df$SegNo),]
                                                      rep.sp<-sum(colSums(prot.species.df[,-1])>=cons.target.top15)})
rep.sp.top25<-sapply(solution.top25,FUN = function(x){prot.species.df<-species.distribution.df[match(x,species.distribution.df$SegNo),]
                                                      rep.sp<-sum(colSums(prot.species.df[,-1])>=cons.target.top25)})
rep.sp.top35<-sapply(solution.top35,FUN = function(x){prot.species.df<-species.distribution.df[match(x,species.distribution.df$SegNo),]
                                                      rep.sp<-sum(colSums(prot.species.df[,-1])>=cons.target.top35)})

summary(rep.sp.top15)
summary(rep.sp.top25)
summary(rep.sp.top35)

#---
# selection frequency
#---
frequency.top15<-data.frame(table(unlist(solution.top15)))
frequency.top25<-data.frame(table(unlist(solution.top25)))
frequency.top35<-data.frame(table(unlist(solution.top35)))

summary(frequency.top15$Freq/100)
summary(frequency.top25$Freq/100)
summary(frequency.top35$Freq/100)

seq.network<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
names(seq.network)

seq.network$"selec_freq"<-0
seq.network$"selec_freq"[match(frequency.top15$Var1,seq.network$SegmentNo)]<-frequency.top15$Freq
#writeLinesShape(seq.network,fn="Data/Shapfile/PCA3 non Mob/Selec freq_top15")

seq.network<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
seq.network$"selec_freq"<-0
seq.network$"selec_freq"[match(frequency.top25$Var1,seq.network$SegmentNo)]<-frequency.top25$Freq
#writeLinesShape(seq.network,fn="Data/Shapfile/PCA3 non Mob/Selec freq_top25")

seq.network<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
seq.network$"selec_freq"<-0
seq.network$"selec_freq"[match(frequency.top35$Var1,seq.network$SegmentNo)]<-frequency.top35$Freq
#writeLinesShape(seq.network,fn="Data/Shapfile/PCA3 non Mob/Selec freq_top35")

#---
# best solution with the minimum objective function value
#---
para.a<-0.4   # position penalty weight
para.b<-0.5  # feature penalty weight

obj.func.top15<-sapply(solution.top15,FUN = function(x){
  prot.species.df<-species.distribution.df[match(x,species.distribution.df$SegNo),]
  RDI<-candidate.df$RDI[match(prot.species.df$SegNo,candidate.df$SegNo)]
  BC<-candidate.df$BC.nor[match(prot.species.df$SegNo,candidate.df$SegNo)]
  feature<-cons.target.top15-colSums(prot.species.df[,-1])
  feature<-sapply(feature,FUN = function(x) {ifelse(x<0,0,x)})
  objective.func<-sum(RDI)+sum(para.a*(1-BC))+sum(para.b*feature)
})
summary(obj.func.top15)

obj.func.top25<-sapply(solution.top25,FUN = function(x){
  prot.species.df<-species.distribution.df[match(x,species.distribution.df$SegNo),]
  RDI<-candidate.df$RDI[match(prot.species.df$SegNo,candidate.df$SegNo)]
  BC<-candidate.df$BC.nor[match(prot.species.df$SegNo,candidate.df$SegNo)]
  feature<-cons.target.top25-colSums(prot.species.df[,-1])
  feature<-sapply(feature,FUN = function(x) {ifelse(x<0,0,x)})
  objective.func<-sum(RDI)+sum(para.a*(1-BC))+sum(para.b*feature)
})
summary(obj.func.top25)

obj.func.top35<-sapply(solution.top35,FUN = function(x){
  prot.species.df<-species.distribution.df[match(x,species.distribution.df$SegNo),]
  RDI<-candidate.df$RDI[match(prot.species.df$SegNo,candidate.df$SegNo)]
  BC<-candidate.df$BC.nor[match(prot.species.df$SegNo,candidate.df$SegNo)]
  feature<-cons.target.top35-colSums(prot.species.df[,-1])
  feature<-sapply(feature,FUN = function(x) {ifelse(x<0,0,x)})
  objective.func<-sum(RDI)+sum(para.a*(1-BC))+sum(para.b*feature)
})
summary(obj.func.top35)

best.solution.top15<-solution.top15[[which(obj.func.top15==min(obj.func.top15))]]
best.solution.top25<-solution.top25[[which(obj.func.top25==min(obj.func.top25))]]
best.solution.top35<-solution.top35[[which(obj.func.top35==min(obj.func.top35))]]

seq.network<-readShapeLines("Data/Shapfile/Threshold of quant 0.5/PCA_Water only")
names(seq.network)

seq.network$"best_top15"<-0
seq.network$"best_top15"[match(best.solution.top15,seq.network$SegmentNo)]<-1

seq.network$"best_top25"<-0
seq.network$"best_top25"[match(best.solution.top25,seq.network$SegmentNo)]<-1

seq.network$"best_top35"<-0
seq.network$"best_top35"[match(best.solution.top35,seq.network$SegmentNo)]<-1

names(seq.network)
#writeLinesShape(seq.network,fn="Data/Shapfile/PCA3 non Mob/Best solution")

# identify reachable streams from best solution
SEQ.Clip<-readShapePoly(fn="Data/Shapfile/SEQ_Clip.shp")
hierarchy<-data.frame(site=SEQ.Clip$SegmentNo,nextds=SEQ.Clip$DWNID1)
sp.mobility<-read.csv("Data/R data/Species mobility.csv")

l.mob<-5
m.mob<-10
h.mob<-15

dam.segment<-c(859803,853302,856476,874709)  # consider the blocking effect of dams

low.reachable.streams.top15<-all.reachable.streams(best.solution.top15,mobility = l.mob)
medium.reachable.streams.top15<-all.reachable.streams(best.solution.top15,mobility = m.mob)
high.reachable.streams.top15<-all.reachable.streams(best.solution.top15,mobility = h.mob)

low.reachable.streams.top25<-all.reachable.streams(best.solution.top25,mobility = l.mob)
medium.reachable.streams.top25<-all.reachable.streams(best.solution.top25,mobility = m.mob)
high.reachable.streams.top25<-all.reachable.streams(best.solution.top25,mobility = h.mob)

low.reachable.streams.top35<-all.reachable.streams(best.solution.top35,mobility = l.mob)
medium.reachable.streams.top35<-all.reachable.streams(best.solution.top35,mobility = m.mob)
high.reachable.streams.top35<-all.reachable.streams(best.solution.top35,mobility = h.mob)


low.reachable.streams.top15<-unique(unlist(low.reachable.streams.top15))
medium.reachable.streams.top15<-unique(unlist(medium.reachable.streams.top15))
high.reachable.streams.top15<-unique(unlist(high.reachable.streams.top15))

low.reachable.streams.top25<-unique(unlist(low.reachable.streams.top25))
medium.reachable.streams.top25<-unique(unlist(medium.reachable.streams.top25))
high.reachable.streams.top25<-unique(unlist(high.reachable.streams.top25))

low.reachable.streams.top35<-unique(unlist(low.reachable.streams.top35))
medium.reachable.streams.top35<-unique(unlist(medium.reachable.streams.top35))
high.reachable.streams.top35<-unique(unlist(high.reachable.streams.top35))

lst_1<-c(low.reachable.streams.top15)
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}

lst_1<-c(medium.reachable.streams.top15)
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}

lst_1<-c(high.reachable.streams.top15)
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}

lst_1<-c(low.reachable.streams.top25)
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}

lst_1<-c(medium.reachable.streams.top25)
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}

lst_1<-c(high.reachable.streams.top25)
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}

lst_1<-c(low.reachable.streams.top35)
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}

lst_1<-c(medium.reachable.streams.top35)
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}

lst_1<-c(high.reachable.streams.top35)
for(i in 1:length(lst_1)){
  cat("\"SegmentNo\"","=",lst_1[i],"OR",file = paste0("Ups-",lst_1[1],".txt"),append = TRUE)
}

#---
# species representation
#---

# sp representation of the best solution: 15.4%,26.2% and 35.9%.
best.species.df<-species.distribution.df[match(best.solution,species.distribution.df$SegNo),]

plot(colSums(best.species.df[,c(2:26)])/n.sp)
mean(colSums(best.species.df[,c(2:26)])/colSums(species.distribution.df[,c(2:26)]))

# sp rep for randomly selected segments in the same number as that systematically selected
#random.seg<-PCA1_SegNo[runif(357,min = 1,max=length(PCA1_SegNo))]
#plot(colSums(species.distribution.df[match(random.seg,species.distribution.df$SegNo),-1])/n.sp)


#****************Optimisation process**************************

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

dam.segment<-c(859803,853302,856476,874709)  # consider the blocking effect of dams

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



