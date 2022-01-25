#---
# this script systematically prioritise refuges for conservation
# the prioritisation incorporates current protected areas.
# author: Songyan Yu
# date created: 24/01/2022
#---

setwd("../../")

# read in candidate data frame
candidate.df <- readRDS("Data/R_data/03_Candidate_df.RDS")
candidate.seg <- candidate.df$SegNo

# read in species distribution
library(maptools)
sdm <- readShapeLines(fn="Data/Shapfile/Species distribution model/PCA_Naive_Species.shp")
names(sdm)
species.distribution.df<-sdm@data[,c(186:215)] #check the col number
species.distribution.df<-species.distribution.df[,-c(14,18,21,26)]  # delete 4 non-selected species
n.sp<-colSums(species.distribution.df[,-1])

# read in stream segments within protected areas
protectedArea.shp <- maptools::readShapeLines("Data/Shapfile/SEQ_networks_strahler02_withinProtectedAreas.shp")
protected.seg <- protectedArea.shp$SegmentNo
unique(protected.seg) # check if any duplicates
sum(protected.seg %in% candidate.seg) / length(protected.seg)

# candidate stream segments that are in the protected areas
fixed.priority.seg <- candidate.seg[candidate.seg %in% protected.seg]
flexible.seg <- setdiff(candidate.seg, fixed.priority.seg)

# species representation by those stream segments in protected areas
prot.species.fixed <- species.distribution.df[match(fixed.priority.seg,species.distribution.df$SegNo),]

# conservation targets
scaling.factor <- c(15,25,35)

# systematic prioritisation
for(m in 1:length(scaling.factor)){
  
  cons.target <- floor(n.sp*(scaling.factor[m]/100))
  
  para.a <- 1   # hydrology penalty weight
  para.b <- 1   # position penalty weight
  para.c <- 1   # feature penalty weight
  rep.sp <- c()
  n.seg <- c()
  solution.lst <- list()
  obj.func <- c()
  int.obj.fun <- c()
  int.obj.fun.lst <- list()

  for(i in 1:100){
    
    remaining.seg <- flexible.seg
    #------
    # 1. start an initial configuration by randomly selecting a stream segment.
    #------
    n <- sample(1:length(remaining.seg), 1)
    
    #------
    # 2. calculate the value of objective function
    #------
    prot.species.df <- species.distribution.df[match(remaining.seg[n],species.distribution.df$SegNo),]
    prot.species.df <- rbind(prot.species.fixed, prot.species.df)
    
    RDI <- candidate.df$RDI[match(prot.species.df$SegNo, candidate.df$SegNo)]
    freq <- candidate.df$Freq.nor[match(prot.species.df$SegNo, candidate.df$SegNo)]
    BC <- candidate.df$BC.nor[match(prot.species.df$SegNo, candidate.df$SegNo)]
    feature <- cons.target - colSums(prot.species.df[,-1])
    feature <- sapply(feature, FUN = function(x) {ifelse(x<0,0,x)})
    
    objective.func <- 
      sum(RDI) + 
      sum(para.a * (1 - freq)) + 
      sum(para.b * (1 - BC)) + 
      sum(para.c * feature)
    
    int.obj.fun <- objective.func
    #------
    # 3. iteratively try a new segment to add to the current configuration
    # 1) randomly select a new segment
    # 2) check whether the new segment can improve/worsen the value of objective function
    # 3) if improving (reducing), keep it in; if worsening (increasing), remove it.
    #------
    iteration.n <- 0
    while(iteration.n <= 1000){
      
      # 1) randomly select a new segment
      n <- sample(1:length(remaining.seg), 1)
      
      # check whether the selected segment is in the current configuration.
      if(remaining.seg[n] %in% prot.species.df$SegNo){
        prot.species.df.2 <- prot.species.df[-match(remaining.seg[n], prot.species.df$SegNo),]
        
        RDI <- candidate.df$RDI[match(prot.species.df.2$SegNo, candidate.df$SegNo)]
        freq <- candidate.df$Freq.nor[match(prot.species.df.2$SegNo, candidate.df$SegNo)]
        BC <- candidate.df$BC.nor[match(prot.species.df.2$SegNo, candidate.df$SegNo)]
        feature <- cons.target - colSums(prot.species.df.2[,-1])
        feature <- sapply(feature, FUN = function(x) {ifelse(x<0,0,x)})
        
        objective.func.2 <- 
          sum(RDI) +
          sum(para.a * (1 - freq)) +
          sum(para.b * (1 - BC)) +
          sum(para.c * feature)  # need to adjust para.b
        
        if(objective.func >= objective.func.2){
          prot.species.df <- prot.species.df.2   # update current configuration
          objective.func <- objective.func.2    # update the value of objective function
          
          int.obj.fun <- append(int.obj.fun, objective.func)
          
          iteration.n <- 0
        }else{
          #cat(iteration.n,"\n")
          iteration.n <- iteration.n + 1
        }
      }else{
        prot.species <- species.distribution.df[match(remaining.seg[n], species.distribution.df$SegNo), ]
        prot.species.df.2 <- rbind(prot.species.df, prot.species)
        
        # 2) check whether the new segment can improve/worsen the value of objective function
        # calculate the value of objective function
        
        RDI <- candidate.df$RDI[match(prot.species.df.2$SegNo, candidate.df$SegNo)]
        freq <- candidate.df$Freq.nor[match(prot.species.df.2$SegNo, candidate.df$SegNo)]
        BC <- candidate.df$BC.nor[match(prot.species.df.2$SegNo, candidate.df$SegNo)]
        feature <- cons.target - colSums(prot.species.df.2[, -1])
        feature <- sapply(feature, FUN = function(x) {ifelse(x<0,0,x)})
        
        objective.func.2 <- 
          sum(RDI) +
          sum(para.a * (1 - freq)) +
          sum(para.b * (1 - BC)) +
          sum(para.c * feature)  # need to adjust para.b
        
        if(objective.func >= objective.func.2){
          prot.species.df <- prot.species.df.2   # update current configuration
          objective.func <- objective.func.2    # update the value of objective function
          
          int.obj.fun <- append(int.obj.fun, objective.func)
          
          iteration.n <- 0
        }
        else{
          iteration.n <- iteration.n + 1
        }
      }
    }
    
    n.seg[i] <- nrow(prot.species.df)
    solution.lst[[i]] <- prot.species.df$SegNo
    rep.sp[i] <- sum(colSums(prot.species.df[, -1]) >= cons.target)
    obj.func[i] <- objective.func
    int.obj.fun.lst[[i]] <- int.obj.fun
    
    cat("new loop ",i,"\n")
  }
  saveRDS(solution.lst, file = paste0("Data/R_data/04_PCA_protectedAreas_solution_0", scaling.factor[m],".RDS"))
  saveRDS(n.seg, file = paste0("Data/R_data/P04_PCA_protectedAreas_size_0", scaling.factor[m], ".RDS"))
  saveRDS(rep.sp, file = paste0("Data/R_data/04_PCA_protectedAreas_repSp_0", scaling.factor[m], ".RDS"))
  saveRDS(obj.func, file = paste0("Data/R_data/04_PCA_protectedAreas_objFunc_0", scaling.factor[m], ".RDS"))
  saveRDS(int.obj.fun.lst, file = paste0("Data/R_data/04_PCA_protectedAreas_intObjFunLst_0", scaling.factor[m],".RDS"))
}

