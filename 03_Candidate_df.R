#--------------------
# Produce "candidate.df" for systematic prioritisation to use.
# Date: 04/07/2019
# update: 25/01/2022
#--------------------

# to get 'candidate.df'
source("02_CandidateRefuges.R")

#------------
# RDI
#------------
SEQ.clip<-readShapePoly("../../Data/Shapfile/SEQ_Clip.shp")
RDI.df<-data.frame(SegNo=SEQ.clip$SegmentNo,RDI=SEQ.clip$RDI)

candidate.df <- 
  candidate.df %>%
  left_join(., RDI.df, by = "SegNo")

sum(is.na(candidate.df$RDI))  # should be "0"

#--------------
# Betweenness centrality (BC)
#--------------
BC.mar<-readRDS("../../Data/R_data/BC.maroochy")
BC.sth<-readRDS("../../Data/R_data/BC.southcoast")
BC.pin<-readRDS("../../Data/R_data/BC.pine")
BC.log<-readRDS("../../Data/R_data/BC.logan")
BC.bne<-readRDS("../../Data/R_data/Betweenness Centrality_BNE")

# To make BC for stream segments from different river networks comparable,
# BC values were normalised within each river network before they were compared.
BC.mar$BC.nor<-(BC.mar$BC-min(BC.mar$BC))/(max(BC.mar$BC)-min(BC.mar$BC))
BC.sth$BC.nor<-(BC.sth$BC-min(BC.sth$BC))/(max(BC.sth$BC)-min(BC.sth$BC))
BC.pin$BC.nor<-(BC.pin$BC-min(BC.pin$BC))/(max(BC.pin$BC)-min(BC.pin$BC))
BC.log$BC.nor<-(BC.log$BC-min(BC.log$BC))/(max(BC.log$BC)-min(BC.log$BC))
BC.bne$BC.nor<-(BC.bne$BC-min(BC.bne$BC))/(max(BC.bne$BC)-min(BC.bne$BC))

BC.SEQ<-rbind(BC.mar,BC.sth,BC.pin,BC.log,BC.bne)

candidate.df <- left_join(candidate.df, BC.SEQ[,-2], by="SegNo")
sum(is.na(candidate.df$BC.nor)) # should be '0'

#----------------
# Plotting BC
#----------------

plot = FALSE

if (plot) {
  BC.SEQ$class[BC.SEQ$BC.nor>=0.71]<-1
  BC.SEQ$class[BC.SEQ$BC.nor>=0.38&BC.SEQ$BC.nor<0.71]<-2
  BC.SEQ$class[BC.SEQ$BC.nor>=0.12&BC.SEQ$BC.nor<0.38]<-3
  BC.SEQ$class[BC.SEQ$BC.nor<0.12]<-4
  
  top5<-quantile(BC.SEQ$BC.nor,probs = 0.95)
  BC.SEQ$top5[BC.SEQ$BC.nor>=top5]<-1
  BC.SEQ$top5[BC.SEQ$BC.nor<top5]<-0
  
  library(dplyr)
  temp@data<-left_join(temp@data,BC.SEQ,by=c("SegmentNo"="SegNo"))
  writeLinesShape(temp,fn = "Data/Shapfile/Betweenness centrality/SEQ_BC")
}

#---
# Hydrological penalty - # year meeting the two criteria
#---
candidate.df <- 
  candidate.df %>%
  mutate(Freq.nor = (Freq - min(Freq))/(max(Freq) - min(Freq))) %>%
  select(-Freq)

saveRDS(candidate.df, file = "../../Data/R_data/03_Candidate_df.RDS")
