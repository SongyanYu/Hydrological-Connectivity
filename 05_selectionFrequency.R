#---
# script to calcualte selection frequency of stream segments
# author: Songyan Yu
# date created: 28/01/2022
#---

# find solution files
solution.files <- list.files("../../Data/R_data/",
                             pattern = "04_PCA_solution",
                             full.names = TRUE)

# extract conservation target info
target <- lapply(solution.files, FUN = function(x){
  strsplit(strsplit(strsplit(x, split = "/")[[1]][5], split = ".RDS")[[1]][1], split = "_")[[1]][4]
})
group <- lapply(solution.files, FUN = function(x){
  strsplit(strsplit(strsplit(x, split = "/")[[1]][5], split = ".RDS")[[1]][1], split = "_")[[1]][5]
})

# calculate selection frequency
freq.lst <- lapply(1:length(solution.files), FUN = function(x){
  solution.seg <- readRDS(solution.files[x])
  freq.df <- as.data.frame(table(unlist(solution.seg)))
  colnames(freq.df) <- c("SegNo", paste0(group[[x]],"_",target[[x]]))
  freq.df$SegNo <- as.character(freq.df$SegNo)
  freq.df
})

#freq.df <- do.call(rbind.data.frame, freq.lst)

# combine selection freq with river networks
SEQ.network <- maptools::readShapeLines(fn = "../../Data/Shapfile/SEQ_networks_strahler02.shp")
terra::plot(SEQ.network)

library(dplyr)

for(i in 1:length(freq.lst)){
  selecfreq.df <- data.frame(SegNo = as.character(SEQ.network$SegmentNo)) %>%
    left_join(., freq.lst[[i]], by = "SegNo")
  selecfreq.df[is.na(selecfreq.df[,2]),2] <- 0
  
  if(identical(selecfreq.df$SegNo, as.character(SEQ.network$SegmentNo))){
    SEQ.network@data <- cbind(SEQ.network@data, selecfreq.df)
  }else{
    cat("Error in combining selection frequency with river network.\n")
  }
  
  maptools::writeLinesShape(SEQ.network, 
                            fn = paste0("../../Data/Shapfile/05_selectionFrequency/SEQ_selecFreq_",
                                        paste0(group[[i]],"_",target[[i]])))
}





