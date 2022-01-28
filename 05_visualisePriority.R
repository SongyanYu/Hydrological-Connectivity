#---
# visualise prioritisation results
# author: Songyan Yu
# date created: 28/01/2022
#---

library(dplyr)

# species representation
rep.sp.files <- list.files("../../Data/R_data/",
                           pattern = "04_PCA_repSp",
                           full.names = TRUE)

# objective function values
objfun.files <- list.files(path = "../../Data/R_data/",
                           pattern = "04_PCA_objFunc",
                           full.names = TRUE)

# solution size
size.files <- list.files("../../Data/R_data/",
                         pattern = "04_PCA_size",
                         full.names = TRUE)

target <- lapply(size.files, FUN = function(x){
  strsplit(strsplit(strsplit(x, split = "/")[[1]][5], split = ".RDS")[[1]][1], split = "_")[[1]][4]
})
group <- lapply(size.files, FUN = function(x){
  strsplit(strsplit(strsplit(x, split = "/")[[1]][5], split = ".RDS")[[1]][1], split = "_")[[1]][5]
})

# combine
priority.lst <- lapply(1:6, FUN = function(x){
  size <- readRDS(size.files[x])
  objfun <- readRDS(objfun.files[x])
  repsp <- readRDS(rep.sp.files[x])
  
  data.frame(size, objfun, repsp) %>%
    mutate(target = target[[x]],
           group = group[[x]],
           no = 1:100)
})

best.solution <- 
  do.call(rbind.data.frame, priority.lst) %>%
  filter(repsp == 25) %>%
  group_by(group, target) %>%
  filter(objfun == min(objfun))

SEQ.network <- maptools::readShapeLines(fn = "../../Data/Shapfile/SEQ_networks_strahler02.shp")
terra::plot(SEQ.network)

solution.files <- list.files("../../Data/R_data/",
                             pattern = "04_PCA_solution",
                             full.names = TRUE)

best.solution.lst <- lapply(1:length(solution.files), FUN = function(x){
  priority.seg <- readRDS(solution.files[x])[[best.solution$no[x]]]
  
  priority.df <- data.frame(SegNo = SEQ.network$SegmentNo,
                            priority = 0)
  
  priority.df$priority[priority.df$SegNo %in% priority.seg] <- 1
  colnames(priority.df)[2] <- paste0(group[[x]],"_",target[[x]])
  priority.df
})

best.solution.df <- do.call(cbind.data.frame, best.solution.lst)[,c(1,2,4,6,8,10,12)]
colnames(best.solution.df)

SEQ.network@data %>%
  left_join(., best.solution.df, by = c("SegmentNo" = "SegNo"))
maptools::writeLinesShape(SEQ.network, fn = "../../Data/Shapfile/05_visualisePriority/SEQ_priority")





