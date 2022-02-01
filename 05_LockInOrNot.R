#---
# script compares priority network between locking in protected areas and not locking in.
# author: Songyan Yu
# date created: 31/01/2022
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

# find which solution is best
best.solution <- 
  do.call(rbind.data.frame, priority.lst) %>%
  filter(repsp == 25) %>%
  group_by(group, target) %>%
  filter(objfun == min(objfun))

# find best solution for each target
solution.files <- list.files("../../Data/R_data/",
                             pattern = "04_PCA_solution",
                             full.names = TRUE)

best.solution.lst <- lapply(1:length(solution.files), FUN = function(x){
  priority.seg <- readRDS(solution.files[x])[[best.solution$no[x]]]
})

# read in protected areas
protectedArea.shp <- maptools::readShapeLines("../../Data/Shapfile/SEQ_networks_strahler02_withinProtectedAreas.shp")
protected.seg <- protectedArea.shp$SegmentNo
length(unique(protected.seg)) == length(protected.seg) # check if any duplicates

# find newly needed segment for each target
best.solution.without.protected.areas.lst <- lapply(best.solution.lst, FUN = function(x){
  setdiff(x, protected.seg)
})

best.solution$addedSegSize = lengths(best.solution.without.protected.areas.lst)

# calculate length of priority network
SEQ.network <- maptools::readShapeLines(fn = "../../Data/Shapfile/SEQ_networks_strahler02.shp")
SEQ.network$RCHLEN

best.solution.length <- sapply(best.solution.lst, FUN = function(x){
  sum(SEQ.network$RCHLEN[match(x, SEQ.network$SegmentNo)], na.rm = TRUE)
})

best.solution.without.protected.areas.length <- sapply(best.solution.without.protected.areas.lst, FUN = function(x){
  sum(SEQ.network$RCHLEN[match(x, SEQ.network$SegmentNo)], na.rm = TRUE)
})

# visualise priority network length
library(ggplot2)
library(reshape2)
library(tidyr)

best.solution %>%
  ungroup %>%
  mutate(length = best.solution.length,
         addedSegLength = best.solution.without.protected.areas.length,
         protSegLength = best.solution.length - best.solution.without.protected.areas.length,
         group = factor(group, levels = c("NO", "protectedAreas"), labels = c("No locked in", "Locked in")),
         target = factor(target, levels = c("015", "025", "035"), labels = c("15%", "25%", "35%"))) %>%
  select(target, addedSegLength, protSegLength, group) %>%
  pivot_longer(cols = c(addedSegLength, protSegLength)) %>%
  mutate(name = factor(name, levels = c("addedSegLength", "protSegLength"),
                       labels = c("Outside protected areas", "Inside protected areas")),
         value = round(value, 0)) %>%
  ggplot(aes(x = group, y = value, fill = name)) +
  geom_bar(position = 'stack', stat = 'identity') +
  facet_grid(~target) +
  theme_bw() +
  ylab("Length of priority network (km)") +
  xlab("") +
  theme(legend.position = 'top', legend.title = element_blank()) +
  geom_text(aes(label = value),size = 3, position = position_stack(vjust = 0.5))
ggsave(filename = "../../Figures/Fig05_PriorityNetworkSize.jpg",
       dpi = 800, width = 6, height = 4)




