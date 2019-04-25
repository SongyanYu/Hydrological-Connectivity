#----------------
# This script aims to calculate correlation between the spatial distribution of each species and that predicted by ROSE
# Author: Songyan Yu
# Date create: 25/04/2019
#----------------

setwd("D:/New folder/Google Drive/PhD at GU/Part 4 Hydrologic connectivity")
dat<-read.csv("Data/R data/Species distribution.csv")

cor.test(dat$Fig._4,dat$Rose)
