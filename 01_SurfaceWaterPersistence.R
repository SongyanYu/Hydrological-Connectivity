#------------------------
# This script aims to quantify annual surface water persistence
# Author: Songyan Yu
# Creat date: 18/12/2018
# update: 25/01/2022
#------------------------

daily.sw.1911.2017<-readRDS("Data/R data/wet50-days_1911-2017")

daily.sw.1911.2017.df<-data.frame(daily.sw.1911.2017)
colnames(daily.sw.1911.2017.df)<-c(1911:2017)
sum(is.na(daily.sw.1911.2017.df))  # should be "0"