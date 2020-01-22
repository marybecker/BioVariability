setwd('/home/mkozlak/Documents/Projects/GitHub/BioVariability')

library(ggplot2)
library(vegan)
library(rgdal)
library(stringr)

CWsum<-read.csv("data/ColdWaterSites_ColdWaterFishSum.csv",header=TRUE)
ftaxa<-read.csv("data/FishSamplesColdWaterSites_012220.csv",header=TRUE)

grpSum<-aggregate(CWsum$MaxOfFishPer100M,by=as.list(CWsum[,c("STA_SEQ","YearGrp")]),FUN=max)
grpWide<-reshape(grpSum,idvar="STA_SEQ",timevar="YearGrp",direction="wide")
colnames(grpWide)<-c("SID","Yr1","Yr2","Yr3")#rename columns
grpWide$Y1Y3<-grpWide$Yr3-grpWide$Yr1
grpWide$Yr1Y2<-grpWide$Yr2-grpWide$Yr1
grpWide$Yr2Yr3<-grpWide$Yr3-grpWide$Yr2
grpWide$ColdYr1<-ifelse(grpWide$Yr1>=10,1,ifelse(is.na(grpWide$Yr1),NA,0))
grpWide$ColdYr2<-ifelse(grpWide$Yr2>=10,1,ifelse(is.na(grpWide$Yr2),NA,0))
grpWide$ColdYr3<-ifelse(grpWide$Yr3>=10,1,ifelse(is.na(grpWide$Yr3),NA,0))
grpWide$ColdY1Y3<-grpWide$ColdYr3-grpWide$ColdYr1
grpWide$ColdYr1Y2<-grpWide$ColdYr2-grpWide$ColdYr1
grpWide$ColdYr2Yr3<-grpWide$ColdYr3-grpWide$ColdYr2
                        
##Percent of samples that decreased from early time period to later time period
t1t3Decs<-dim(grpWide[which(grpWide[,5]<0),])[1]/dim(grpWide[complete.cases(grpWide[,5]),])[1]
t1t2Decs<-dim(grpWide[which(grpWide[,6]<0),])[1]/dim(grpWide[complete.cases(grpWide[,6]),])[1]
t2t3Decs<-dim(grpWide[which(grpWide[,7]<0),])[1]/dim(grpWide[complete.cases(grpWide[,7]),])[1]

##Changes in Cold Water Category Across Different Time Periods
t1t3DescCold<-dim(grpWide[which(grpWide[,11]==-1),])[1]/dim(grpWide[complete.cases(grpWide[,11]),])[1]
t1t3IncCold<-dim(grpWide[which(grpWide[,11]==1),])[1]/dim(grpWide[complete.cases(grpWide[,11]),])[1]
t1t3Stable<-dim(grpWide[which(grpWide[,11]==0),])[1]/dim(grpWide[complete.cases(grpWide[,11]),])[1]

t1t2DescCold<-dim(grpWide[which(grpWide[,12]==-1),])[1]/dim(grpWide[complete.cases(grpWide[,12]),])[1]
t1t2IncCold<-dim(grpWide[which(grpWide[,12]==1),])[1]/dim(grpWide[complete.cases(grpWide[,12]),])[1]
t1t2Stable<-dim(grpWide[which(grpWide[,12]==0),])[1]/dim(grpWide[complete.cases(grpWide[,12]),])[1]

t2t3DescCold<-dim(grpWide[which(grpWide[,13]==-1),])[1]/dim(grpWide[complete.cases(grpWide[,13]),])[1]
t2t3IncCold<-dim(grpWide[which(grpWide[,13]==1),])[1]/dim(grpWide[complete.cases(grpWide[,13]),])[1]
t2t3Stable<-dim(grpWide[which(grpWide[,13]==0),])[1]/dim(grpWide[complete.cases(grpWide[,13]),])[1]









