setwd('P:/Projects/GitHub_Prj/BioVariability')

library(ggplot2)
library(vegan)
library(rgdal)
library(stringr)

CWsum<-read.csv("data/ColdWaterSites_ColdWaterFishSum.csv",header=TRUE)
ftaxa<-read.csv("data/FishSamplesColdWaterSites_012220.csv",header=TRUE,stringsAsFactors=FALSE)
pt<-read.csv("data/phylo_tree.csv",header=TRUE,stingsAsFactors=FALSE)

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
                        
##Percent of samples that decreased in abundance from early time period to later time period
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

###Format data for PhisViz##########
sample_data<-merge(ftaxa,pt,by.x="OTU",by.y="Taxon_name")
sample_data<-sample_data[sample_data$IsStocked=='FALSE',]
sites<-unique(sample_data[c("STA_SEQ","Station_Name","YLat","XLong")])
colnames(sites)[1]<-c("SID")
write.csv(sites,"sites.csv",row.names=FALSE,quote=FALSE)


sites<-read.csv("sites.csv",header=TRUE)
#Transform coordinates to numeric
sites$YLat  <- as.numeric(sites$YLat)
sites$XLong  <- as.numeric(sites$XLong)
sites.SP  <- SpatialPointsDataFrame(sites[,c(4,3)],
                                    sites[,-c(4,3)])
proj4string(sites.SP) <- CRS("+proj=utm +zone=18 +datum=WGS84") 
#proj4string(dataMap.SP) <- CRS("+init=epsg:4326") #WGS 84

str(sites.SP) # Now is class SpatialPointsDataFrame

#Write as geojson
writeOGR(sites.SP,"sitesfish",layer="sites", driver='GeoJSON',overwrite_layer = TRUE)

#Format phylo_tree
sample_data<-sample_data[,c(2,8,1,15)]
colnames(sample_data)<-c("SID","Collection_Date","Taxon_name","RelAbund")
sample_data$dateLen<-str_length(sample_data$Collection_Date)
sample_data$Year<-str_sub(sample_data$Collection_Date,sample_data$dateLen-3,sample_data$dateLen)
sample_data$Month<-0
sample_data$Day<-0

for (i in 1:dim(sample_data)[1]){
  sample_data$Day[i]<-str_sub(sample_data$Collection_Date[i],
                           str_locate(sample_data$Collection_Date[i],"/")[1]+1,
                           sample_data$dateLen[i]-5)
}

for (i in 1:dim(sample_data)[1]){
  sample_data$Month[i]<-str_sub(sample_data$Collection_Date[i],1,
                             str_locate(sample_data$Collection_Date[i],"/")[1]-1)
}

sample_data$Month<- ifelse(str_length(sample_data$Month)==2,sample_data$Month,paste0("0",sample_data$Month))
sample_data$Day<- ifelse(str_length(sample_data$Day)==2,sample_data$Day,paste0("0",sample_data$Day))
sample_data$Collection_Date<-paste0(sample_data$Year,"/",sample_data$Month,"/",sample_data$Day)
sample_data<-sample_data[,1:4]
sample_data<-sample_data[order(sample_data$SID,sample_data$Collection_Date),]

write.csv(sample_data,"sample_data_PV.csv",row.names=FALSE,quote=FALSE)











