setwd("/home/mkozlak/Documents/Projects/GitHub/BioVariability")

library(ggplot2)
library(vegan)
library(rgdal)
library(stringr)

tier<-5
attr<-4

sites<-read.csv("data/Stations.csv",header=TRUE)
sites<-sites[,c(1,5,6,9,10)]
BCG<-read.csv(paste0("data/bugBCG",tier,"_MutliYr_02112020.csv"),header=TRUE)
#BCG$YrGrp<-ifelse(BCG$SampYr< 1999,1,ifelse(BCG$SampYr> 2008,3,2))
BCG$YrGrp<-ifelse(BCG$SampYr< 1996,1,
                  ifelse(BCG$SampYr>1995&BCG$SampYr<2001,2,
                         ifelse(BCG$SampYr>2000&BCG$SampYr<2006,3,
                                ifelse(BCG$SampYr>2005&BCG$SampYr<2011,4,
                                       ifelse(BCG$SampYr>2010&BCG$SampYr<2016,5,6)))))
taxa<-read.csv(paste0("data/taxaData_Tier",tier,".csv"),header=TRUE)
# mtaxa<-read.csv("data/masterTaxalist.csv",header=TRUE)
# colnames(mtaxa)[2]<-"TaxonNameCurrent"
# taxa<-merge(taxa,mtaxa[,c(2,5)],by="TaxonNameCurrent")
# taxa$EPT<-ifelse(taxa$ORDER=="Ephemeroptera"|taxa$ORDER=="Plecoptera"|taxa$ORDER=="Trichoptera",1,0)



#cnttaxaFam<-unique(taxaFam[c("STA_SEQ","Station_Name","SampYr")])
#aggregate(SampYr~STA_SEQ+Station_Name,data=cnttaxaFam,FUN=length)

#taxa_agg<-taxa[which(taxa$BCG_Attribute==2|taxa$BCG_Attribute==3),]
taxa_agg<-taxa[which(taxa$BCG_Attribute==4),]
taxaFam<-unique(taxa_agg[c("STA_SEQ","Station_Name","SampYr","FAMILY")])
taxaFam<-aggregate(FAMILY~STA_SEQ+Station_Name+SampYr,data=taxaFam,FUN=length)
#taxa_agg<-taxa[which(taxa$EPT==1),]
taxa_agg<-aggregate(taxa_agg$ABUNDANCE,by=as.list(taxa_agg[,c("STA_SEQ","SampYr")]),FUN=sum)
# taxa_agg<-aggregate(taxaFam$EPT,by=as.list(taxaFam[,c("STA_SEQ","SampYr")]),FUN=sum)
BCG<-merge(BCG,taxa_agg,by=c("STA_SEQ","SampYr"),all.x=TRUE)
BCG$x[is.na(BCG$x)]<-0
BCG$RelAbund<-BCG$x/BCG$SumOfABUNDANCE
BCG<-merge(BCG,taxaFam,by=c("STA_SEQ","Station_Name","SampYr"),all.x=TRUE)


strNames<-unique(BCG[c("STA_SEQ","Station_Name","CountOfSampYr")])


##Summary of Differences by Time Period
grpSum<-aggregate(BCG$RelAbund~STA_SEQ+YrGrp,data=BCG,FUN=mean)
grpWide<-reshape(grpSum,idvar="STA_SEQ",timevar="YrGrp",direction="wide")
#colnames(grpWide)<-c("STA_SEQ","Yr1","Yr2","Yr3")#rename columns
colnames(grpWide)<-c("STA_SEQ","Yr1","Yr2","Yr3","Yr4","Yr5","Yr6")#rename columns


grpWide<-merge(grpWide,strNames,by="STA_SEQ")
grpWide$cntNA<- apply(grpWide, 1, function(x) sum(is.na(x)))

ModTaxaOverTime<-grpWide[grpWide$cntNA<2,]
ModTaxaOverTime<-merge(ModTaxaOverTime,sites,by="STA_SEQ")

write.csv(ModTaxaOverTime,paste0("BCG",tier,"TaxaTier",attr,"_OverTime.csv"),row.names=FALSE)

#grpWide[complete.cases(grpWide),]
#grpWide[grpWide$CountOfSampYr>=10,]


