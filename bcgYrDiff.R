setwd('P:/Projects/GitHub_Prj/BioVariability')

library(ggplot2)
library(vegan)

bcg<-read.csv("data/bug_bcg_5Yrs_Tier2.csv",header=TRUE)
taxa<-read.csv("data/taxaData_5Yr_Tier2.csv",header=TRUE,stringsAsFactors = FALSE)
taxa$BCG_Attribute<-ifelse(is.na(taxa$BCG_Attribute),0,taxa$BCG_Attribute)





##########Raw Taxa ###########################################################
##############################################################################
taxa$OTU<-taxa$GENUS
for (i in 1:dim(taxa)[1]){
  if(taxa$GENUS[i]=='Na'&taxa$FAMILY[i]=='Na'){
    taxa$OTU[i]<-taxa$TaxonNameCurrent[i]
  }
  if(taxa$FAMILY[i]!='Na'&taxa$GENUS[i]=='Na'){
    taxa$OTU[i]<-taxa$FAMILY[i]
  }
  
}
taxa<-taxa[,c(1,2,11,12,16,17,18,10)]
taxa<-aggregate(taxa$ABUNDANCE,by=as.list(taxa[,c(1:7)]),FUN=sum)
taxaSite<-unique(taxa$STA_SEQ)
sentaxa<-taxa[taxa$BCG_Attribute==2.0,]
toltaxa<-taxa[taxa$BCG_Attribute>=3.0,]

###Cnt of total individuals
samplecnt<-taxa[,c(1,2,4,8)]
samplecnt<-aggregate(samplecnt$x,by=as.list(samplecnt[,c(1:3)]),FUN=sum)
samplecntSen<-sentaxa[,c(1,2,4,8)]
samplecntSen<-aggregate(samplecntSen$x,by=as.list(samplecntSen[,c(1:3)]),FUN=sum)

samplecnt<-merge(samplecnt,samplecntSen,by=c("STA_SEQ","Station_Name","SampYr"))
colnames(samplecnt)[4:5]<-c("Total","T2T3")
samplecnt$Pct<-samplecnt$T2T3/samplecnt$Total

aggregate(samplecnt$Pct,by=as.list(samplecnt[,1:2]),FUN=mean)



samplecnt$YrGrp<-ifelse(samplecnt$SampYr<=1999,"1989-1999",
                  ifelse (samplecnt$SampYr <=2009,"2000-2009","2010-2017"))
samplecntAvg<- aggregate(samplecnt$Pct,by=as.list(samplecnt[,c(1,2,7)]),FUN=mean)

for (i in 1:length(taxaSite)){
  site<-(samplecntAvg[samplecntAvg$STA_SEQ==taxaSite[i],])
  plot<- ggplot(site,aes(YrGrp,x))+
    geom_col()+
    labs(title=paste(site$Station_Name,site$STA_SEQ),x=NULL,y="Avg Pct Individuals Tier 2 3")+
    ylim(0,max(samplecntAvg$x))
  ggsave(paste0(site$STA_SEQ,"_",site$Station_Name,"_IndCnt.jpg"),device="jpeg",width = 4, height = 4,)
}

samplecnt<-taxa[,c(1,2,4,8)]
samplecnt<-aggregate(samplecnt$x,by=as.list(samplecnt[,c(1:3)]),FUN=sum)
samplecntTol<-toltaxa[,c(1,2,4,8)]
samplecntTol<-aggregate(samplecntTol$x,by=as.list(samplecntTol[,c(1:3)]),FUN=sum)
samplecntTol<-merge(samplecnt,samplecntTol,by=c("STA_SEQ","Station_Name","SampYr"))
colnames(samplecntTol)[4:5]<-c("Total","T5")
samplecntTol$Pct<-samplecntTol$T5/samplecntTol$Total
samplecntTol$YrGrp<-ifelse(samplecntTol$SampYr<=1999,"1989-1999",
                        ifelse (samplecntTol$SampYr <=2009,"2000-2009","2010-2017"))
samplecntTolAvg<- aggregate(samplecntTol$Pct,by=as.list(samplecntTol[,c(1,2,7)]),FUN=mean)

for (i in 1:length(taxaSite)){
  site<-(samplecntTolAvg[samplecntTolAvg$STA_SEQ==taxaSite[i],])
  plot<- ggplot(site,aes(YrGrp,x))+
    geom_col()+
    labs(title=paste(site$Station_Name,site$STA_SEQ),x=NULL,y="Avg Pct Individuals Tier 5")+
    ylim(0,max(samplecntTolAvg$x))
  ggsave(paste0(site$STA_SEQ,"_",site$Station_Name,"_IndTolPct.jpg"),device="jpeg",width = 4, height = 4,)
}

###Cnt of OTU
cntOTUSen<-unique(sentaxa[c("STA_SEQ","Station_Name","SampYr","OTU")])
cntOTUSen<-aggregate(cntOTUSen$OTU,by=as.list(cntOTUSen[,c(1:3)]),FUN=length)
cntOTU<-unique(taxa[c("STA_SEQ","Station_Name","SampYr","OTU")])
cntOTU<-aggregate(cntOTU$OTU,by=as.list(cntOTU[,c(1:3)]),FUN=length)
cntOTU<-merge(cntOTU,cntOTUSen,by=c("STA_SEQ","Station_Name","SampYr"))
colnames(cntOTU)[4:5]<-c("Total","T2T3")
cntOTU$Pct<-cntOTU$T2T3/cntOTU$Total

cntOTU$YrGrp<-ifelse(cntOTU$SampYr<=1999,"1989-1999",
                        ifelse (cntOTU$SampYr <=2009,"2000-2009","2010-2017"))
cntOTUAvg<- aggregate(cntOTU$Pct,by=as.list(cntOTU[,c(1,2,7)]),FUN=mean)

for (i in 1:length(taxaSite)){
  site<-(cntOTUAvg[cntOTUAvg$STA_SEQ==taxaSite[i],])
  plot<- ggplot(site,aes(YrGrp,x))+
    geom_col()+
    labs(title=paste(site$Station_Name,site$STA_SEQ),x=NULL,y="Pct OTU 2-3")+
    ylim(0,max(cntOTUAvg$x))
  ggsave(paste0(site$STA_SEQ,"_",site$Station_Name,"_CntOTU.jpg"),device="jpeg",width = 4, height = 4,)
}


##Reshape taxa to matrix format for taxa analysis
taxaWide<-taxa[,c(1,2,3,4,7,8)]
taxaWide$SampID<-paste0(taxaWide$STA_SEQ,"YR",taxaWide$SampYr)
taxaWide<-taxaWide[,c(7,5,6)]
taxaWide<-reshape(taxaWide,idvar="SampID",timevar="OTU",direction="wide")
taxaWide[is.na(taxaWide)]<-0


head(taxa[order(taxa$STA_SEQ,taxa$SampYr),],10)
samples<-unique(taxa[c("STA_SEQ","Station_Name","SampYr")])


x<=0
bcgYrDiff<-data.frame(STA_SEQ=integer(),SampYr=integer(),TaxaRich=)#Create empty dataframe to store combinations
for i in dim(samples)[1]{
  dim(taxa[taxa$STA_SEQ==samples$STA_SEQ[i]&taxa$SampYr==samples$SampYr[i],])[1]
}


##########BCG and Flow #######################################################
##############################################################################                 
flow<-read.csv("data/findex_1989_2017.csv",header=TRUE)
flow<-flow[,c(1,14)]
colnames(flow)[1]<-"tripdate"
bcg<-merge(bcg,flow,by.x="tripdate",all.x=TRUE)
sites<-unique(bcg$STA_SEQ)
sampleflow<-unique(bcg[c("STA_SEQ","SampYr","index")])
bcgAvgIndex<-aggregate(bcg$index,by=as.list(bcg[,c(2,3,5)]),FUN=mean)
ggplot(bcgAvgIndex[bcgAvgIndex$STA_SEQ=='14706',],aes(as.factor(NominalTier),x))+
  geom_col()


bcg$YrGrp<-ifelse(bcg$SampYr<=1999,"1989-1999",
                  ifelse (bcg$SampYr <=2009,"2000-2009","2010-2017"))

bcg<-bcg[,c(2,3,5,11)]
#bcgAvgSite<- aggregate(bcg$NominalTier,by=as.list(bcg[,c(1,2,4)]),FUN=function(x) c(mean=mean(x),n=length(x)))
bcgAvgSite<- aggregate(bcg$NominalTier,by=as.list(bcg[,c(1,2,4)]),FUN=mean)

for (i in 1:length(sites)){
  site<-(bcgAvgSite[bcgAvgSite$STA_SEQ==sites[i],])
  plot<- ggplot(site,aes(YrGrp,x))+
            geom_col()+
            labs(title=paste(site$Station_Name,site$STA_SEQ),x=NULL,y="BCG Avg Score")+
            ylim(0,6)
  ggsave(paste0(site$STA_SEQ,"_",site$Station_Name,".jpg"),device="jpeg",width = 4, height = 4,)
}

bcgYrG1<-bcgAvgSite[bcgAvgSite$YrGrp=="1989-1999",]
bcgYrG2<-bcgAvgSite[bcgAvgSite$YrGrp=="2000-2009",]
bcgYrG3<-bcgAvgSite[bcgAvgSite$YrGrp=="2010-2017",]

bcgT1T2<-merge(bcgYrG1,bcgYrG3,by=c("STA_SEQ","Station_Name"),all.x=TRUE)
colnames(bcgT1T2)[3:6]<-c("T1","BCG1","T2","BCG2")
bcgT1T2$Delta<-bcgT1T2$BCG1-bcgT1T2$BCG2
bcgT1T2$Trend<-as.factor(ifelse(bcgT1T2$Delta<0,"Neg",ifelse(bcgT1T2$Delta>0,"Pos","None")))

#####Change in Delta BCG between decades###
ggplot(bcgT1T2,aes(as.factor(STA_SEQ),Delta,fill=Trend))+
  geom_col(color="black")+
  ylim(-2,2)+
  labs(x="",y="Change",title=paste("Average Change in BCG Score between",bcgT1T2$T1,"and",bcgT1T2$T2))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0),legend.position = "none")

####Plot 1 site all years#####
bcgplotsite<-bcg[bcg$STA_SEQ==16119,]
ggplot(bcgplotsite,aes(SampYr,NominalTier))+
  geom_col()+
  labs(title=paste(bcgplotsite$Station_Name,bcgplotsite$STA_SEQ),x=NULL,y="BCG Score")+
  ylim(0,6)

#########SD BCG#############################################################################
bcgYrDiff<-data.frame(STA_SEQ=integer(),SampYr1=integer(),
                      bcgYr1=double(),SampYr2=integer(),
                      bcgYr2=double())#Create empty dataframe to store combinations

for (i in 1:length(sites)){
  bcgSite<-bcg[bcg$STA_SEQ == sites[i],c(2,4:5)]
  comb<-as.data.frame(t(combn(bcgSite$SampYr,2)))#create a dataframe of all possible combinations
  bcgSiteDiff<-merge(bcgSite,comb,by.x="SampYr",by.y="V1",all.y=TRUE)
  bcgSiteDiff<-merge(bcgSite,bcgSiteDiff,by.x="SampYr",by.y="V2",all.y=TRUE)
  bcgSiteDiff<-bcgSiteDiff[,c(2,1,3,4,6)]
  colnames(bcgSiteDiff)<-c("STA_SEQ","SampYr1","bcgYr1","SampYr2","bcgYr2")
  bcgYrDiff<-rbind(bcgYrDiff,bcgSiteDiff)
}



bcgYrDiff$bcgDiff<-abs(bcgYrDiff$bcgYr1-bcgYrDiff$bcgYr2)
bcgYrDiff$yrDiff<-abs(bcgYrDiff$SampYr1-bcgYrDiff$SampYr2)

bcg5Yr<-bcgYrDiff[bcgYrDiff$yrDiff==1,]
sd(bcg5Yr$bcgDiff)/mean(bcg5Yr$bcgDiff)


##Plot of typical differences in BCG at same site - taken 1 year apart##
typDiff<-(dim(bcg5Yr[bcg5Yr$bcgDiff<=1,])[1])/dim(bcg5Yr)[1]
ggplot(bcg5Yr,aes(bcgDiff))+
  geom_bar(aes(y = (..count..)/sum(..count..)*100))+
  ylim(0,100)+
  labs(y="Percent",x="Difference in BCG",title="Change in BCG for Samples Taken 1 Year Apart at the Same Site")


bcgDiff<-bcgYrDiff[,c(1,6)]
aggregate(bcgDiff$bcgDiff,by=list(STA_SEQ=bcgDiff$STA_SEQ),FUN=mean)

bcgYrDiff<-merge(bcgYrDiff,sampleflow,by.x=c("STA_SEQ","SampYr1"),by.y=c("STA_SEQ","SampYr"),all.x=TRUE)
colnames(bcgYrDiff)[8]<-"flow1"
bcgYrDiff<-merge(bcgYrDiff,sampleflow,by.x=c("STA_SEQ","SampYr2"),by.y=c("STA_SEQ","SampYr"),all.x=TRUE)
colnames(bcgYrDiff)[9]<-"flow2"
bcgYrDiff$flowDiff<-bcgYrDiff$flow1-bcgYrDiff$flow2

dim(bcgYrDiff[bcgYrDiff$bcgDiff>=1,])[1]/dim(bcgYrDiff)[1]


ggplot(bcgYrDiff[bcgYrDiff$bcgDiff<0,],aes(flowDiff,bcgDiff))+
  geom_point()



#########Pct 2 BCG#############################################################################

bcgPct<-data.frame(STA_SEQ=integer(),Station_Name=factor(),Pct2=numeric(),Years=integer(),Years2=integer(),AvgMem2=integer())

for (i in 1:length(sites)){
  bcgSite<-bcg[bcg$STA_SEQ == sites[i],c(2:5)]
  Years<-dim(bcgSite)[1]
  Years2<-dim(bcgSite[bcgSite$NominalTier<=2,])[1]
  Pct2<-Years2/Years
  bcg2<-bcg[bcg$NominalTier==2,]
  bcg2Site<-bcg2[bcg2$STA_SEQ==sites[i],]
  AvgMem2<-mean(bcg2Site$NominalMem)
  siteID<-unique(bcgSite$STA_SEQ)
  sitename<-unique(bcgSite$Station_Name)
  bcgSitePct<-data.frame(siteID,sitename,Pct2,Years,Years2,AvgMem2)
  colnames(bcgSitePct)<-c("STA_SEQ","Station_Name","Pct2","Years","Years2","AvgMem2")
  bcgPct<-rbind(bcgPct,bcgSitePct)
}

bcgPct<-bcgPct[order(bcgPct$Pct2),]
bcgPct$C2<-ifelse(bcgPct$Pct2>=0.5,1,0)
aggregate(bcgPct$AvgMem2)

##########Frag Forest Change################################################
##############################################################################
frag90<-read.csv("data/1990_ForestFragSum_RefSiteCatch.csv",header=TRUE)
frag90<-frag90[,6:14]
fragSum90<-aggregate(.~Str_Drain, frag90, sum)
fragSum90$CoreSum<-(fragSum90$X4+fragSum90$X5+fragSum90$X6)/fragSum90$Sum
fragSum90$FragSum<-(fragSum90$X1+fragSum90$X2+fragSum90$X3)/fragSum90$Sum
fragSum90$OtherSum<-(fragSum90$X0)/fragSum90$Sum
fragSum90<- fragSum90[,c(1,10:12)]

frag15<-read.csv("data/2015_ForestFragSum_RefSiteCatch.csv",header=TRUE)
frag15<-frag15[,6:14]
fragSum15<-aggregate(.~Str_Drain, frag15, sum)
fragSum15$CoreSum<-(fragSum15$X4+fragSum15$X5+fragSum15$X6)/fragSum15$Sum
fragSum15$FragSum<-(fragSum15$X1+fragSum15$X2+fragSum15$X3)/fragSum15$Sum
fragSum15$OtherSum<-(fragSum15$X0)/fragSum15$Sum
fragSum15<- fragSum15[,c(1,10:12)]

fragSum<-merge(fragSum90,fragSum15,by="Str_Drain")
fragSum$fragSumCoreDiff<-fragSum[,2]-fragSum[,5]
fragSum$fragSumFragDiff<-fragSum$FragSum.x-fragSum$FragSum.y
fragSum$fragSumOtherDiff<-fragSum$OtherSum.x-fragSum$OtherSum.y



