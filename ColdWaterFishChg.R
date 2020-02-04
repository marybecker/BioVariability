setwd('P:/Projects/GitHub_Prj/BioVariability')

library(ggplot2)
library(vegan)
library(rgdal)
library(stringr)

CWsum<-read.csv("data/ColdWaterSites_ColdWaterFishSum.csv",header=TRUE)
ftaxa<-read.csv("data/FishSamplesColdWaterSites_012220.csv",header=TRUE,stringsAsFactors=FALSE)
ftaxa$OTU<-ifelse(ftaxa$IsStocked==TRUE,paste0(ftaxa$OTU,"-STK"),ftaxa$OTU)
pt<-read.csv("data/phylo_tree.csv",header=TRUE)
temp<-read.csv("data/CWSite_HOBOdata.csv",header=TRUE)

taxa<-ftaxa
taxa$SID<-paste0(taxa$STA_SEQ,"_",taxa$SampleYr)
taxa<-taxa[,c("SID","OTU","FishPer100M")]
taxa<-reshape(taxa,idvar="SID",timevar="OTU",direction="wide")
colnames(taxa)[2:dim(taxa)[2]]<-str_sub(colnames(taxa)[2:53],13,-1)
taxa[is.na(taxa)]<-0
row.names(taxa)<-taxa$SID
taxa<-taxa[,2:53]

taxa<-taxa[rowSums(taxa)>0,]
vegdist(taxa[c("16933_2015","16933_1998"),],"bray")

apply(taxa>0,1,sum)

###################################################################################################
TPColors=c("#e66101","#fdb863","#b2abd2","#49006a")

ggplot(CWsum,aes(YearGrp,MaxOfFishPer100M))+
  geom_boxplot(fill=TPColors[1:3])+
  scale_y_sqrt()+
  labs(y="Fish Per 100 M")+
  scale_x_discrete(labels=c("1990s","2000s","2010s"))+
  theme(axis.title.x=element_blank())

wilcox.test(MaxOfFishPer100M~YearGrp, data=CWsum[CWsum$YearGrp=="G1"|CWsum$YearGrp=="G3",],conf.int=TRUE)

##Identify all possible combinations of the same site in different years for comparison
CW<-CWsum[,1:3]
sites<-unique(CW$STA_SEQ)
CWYrDiff<-data.frame(STA_SEQ=integer(),SampYr1=integer(),
                      CWFishYr1=double(),SampYr2=integer(),
                      CWFishYr2=double())#Create empty dataframe to store combinations

for (i in 1:length(sites)){
  CWSite<-CW[CW$STA_SEQ == sites[i],]
  comb<-as.data.frame(t(combn(CWSite$SampleYear,2)))#create a dataframe of all possible combinations
  CWSiteDiff<-merge(CWSite,comb,by.x="SampleYear",by.y="V1",all.y=TRUE)
  CWSiteDiff<-merge(CWSite,CWSiteDiff,by.x="SampleYear",by.y="V2",all.y=TRUE)
  CWSiteDiff<-CWSiteDiff[,c(2,1,3,4,6)]
  colnames(CWSiteDiff)<-c("STA_SEQ","SampYr1","CWFishYr1","SampYr2","CWFishYr2")
  CWYrDiff<-rbind(CWYrDiff,CWSiteDiff)
}


CWYrDiff$FishDiff<-abs(CWYrDiff$CWFishYr1-CWYrDiff$CWFishYr2)#Difference in FishPer100M between samples
CWYrDiff$YrDiff<-abs(CWYrDiff$SampYr1-CWYrDiff$SampYr2)#Number of years between samples
CWYrDiff$y1Cold<-ifelse(CWYrDiff$CWFishYr1>=10,1,0)#If classified cold Samp1
CWYrDiff$y2Cold<-ifelse(CWYrDiff$CWFishYr2>=10,1,0)#If classified cold Samp1
CWYrDiff$ColdDiff<-CWYrDiff$y1Cold+CWYrDiff$y2Cold#Note if classification has changed
CWYrDiff$YrMax<-ifelse(CWYrDiff$SampYr1>CWYrDiff$SampYr2,CWYrDiff$SampYr1,CWYrDiff$SampYr2)#Most recent Yr
CWYrDiff$YrMaxDiff<-ifelse(CWYrDiff$SampYr1==CWYrDiff$YrMax,
                           CWYrDiff$CWFishYr1-CWYrDiff$CWFishYr2,
                           CWYrDiff$CWFishYr2-CWYrDiff$CWFishYr1)#Chg in FishPer100M most recent - earlier Samp
CWYrDiff$YrMaxColdDiff<-ifelse(CWYrDiff$SampYr1==CWYrDiff$YrMax,
                           CWYrDiff$y1Cold-CWYrDiff$y2Cold,
                           CWYrDiff$y2Cold-CWYrDiff$y1Cold)#Chg in Cold Class Most recent - earlier samp
CWYrDiff$CWFishMaxYr<-ifelse(CWYrDiff$SampYr1==CWYrDiff$YrMax,
               CWYrDiff$CWFishYr1,
               CWYrDiff$CWFishYr2)
CWYrDiff$CWFishMinYr<-ifelse(CWYrDiff$SampYr1==CWYrDiff$YrMax,
                             CWYrDiff$CWFishYr2,
                             CWYrDiff$CWFishYr1)
CWYrDiff$Dist<-NA

for (i in 1:dim(CWYrDiff)[1]){
  s1<-paste0(CWYrDiff[i,"STA_SEQ"],"_",CWYrDiff[i,"SampYr1"])
  s2<-paste0(CWYrDiff[i,"STA_SEQ"],"_",CWYrDiff[i,"SampYr2"])
  ifelse((any(row.names(taxa)==s1)==TRUE&any(row.names(taxa)==s2)==TRUE),
         CWYrDiff$Dist[i]<-vegdist(taxa[c(s1,s2),],"bray",binary=TRUE), CWYrDiff$Dist[i]<-1)
}

##ggplot(CWYrDiff,(aes(YrDiff,Dist)))+geom_point()

#Summary stats of differences in samples taken within 5 Years apart
dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]#n combinations of samples within 5 Years apart
summary(CWYrDiff[CWYrDiff$YrDiff<=5,6])
quantile(CWYrDiff[CWYrDiff$YrDiff<=5,6],c(0.05,0.95))

#Pct of samples that do not change cold/not cold category between combinations of samples
dim(CWYrDiff[CWYrDiff$ColdDiff==0&CWYrDiff$YrDiff<=5|CWYrDiff$ColdDiff==2&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]
dim(CWYrDiff[CWYrDiff$ColdDiff==0|CWYrDiff$ColdDiff==2,])[1]/dim(CWYrDiff)[1]

ggplot(CWYrDiff[CWYrDiff$YrDiff<=5,],aes(x=YrMaxDiff,y=(..count..)/sum(..count..)))+
  geom_histogram(binwidth=15,fill=TPColors[4],alpha=0.7)+
  labs(x="Difference in FishPer100M Between Two Years",y="Percent of Samples",
       title="Distribution of FishPer100M Differences Between Samples Taken 5 or Less Years Apart (n=409)")+
  xlim(-300,300)



dim(CWYrDiff[CWYrDiff$YrMaxDiff<(0)&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]
dim(CWYrDiff[CWYrDiff$YrMaxDiff<(-14)&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]

dim(CWYrDiff[CWYrDiff$YrMaxColdDiff==-1&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]
dim(CWYrDiff[CWYrDiff$YrMaxColdDiff==1&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]
dim(CWYrDiff[CWYrDiff$YrMaxColdDiff==0&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]


####################################################################################################


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

##Comparison of Distributions between samples taken within 5 year of each & Stream Survey to 2010 - 2018
ggplot()+
  geom_histogram(data=grpWide[complete.cases(grpWide[,5]),],aes(x=Y1Y3,y=(..count..)/sum(..count..)),binwidth=15,fill=TPColors[1],alpha=0.9)+
  geom_histogram(data=CWYrDiff[CWYrDiff$YrDiff<=5,],aes(x=YrMaxDiff,y=(..count..)/sum(..count..)),binwidth=15,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in FishPer100M Between Two Samples",y="Percent of Samples")+
  xlim(-300,300)

ggplot()+
  geom_histogram(data=grpWide[complete.cases(grpWide[,6]),],aes(x=Yr1Y2,y=(..count..)/sum(..count..)),binwidth=15,fill=TPColors[2],alpha=0.9)+
  geom_histogram(data=CWYrDiff[CWYrDiff$YrDiff<=5,],aes(x=YrMaxDiff,y=(..count..)/sum(..count..)),binwidth=15,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in FishPer100M Between Two Samples",y="Percent of Samples")+
  xlim(-300,300)
       #title="Distributions of FishPer100M Differences Between Years")

ggplot()+
  geom_histogram(data=grpWide[complete.cases(grpWide[,7]),],aes(x=Yr2Yr3,y=(..count..)/sum(..count..)),binwidth=15,fill=TPColors[3],alpha=0.9)+
  geom_histogram(data=CWYrDiff[CWYrDiff$YrDiff<=5,],aes(x=YrMaxDiff,y=(..count..)/sum(..count..)),binwidth=15,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in FishPer100M Between Two Samples",y="Percent of Samples")+
  xlim(-300,300)
       #title="Distributions of FishPer100M Differences Between Years")

summary(grpWide[complete.cases(grpWide[,5:7]),2])
summary(grpWide[complete.cases(grpWide[,5:7]),3])
summary(grpWide[complete.cases(grpWide[,5:7]),4])

AllYrGrp<-grpWide[complete.cases(grpWide[,5:7]),]
AllYrGrp$AllCold<-AllYrGrp$ColdYr1+AllYrGrp$ColdYr2+AllYrGrp$ColdYr3

                        
##Percent of samples that decreased in abundance from early time period to later time period
n<-(0)

t1t3Desc<-dim(grpWide[which(grpWide[,5]<n),])[1]/dim(grpWide[complete.cases(grpWide[,5]),])[1]
t1t2Desc<-dim(grpWide[which(grpWide[,6]<n),])[1]/dim(grpWide[complete.cases(grpWide[,6]),])[1]
t2t3Desc<-dim(grpWide[which(grpWide[,7]<n),])[1]/dim(grpWide[complete.cases(grpWide[,7]),])[1]
BLDesc<-dim(CWYrDiff[CWYrDiff$YrMaxDiff<n&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]

CWFishDesc<-data.frame(TimeP=factor(c("1990s-2010s","1990s-2000s","2000s-2010s","Baseline"),levels=c("1990s-2010s","1990s-2000s","2000s-2010s","Baseline")),
                       pctDesc=c(t1t3Desc,t1t2Desc,t2t3Desc,BLDesc))

names(TPColors)<-CWFishDesc$TimeP

ggplot(CWFishDesc,aes(TimeP,pctDesc))+
  geom_col(aes(fill=TimeP,alpha=0.6))+
  labs(y="Percent Decreasing")+
  ylim(0,1)+
  scale_fill_manual(values=TPColors)+
  theme(axis.title.x=element_blank(),legend.position="none")

##Percent of samples that decreased to zero CW Fish from early time period to later time period
t1t3Zero<-dim(grpWide[which(grpWide$Yr1>0&grpWide$Yr3==0),])[1]/dim(grpWide[complete.cases(grpWide[,5]),])[1]
t1t2Zero<-dim(grpWide[which(grpWide$Yr1>0&grpWide$Yr2==0),])[1]/dim(grpWide[complete.cases(grpWide[,6]),])[1]
t2t3Zero<-dim(grpWide[which(grpWide$Yr2>0&grpWide$Yr3==0),])[1]/dim(grpWide[complete.cases(grpWide[,7]),])[1]
BaseZero<-dim(CWYrDiff[which(CWYrDiff$CWFishMinYr>0&CWYrDiff$CWFishMaxYr==0&CWYrDiff$YrDiff<=5),])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]
CWFishZero<-data.frame(TimeP=factor(c("1990s-2010s","1990s-2000s","2000s-2010s","Baseline"),levels=c("1990s-2010s","1990s-2000s","2000s-2010s","Baseline")),
                       pctZero=c(t1t3Zero,t1t2Zero,t2t3Zero,BaseZero))

ggplot(CWFishZero,aes(TimeP,pctZero))+
  geom_col(aes(fill=TimeP,alpha=0.2))+
  labs(y="Percent Zero")+
  ylim(0,0.25)+
  scale_fill_manual(values=TPColors)+
  theme(axis.title.x=element_blank(),legend.position="none")


##Changes in Cold Water Category Across Different Time Periods
t1t3DescCold<-dim(grpWide[which(grpWide[,11]==-1),])[1]/dim(grpWide[complete.cases(grpWide[,11]),])[1]
t1t3IncCold<-dim(grpWide[which(grpWide[,11]==1),])[1]/dim(grpWide[complete.cases(grpWide[,11]),])[1]
t1t3Stable<-dim(grpWide[which(grpWide[,11]==0),])[1]/dim(grpWide[complete.cases(grpWide[,11]),])[1]
cwChgT1T3<-data.frame(chg=c("decreasing","increasing","stable"),
           pct=c(t1t3DescCold,t1t3IncCold,t1t3Stable))

chgP1<- ggplot(cwChgT1T3,aes(chg,pct))+
          geom_col(fill=TPColors[1])+
          ylim(0,1)+
          labs(y="Percent of Samples")+
          theme(axis.title.x=element_blank())

ggsave(plot=chgP1,"fishPlots/cwChgT1T3.jpg",width=4,height=4,units="in")

t1t2DescCold<-dim(grpWide[which(grpWide[,12]==-1),])[1]/dim(grpWide[complete.cases(grpWide[,12]),])[1]
t1t2IncCold<-dim(grpWide[which(grpWide[,12]==1),])[1]/dim(grpWide[complete.cases(grpWide[,12]),])[1]
t1t2Stable<-dim(grpWide[which(grpWide[,12]==0),])[1]/dim(grpWide[complete.cases(grpWide[,12]),])[1]
cwChgT1T2<-data.frame(chg=c("decreasing","increasing","stable"),
                      pct=c(t1t2DescCold,t1t2IncCold,t1t2Stable))

chgP2<-  ggplot(cwChgT1T2,aes(chg,pct))+
          geom_col(fill=TPColors[2])+
          ylim(0,1)+
          labs(y="Percent of Samples")+
          theme(axis.title.x=element_blank())

ggsave(plot=chgP2,"fishPlots/cwChgT1T2.jpg",width=4,height=4,units="in")

t2t3DescCold<-dim(grpWide[which(grpWide[,13]==-1),])[1]/dim(grpWide[complete.cases(grpWide[,13]),])[1]
t2t3IncCold<-dim(grpWide[which(grpWide[,13]==1),])[1]/dim(grpWide[complete.cases(grpWide[,13]),])[1]
t2t3Stable<-dim(grpWide[which(grpWide[,13]==0),])[1]/dim(grpWide[complete.cases(grpWide[,13]),])[1]
cwChgT2T3<-data.frame(chg=c("decreasing","increasing","stable"),
                      pct=c(t2t3DescCold,t2t3IncCold,t2t3Stable))
ggplot(cwChgT2T3,aes(chg,pct))+
  geom_col()

chgP3<- ggplot(cwChgT2T3,aes(chg,pct))+
          geom_col(fill=TPColors[3])+
          ylim(0,1)+
          labs(y="Percent of Samples")+
          theme(axis.title.x=element_blank())

ggsave(plot=chgP3,"fishPlots/cwChgT2T3.jpg",width=4,height=4,units="in")

cwChgBase<-data.frame(chg=c("decreasing","increasing","stable"),
                      pct=c(dim(CWYrDiff[CWYrDiff$YrMaxColdDiff==-1&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1],
                            dim(CWYrDiff[CWYrDiff$YrMaxColdDiff==1&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1],
                            dim(CWYrDiff[CWYrDiff$YrMaxColdDiff==0&CWYrDiff$YrDiff<=5,])[1]/dim(CWYrDiff[CWYrDiff$YrDiff<=5,])[1]))

chgP4<- ggplot(cwChgBase,aes(chg,pct))+
          geom_col(fill=TPColors[4])+
          ylim(0,1)+
          labs(y="Percent of Samples")+
          theme(axis.title.x=element_blank())

ggsave(plot=chgP4,"fishPlots/cwChgBase.jpg",width=4,height=4,units="in")


ChgPct<-rbind(cwChgT1T3,cwChgT1T2,cwChgT2T3,cwChgBase)
ChgPct$TP<-c("1990s-2010s","1990s-2010s","1990s-2010s","1990s-2000s","1990s-2000s",
             "1990s-2000s","2000s-2010s","2000s-2010s","2000s-2010s","Baseline","Baseline","Baseline")
ChgPct$TP<-factor(ChgPct$TP,levels=(c("1990s-2010s","1990s-2000s","2000s-2010s","Baseline")))


ggplot(ChgPct[ChgPct$chg=="decreasing",],aes(TP,pct))+
  geom_col(fill=TPColors)+
  ylim(0,0.4)+
  labs(y="Percent of Samples")+
  theme(axis.title.x=element_blank())





########################################################################################################
###Multi-Year Temperature data #########################################################################
########################################################################################################

##Identify all possible combinations of the same site in different years for comparison

tempCat<-temp[,c(1:3)]
tempSites<-unique(tempCat$SID)
TempYrDiff<-data.frame(STA_SEQ=integer(),SampYr1=integer(),
                     TempYr1=double(),SampYr2=integer(),
                     TempYr2=double())#Create empty dataframe to store combinations

for (i in 1:length(tempSites)){
  tempSite<-tempCat[tempCat$SID == tempSites[i],]
  comb<-as.data.frame(t(combn(tempSite$Year,2)))#create a dataframe of all possible combinations
  TempSiteDiff<-merge(tempSite,comb,by.x="Year",by.y="V1",all.y=TRUE)
  TempSiteDiff<-merge(tempSite,TempSiteDiff,by.x="Year",by.y="V2",all.y=TRUE)
  TempSiteDiff<-TempSiteDiff[,c(2,1,3,4,6)]
  colnames(TempSiteDiff)<-c("STA_SEQ","SampYr1","TempYr1","SampYr2","TempYr2")
  TempYrDiff<-rbind(TempYrDiff,TempSiteDiff)
}


TempYrDiff$TempDiff<-abs(TempYrDiff$TempYr1-TempYrDiff$TempYr2)#Difference in FishPer100M between samples
TempYrDiff$YrDiff<-abs(TempYrDiff$SampYr1-TempYrDiff$SampYr2)#Number of years between samples
TempYrDiff$y1Cold<-ifelse(TempYrDiff$TempYr1<18.29,1,0)#If classified cold Samp1
TempYrDiff$y2Cold<-ifelse(TempYrDiff$TempYr2<18.29,1,0)#If classified cold Samp1
TempYrDiff$ColdDiff<-TempYrDiff$y1Cold+TempYrDiff$y2Cold#Note if classification has changed
TempYrDiff$YrMax<-ifelse(TempYrDiff$SampYr1>TempYrDiff$SampYr2,TempYrDiff$SampYr1,TempYrDiff$SampYr2)#Most recent Yr
TempYrDiff$YrMaxDiff<-ifelse(TempYrDiff$SampYr1==TempYrDiff$YrMax,
                           TempYrDiff$TempYr1-TempYrDiff$TempYr2,
                           TempYrDiff$TempYr2-TempYrDiff$TempYr1)#Chg in FishPer100M most recent - earlier Samp
TempYrDiff$YrMaxColdDiff<-ifelse(TempYrDiff$SampYr1==TempYrDiff$YrMax,
                               TempYrDiff$y1Cold-TempYrDiff$y2Cold,
                               TempYrDiff$y2Cold-TempYrDiff$y1Cold)#Chg in Cold Class Most recent - earlier samp
TempYrDiff$TempMaxYr<-ifelse(TempYrDiff$SampYr1==TempYrDiff$YrMax,
                             TempYrDiff$TempYr1,
                             TempYrDiff$TempYr2)
TempYrDiff$TempMinYr<-ifelse(TempYrDiff$SampYr1==TempYrDiff$YrMax,
                             TempYrDiff$TempYr2,
                             TempYrDiff$TempYr1)

#Summary stats of differences in samples taken within 5 Years apart
dim(TempYrDiff[TempYrDiff$YrDiff<=5,])[1]#n combinations of samples within 5 Years apart
summary(TempYrDiff[TempYrDiff$YrDiff<=5,6])
quantile(TempYrDiff[TempYrDiff$YrDiff<=5,6],c(0.05,0.95))

#Pct of samples that remain cold between combinations of samples
dim(TempYrDiff[TempYrDiff$ColdDiff==0&TempYrDiff$YrDiff<=5|TempYrDiff$ColdDiff==2&TempYrDiff$YrDiff<=5,])[1]/dim(TempYrDiff[TempYrDiff$YrDiff<=5,])[1]
dim(TempYrDiff[TempYrDiff$ColdDiff==0|TempYrDiff$ColdDiff==2,])[1]/dim(TempYrDiff)[1]

dim(TempYrDiff[TempYrDiff$YrMaxColdDiff==(1)&TempYrDiff$YrDiff<=5,])[1]/dim(TempYrDiff[TempYrDiff$YrDiff<=5,])[1]

ggplot(TempYrDiff[TempYrDiff$YrDiff<=5,],aes(x=YrMaxDiff,y=(..count..)/sum(..count..)))+
  geom_histogram(binwidth=0.5,fill="#5ab4ac",alpha=0.7)+
  labs(x="Difference in Temp (Degree C) Between Two Years",y="Percent of Samples",
       title="Distribution of Temp Differences Between Samples Taken 5 or Less Years Apart (n=723)")



dim(TempYrDiff[TempYrDiff$YrMaxDiff<(0)&TempYrDiff$YrDiff<=5,])[1]/dim(TempYrDiff[TempYrDiff$YrDiff<=5,])[1]
dim(TempYrDiff[TempYrDiff$YrMaxDiff<(-1.2)&TempYrDiff$YrDiff<=5,])[1]/dim(TempYrDiff[TempYrDiff$YrDiff<=5,])[1]

########################################################################################################
###Prepare geospatial data##############################################################################
########################################################################################################

grpWide[complete.cases(grpWide[,5]),]


########################################################################################################
###Format data for PhisViz / Mapping####################################################################
########################################################################################################
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


###Prepare geospatial data##############################################################################
Y1Y3Sites<-merge(grpWide[complete.cases(grpWide[,5]),],sites,by="SID")

Y1Y3Sites$YLat  <- as.numeric(Y1Y3Sites$YLat)
Y1Y3Sites$XLong  <- as.numeric(Y1Y3Sites$XLong)
Y1Y3Sites.SP  <- SpatialPointsDataFrame(Y1Y3Sites[,c(16,15)],
                                        Y1Y3Sites[,-c(16,15)])
proj4string(Y1Y3Sites.SP) <- CRS("+proj=utm +zone=18 +datum=WGS84") 
#proj4string(dataMap.SP) <- CRS("+init=epsg:4326") #WGS 84

str(Y1Y3Sites.SP) # Now is class SpatialPointsDataFrame

#Write as geojson
writeOGR(Y1Y3Sites.SP,"sitesY1Y3fish",layer="Y1Y3Sites", driver='GeoJSON',overwrite_layer = TRUE)













