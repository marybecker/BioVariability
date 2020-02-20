setwd("P:/Projects/GitHub_Prj/BioVariability")

library(ggplot2)
library(vegan)
library(rgdal)
library(stringr)

BCG<-read.csv("data/bugBCG2_MutliYr_02112020.csv",header=TRUE)
BCG$YrGrp<-ifelse(BCG$SampYr< 2000,1,ifelse(BCG$SampYr< 2010,3,2))
taxa<-read.csv("data/taxaData_Tier2.csv",header=TRUE)

taxa_agg<-taxa[which(taxa$BCG_Attribute==2|taxa$BCG_Attribute==3),]
#taxa_agg<-taxa[which(taxa$BCG_Attribute>=4),]
taxa_agg<-aggregate(taxa_agg$ABUNDANCE,by=as.list(taxa_agg[,c("STA_SEQ","SampYr")]),FUN=sum)
BCG<-merge(BCG,taxa_agg,by=c("STA_SEQ","SampYr"),all.x=TRUE)
BCG$x[is.na(BCG$x)]<-0
BCG$RelAbund<-BCG$x/BCG$SumOfABUNDANCE

###################################################################################################
TPColors=c("#e66101","#fdb863","#636363","#49006a")

# BCGBP<- ggplot(BCG,aes(as.factor(YrGrp),RelAbund))+
#   geom_boxplot(fill=TPColors[1:3])+
#   scale_y_sqrt()+
#   labs(y="Rel Abund")+
#   scale_x_discrete(labels=c("1990s","2000s","2010s"))+
#   theme(axis.title.x=element_blank())


BugC<-BCG[,c(1,2,14)]
BCGC<-BCG[,c(1,2,6)]
sites<-unique(BCGC$STA_SEQ)
BugYrDiff<-data.frame(STA_SEQ=integer(),SampYr1=integer(),
                     BugYr1=double(),SampYr2=integer(),
                     BugYr2=double())#Create empty dataframe to store combinations

for (i in 1:length(sites)){
  bugSite<-BugC[BugC$STA_SEQ == sites[i],]
  comb<-as.data.frame(t(combn(bugSite$SampYr,2)))#create a dataframe of all possible combinations
  bugSiteDiff<-merge(bugSite,comb,by.x="SampYr",by.y="V1",all.y=TRUE)
  bugSiteDiff<-merge(bugSite,bugSiteDiff,by.x="SampYr",by.y="V2",all.y=TRUE)
  bugSiteDiff<-bugSiteDiff[,c(2,1,3,4,6)]
  colnames(bugSiteDiff)<-c("STA_SEQ","SampYr1","BugYr1","SampYr2","BugYr2")
  BugYrDiff<-rbind(BugYrDiff,bugSiteDiff)
}

BCGYrDiff<-data.frame(STA_SEQ=integer(),SampYr1=integer(),
                      BCGYr1=double(),SampYr2=integer(),
                      BCGYr2=double())#Create empty dataframe to store combinations

for (i in 1:length(sites)){
  BCGSite<-BCGC[BCGC$STA_SEQ == sites[i],]
  comb<-as.data.frame(t(combn(BCGSite$SampYr,2)))#create a dataframe of all possible combinations
  BCGSiteDiff<-merge(BCGSite,comb,by.x="SampYr",by.y="V1",all.y=TRUE)
  BCGSiteDiff<-merge(BCGSite,BCGSiteDiff,by.x="SampYr",by.y="V2",all.y=TRUE)
  BCGSiteDiff<-BCGSiteDiff[,c(2,1,3,4,6)]
  colnames(BCGSiteDiff)<-c("STA_SEQ","SampYr1","BCGYr1","SampYr2","BCGYr2")
  BCGYrDiff<-rbind(BCGYrDiff,BCGSiteDiff)
}

BCGYrDiff<-merge(BCGYrDiff,BugYrDiff,by=c("STA_SEQ","SampYr1","SampYr2"))


BCGYrDiff$BugDiff<-BCGYrDiff$BugYr1-BCGYrDiff$BugYr2
BCGYrDiff$YrDiff<-BCGYrDiff$SampYr1-BCGYrDiff$SampYr2
BCGYrDiff$BCGDiff<-BCGYrDiff$BCGYr1-BCGYrDiff$BCGYr2
BCGYrDiff$ABSBugDiff<-abs(BCGYrDiff$BugYr1-BCGYrDiff$BugYr2)

summary(BCGYrDiff[BCGYrDiff$YrDiff<=5,11])
####NEEDS WORK**************################################################
DescBase<- dim(BCGYrDiff[BCGYrDiff$BugDiff<(0)&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
MedDescBase<- dim(BCGYrDiff[BCGYrDiff$BugDiff<(-0.1)&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
BCGDBase<- dim(BCGYrDiff[BCGYrDiff$BCGDiff<(0)&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
BCGIBase<- dim(BCGYrDiff[BCGYrDiff$BCGDiff>1&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
BCGSBase<- dim(BCGYrDiff[BCGYrDiff$BCGDiff==0&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]


Desc20<- dim(BCGYrDiff[BCGYrDiff$BugDiff<(0)&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]
MedDesc20<- dim(BCGYrDiff[BCGYrDiff$BugDiff<(-0.1)&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]
BCGD20<- dim(BCGYrDiff[BCGYrDiff$BCGDiff<(0)&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]
BCGI20<- dim(BCGYrDiff[BCGYrDiff$BCGDiff>1&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]
BCGS20<- dim(BCGYrDiff[BCGYrDiff$BCGDiff==0&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]


Desc10<- dim(BCGYrDiff[BCGYrDiff$BugDiff<(0)&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]
MedDesc10<- dim(BCGYrDiff[BCGYrDiff$BugDiff<(-0.1)&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]
BCGD10<- dim(BCGYrDiff[BCGYrDiff$BCGDiff<0&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]
BCGI10<- dim(BCGYrDiff[BCGYrDiff$BCGDiff>1&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]
BCGS10<- dim(BCGYrDiff[BCGYrDiff$BCGDiff==0&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]


Desc5<- dim(BCGYrDiff[BCGYrDiff$BugDiff<(0)&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]
MedDesc5<- dim(BCGYrDiff[BCGYrDiff$BugDiff<(-0.1)&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]
BCGD5<- dim(BCGYrDiff[BCGYrDiff$BCGDiff<0&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]
BCGI5<- dim(BCGYrDiff[BCGYrDiff$BCGDiff>1&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]
BCGS5<- dim(BCGYrDiff[BCGYrDiff$BCGDiff==0&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]


BCGMet <- data.frame(TimeP=factor(c("20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs"),levels=c("20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs")),
                       Desc=c(Desc20,Desc10,Desc5,DescBase),
                       MedDesc=c(MedDesc20,MedDesc10,MedDesc5,MedDescBase),
                       BCGD=c(BCGD20,BCGD10,BCGD5,BCGDBase),
                       BCGI=c(BCGI20,BCGI10,BCGI5,BCGIBase),
                       BCGS=c(BCGS20,BCGS10,BCGS5,BCGSBase))

names(TPColors)<-BCGMet$TimeP

DescP<-ggplot(BCGMet,aes(TimeP,Desc))+
  geom_col(aes(fill=TimeP,alpha=0.8))+
  labs(y="Percent Decreasing")+
  #ylim(0,1)+
  scale_fill_manual(values=TPColors)+
  theme(axis.title.x=element_blank(),legend.position="none")

MedDescP<- ggplot(BCGMet,aes(TimeP,MedDesc))+
  geom_col(aes(fill=TimeP,alpha=0.8))+
  labs(y="Percent Decreasing Greater Than Median")+
  #ylim(0,1)+
  scale_fill_manual(values=TPColors)+
  theme(axis.title.x=element_blank(),legend.position="none")

BCGDP<-ggplot(BCGMet,aes(TimeP,BCGI))+
  geom_col(aes(fill=TimeP,alpha=0.8))+
  labs(y="Percent Change of BCG Increase")+
  #ylim(0,1)+
  scale_fill_manual(values=TPColors)+
  theme(axis.title.x=element_blank(),legend.position="none")


DiffHistBase<-  ggplot(BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BugDiff,y=(..count..)/sum(..count..)))+
  geom_histogram(binwidth=0.1,fill=TPColors[4],alpha=0.7)+
  labs(x="Difference in Taxa Rel Abundance Between Two Years",y="Percent of Samples")+
  xlim(-0.6,0.6)

DiffHist20<-  ggplot()+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff>=20,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[1],alpha=0.9)+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in Taxa Rel Abundance Between Two Samples",y="Percent of Samples")+
  xlim(-0.6,0.6)

DiffHist10<-  ggplot()+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[2],alpha=0.9)+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in Taxa Rel Abundance Between Two Samples",y="Percent of Samples")+
  xlim(-0.6,0.6)

DiffHist5<- ggplot()+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.1,fill=TPColors[3],alpha=0.9)+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.1,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in Taxa Rel Abundance Between Two Samples",y="Percent of Samples")+
  xlim(-0.6,0.6)

#Pct of samples that do not change cold/not cold category between combinations of samples
dim(BCGYrDiff[BCGYrDiff$BugDiff==0&BCGYrDiff$YrDiff<=5|BCGYrDiff$BugDiff==2&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
dim(BCGYrDiff[BCGYrDiff$BugDiff==0|BCGYrDiff$BugDiff==2,])[1]/dim(BCGYrDiff)[1]


###SlopeGraph##############################################################################

BCGMet2 <-data.frame(TimeP=factor(c("20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs",
                                       "20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs",
                                       "20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs"),
                                     levels=c("20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs")),
                        Lab=c("Decrease","Decrease","Decrease","Decrease",
                              "Increase","Increase","Increase","Increase",
                              "Stable", "Stable", "Stable", "Stable"),
                        Met=c(BCGD20,BCGD10,BCGD5,BCGDBase,BCGI20,BCGI10,BCGI5,BCGIBase,BCGS20,BCGS10,BCGS5,BCGSBase))

ggplot(BCGMet2,aes(x=TimeP,y=Met,fill=Lab))+
  geom_bar(position="fill",stat="identity")+
  scale_y_continuous(labels=scales::percent_format())+
  scale_fill_manual(values=c("#a6611a","#404040","#018571"))+
  coord_flip()+
  annotate("text",x=c(4.5,4.5,4.5),y=c(0.1,0.77,0.93),label=c("Stable","Increase","Decrease"))+
  theme(legend.position = "none",panel.background = element_rect(fill = "white", colour = "white"),
        axis.title=element_blank(),text=element_text(size=12,family="Sans"))

ggplot(BCGMet2)+
  geom_line(aes(TimeP,Met*100,group=Lab,colour=Lab),size=2)+
  geom_point(aes(TimeP,Met*100,group=Lab),colour="white",size=6)+
  geom_text(data=BCGMet2[BCGMet2$TimeP=="10-20 Yrs",],
            aes(x=0,y=Met*100,label=c("Cold to Not Cold","Not Cold to Cold","Stable"),colour=Lab),hjust=0)+
  geom_text(data=CWFishMet2[CWFishMet2$Lab=="Decrease"|
                              CWFishMet2$Lab=="Stable",],aes(x=TimeP,y=Met*100,label=round(Met*100,0)),vjust=-1)+
  geom_text(data=CWFishMet2[CWFishMet2$Lab=="Increase",],aes(x=TimeP,y=Met*100,label=round(Met*100,0)),vjust=1.5)+
  scale_x_discrete(position = "top")+ 
  scale_color_manual(values=c("#a6611a","#404040","#018571"))+
  labs(x="Number of Years Between Samples\n")+
  lims(y=c(0,80))+
  theme(legend.position = "none",panel.background = element_rect(fill = "white", colour = "white"),
        axis.title.y=element_blank(),text=element_text(size=16,family="Sans"),
        axis.text.x=element_text(colour="black"),axis.ticks.length = unit(.5, "cm"),
        axis.ticks.y=element_blank(),axis.text.y=element_blank())