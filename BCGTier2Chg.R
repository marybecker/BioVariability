setwd("C:/Users/deepuser/Documents/Projects/ProgramDev/BioVariability")

library(ggplot2)
library(vegan)
library(rgdal)
library(stringr)

tier<-2
attr<-4

BCG<-read.csv(paste0("data/bugBCG",tier,"_MutliYr_02112020.csv"),header=TRUE)
BCG$YrGrp<-ifelse(BCG$SampYr< 2000,1,ifelse(BCG$SampYr< 2010,3,2))
taxa<-read.csv(paste0("data/taxaData_Tier",tier,".csv"),header=TRUE)


#taxa_agg<-taxa[which(taxa$BCG_Attribute==2|taxa$BCG_Attribute==3),]
taxa_agg<-taxa[which(taxa$BCG_Attribute==4),]
taxa_agg<-aggregate(taxa_agg$ABUNDANCE,by=as.list(taxa_agg[,c("STA_SEQ","SampYr")]),FUN=sum)
BCG<-merge(BCG,taxa_agg,by=c("STA_SEQ","SampYr"),all.x=TRUE)
BCG$x[is.na(BCG$x)]<-0
BCG$RelAbund<-BCG$x/BCG$SumOfABUNDANCE

###################################################################################################
TPColors=c("#a6611a","#404040","#018571","#49006a")

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

#############SUMMARY STATS############################################################3
cln<-10

#Summary stats of differences in samples taken within 5 Years apart
nbase<-dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]#n combinations of samples within 5 Years apart
statbase<-summary(BCGYrDiff[BCGYrDiff$YrDiff<=5,cln])
quantbase<-quantile(BCGYrDiff[BCGYrDiff$YrDiff<=5,cln],c(0.05,0.95))
stdevbase<-sd(BCGYrDiff[BCGYrDiff$YrDiff<=5,cln])

n20<-dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]#n combinations of samples 20 Years apart
stat20<-summary(BCGYrDiff[BCGYrDiff$YrDiff>=20,cln])
quant20<-quantile(BCGYrDiff[BCGYrDiff$YrDiff>=20,cln],c(0.05,0.95))
stdev20<-sd(BCGYrDiff[BCGYrDiff$YrDiff>=20,cln])

n10<-dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]#n combinations of samples 10 Years apart
stat10<-summary(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,cln])
quant10<-quantile(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,cln],c(0.05,0.95))
stdev10<-sd(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,cln])

n5<-dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]#n combinations of samples 5-10 Years apart
stat5<-summary(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,cln])
quant5<-quantile(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,cln],c(0.05,0.95))
stdev5<-sd(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,cln])

BCGStat <-as.data.frame(rbind(statbase,stat20,stat10,stat5))
BCGStat$TimeP<-c("5 Yrs","20-30 Yrs","10-20 Yrs","5-10 Yrs")
BCGStat$N<-c(nbase,n20,n10,n5)
quant<-rbind(quantbase,quant20,quant10,quant5)
stdev<-rbind(stdevbase,stdev20,stdev10,stdev5)
BCGStat<-cbind(BCGStat,quant,stdev)

#write.csv(BCGStat,paste0("bugPlots/BCGStatDiff",tier,"Attr",attr,".csv"),row.names=FALSE)
write.csv(BCGStat,paste0("bugPlots/BCGStatDiff",tier,".csv"),row.names=FALSE)

####BCGAttrMets##################################################################################################
DescBase<- dim(BCGYrDiff[BCGYrDiff$BugDiff>(0)&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
MedDescBase<- dim(BCGYrDiff[BCGYrDiff$BugDiff>((stdevbase))&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
BCGDBase<- dim(BCGYrDiff[BCGYrDiff$BCGDiff<(0)&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
BCGIBase<- dim(BCGYrDiff[BCGYrDiff$BCGDiff>=1&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
BCGSBase<- dim(BCGYrDiff[BCGYrDiff$BCGDiff==0&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]


Desc20<- dim(BCGYrDiff[BCGYrDiff$BugDiff>(0)&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]
MedDesc20<- dim(BCGYrDiff[BCGYrDiff$BugDiff>((stdevbase))&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]
BCGD20<- dim(BCGYrDiff[BCGYrDiff$BCGDiff<(0)&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]
BCGI20<- dim(BCGYrDiff[BCGYrDiff$BCGDiff>=1&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]
BCGS20<- dim(BCGYrDiff[BCGYrDiff$BCGDiff==0&BCGYrDiff$YrDiff>=20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=20,])[1]


Desc10<- dim(BCGYrDiff[BCGYrDiff$BugDiff>(0)&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]
MedDesc10<- dim(BCGYrDiff[BCGYrDiff$BugDiff>((stdevbase))&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]
BCGD10<- dim(BCGYrDiff[BCGYrDiff$BCGDiff<0&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]
BCGI10<- dim(BCGYrDiff[BCGYrDiff$BCGDiff>=1&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]
BCGS10<- dim(BCGYrDiff[BCGYrDiff$BCGDiff==0&BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,])[1]


Desc5<- dim(BCGYrDiff[BCGYrDiff$BugDiff>(0)&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]
MedDesc5<- dim(BCGYrDiff[BCGYrDiff$BugDiff>((stdevbase))&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]
BCGD5<- dim(BCGYrDiff[BCGYrDiff$BCGDiff<0&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]
BCGI5<- dim(BCGYrDiff[BCGYrDiff$BCGDiff>=1&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]
BCGS5<- dim(BCGYrDiff[BCGYrDiff$BCGDiff==0&BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,])[1]


BCGMet <- data.frame(TimeP=factor(c("20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs"),levels=c("20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs")),
                       Desc=c(Desc20,Desc10,Desc5,DescBase),
                       MedDesc=c(MedDesc20,MedDesc10,MedDesc5,MedDescBase),
                       BCGD=c(BCGD20,BCGD10,BCGD5,BCGDBase),
                       BCGI=c(BCGI20,BCGI10,BCGI5,BCGIBase),
                       BCGS=c(BCGS20,BCGS10,BCGS5,BCGSBase))

names(TPColors)<-BCGMet$TimeP

DescP<-ggplot(BCGMet,aes(TimeP,Desc))+
  geom_col(aes(fill=TimeP,alpha=1))+
  labs(y="Proportion Increasing")+
  geom_text(aes(label=round(Desc,2)),position=position_dodge(width=0.9),vjust=-0.25)+
  ylim(0,1)+
  scale_fill_manual(values=TPColors)+
  theme(axis.title.x=element_blank(),legend.position="none",
        panel.background = element_rect(fill = "white", colour = "grey"))

MedDescP<-ggplot(BCGMet,aes(TimeP,MedDesc))+
  geom_col(aes(fill=TimeP,alpha=1))+
  labs(y="Proportion Increasing Greater Than One Stdev")+
  geom_text(aes(label=round(MedDesc,2)),position=position_dodge(width=0.9),vjust=-0.25)+
  ylim(0,1)+
  scale_fill_manual(values=TPColors)+
  theme(axis.title.x=element_blank(),legend.position="none",
        panel.background = element_rect(fill = "white", colour = "grey"))

ggsave(plot=DescP,paste0("bugPlots/PctTaxaAttChg",tier,"Att",attr,".jpg"),width=5,height=5,units="in")
ggsave(plot=MedDescP,paste0("bugPlots/PctTaxaAttChgStdDEV",tier,"Att",attr,".jpg"),width=5,height=5,units="in")

# MedDescP<- ggplot(BCGMet,aes(TimeP,MedDesc))+
#   geom_col(aes(fill=TimeP,alpha=0.8))+
#   labs(y="Percent Decreasing Greater Than Median")+
#   #ylim(0,1)+
#   scale_fill_manual(values=TPColors)+
#   theme(axis.title.x=element_blank(),legend.position="none")
# 
# BCGDP<-ggplot(BCGMet,aes(TimeP,BCGI))+
#   geom_col(aes(fill=TimeP,alpha=0.8))+
#   labs(y="Percent Change of BCG Increase")+
#   #ylim(0,1)+
#   scale_fill_manual(values=TPColors)+
#   theme(axis.title.x=element_blank(),legend.position="none")


DiffHistBase<-  ggplot(BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BugDiff,y=(..count..)/sum(..count..)))+
  geom_histogram(binwidth=0.05,fill=TPColors[4],alpha=0.7)+
  labs(x="Difference in Taxa Rel Abundance Between Two Years",y="Proportion of Samples")+
  xlim(-0.6,0.6)+
  ylim(0,0.6)+
  theme(panel.background = element_rect(fill = "white", colour = "grey"))

DiffHist20<-  ggplot()+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff>=20,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[1],alpha=0.9)+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in Taxa Rel Abundance Between Two Samples",y="Proportion of Samples")+
  xlim(-0.6,0.6)+
  ylim(0,0.6)+
  theme(panel.background = element_rect(fill = "white", colour = "grey"))

DiffHist10<-  ggplot()+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[2],alpha=0.9)+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in Taxa Rel Abundance Between Two Samples",y="Proportion of Samples")+
  xlim(-0.6,0.6)+
  ylim(0,0.6)+
  theme(panel.background = element_rect(fill = "white", colour = "grey"))

DiffHist5<- ggplot()+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[3],alpha=0.9)+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BugDiff,y=(..count..)/sum(..count..)),binwidth=0.05,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in Taxa Rel Abundance Between Two Samples",y="Proportion of Samples")+
  xlim(-0.6,0.6)+
  ylim(0,0.6)+
  theme(panel.background = element_rect(fill = "white", colour = "grey"))

##histogram plots
ggsave(plot=DiffHistBase,paste0("bugPlots/DiffHistBaseTier",tier,"Att",attr,".jpg"),width=5,height=5,units="in")
ggsave(plot=DiffHist20,paste0("bugPlots/DiffHist20Tier",tier,"Att",attr,".jpg"),width=5,height=5,units="in")
ggsave(plot=DiffHist10,paste0("bugPlots/DiffHist10Tier",tier,"Att",attr,".jpg"),width=5,height=5,units="in")
ggsave(plot=DiffHist5,paste0("bugPlots/DiffHist5Tier",tier,"Att",attr,".jpg"),width=5,height=5,units="in")

#Pct of samples that do not change between combinations of samples
dim(BCGYrDiff[BCGYrDiff$BugDiff==0&BCGYrDiff$YrDiff<=5|BCGYrDiff$BugDiff==2&BCGYrDiff$YrDiff<=5,])[1]/dim(BCGYrDiff[BCGYrDiff$YrDiff<=5,])[1]
dim(BCGYrDiff[BCGYrDiff$BugDiff==0|BCGYrDiff$BugDiff==2,])[1]/dim(BCGYrDiff)[1]

#Stat diff between distributions

BCGYrDiff$Grp<-ifelse(BCGYrDiff$YrDiff<=5,"Base",
                     ifelse(BCGYrDiff$YrDiff>=20,"Yr20",
                            ifelse(BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,"Yr10","Yr5")))

wilcox.test(BugDiff~Grp, data=BCGYrDiff[BCGYrDiff$Grp=="Base"|BCGYrDiff$Grp=="Yr5",],conf.int=TRUE)

###SlopeGraph##############################################################################

BCGMet2 <-data.frame(TimeP=factor(c("20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs",
                                       "20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs",
                                       "20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs"),
                                     levels=c("20-30 Yrs","10-20 Yrs","5-10 Yrs","5 Yrs")),
                        Lab=c("Decrease","Decrease","Decrease","Decrease",
                              "Increase","Increase","Increase","Increase",
                              "Stable", "Stable", "Stable", "Stable"),
                        Met=c(BCGD20,BCGD10,BCGD5,BCGDBase,BCGI20,BCGI10,BCGI5,BCGIBase,BCGS20,BCGS10,BCGS5,BCGSBase))

# ggplot(BCGMet2,aes(x=TimeP,y=Met,fill=Lab))+
#   geom_bar(position="fill",stat="identity")+
#   scale_y_continuous(labels=scales::percent_format())+
#   scale_fill_manual(values=c("#a6611a","#404040","#018571"))+
#   coord_flip()+
#   annotate("text",x=c(4.5,4.5,4.5),y=c(0.1,0.77,0.93),label=c("Stable","Increase","Decrease"))+
#   theme(legend.position = "none",panel.background = element_rect(fill = "white", colour = "white"),
#         axis.title=element_blank(),text=element_text(size=12,family="Sans"))

BCGmet2<- ggplot(BCGMet2)+
            geom_line(aes(TimeP,Met*100,group=Lab,colour=Lab),size=2)+
            geom_point(aes(TimeP,Met*100,group=Lab),colour="white",size=6)+
            geom_text(data=BCGMet2[BCGMet2$TimeP=="20-30 Yrs",],
                      aes(x=0,y=Met*100,label=c("Increase Tier","Decrease Tier","No Change"),colour=Lab),hjust=0)+
            geom_text(data=BCGMet2[BCGMet2$Lab=="Decrease"|
                                        BCGMet2$Lab=="Stable"|BCGMet2$Lab=="Increase",],aes(x=TimeP,y=Met*100,label=round(Met*100,0)),vjust=0)+
            #geom_text(data=BCGMet2[BCGMet2$Lab=="Increase",],aes(x=TimeP,y=Met*100,label=round(Met*100,0)),vjust=1.5)+
            scale_x_discrete(position = "top")+ 
            scale_color_manual(values=c("#a6611a","#404040","#018571"))+
            labs(x="Number of Years Between Samples\n")+
            lims(y=c(0,80))+
            theme(legend.position = "none",panel.background = element_rect(fill = "white", colour = "white"),
                  axis.title.y=element_blank(),text=element_text(size=16,family="Sans"),
                  axis.text.x=element_text(colour="black"),axis.ticks.length = unit(.5, "cm"),
                  axis.ticks.y=element_blank(),axis.text.y=element_blank())

ggsave(plot=BCGmet2,paste0("bugPlots/BCGmet2Tier",tier,".jpg"),width=8,height=5,units="in")

##########HISTOGRAM BCG###################################################################################
###########################################################################################################

DiffHistBaseBCG<-  ggplot(BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BCGDiff,y=(..count..)/sum(..count..)))+
  geom_histogram(binwidth=1,fill=TPColors[4],alpha=0.7)+
  labs(x="Difference in BCG Between Two Years",y="Proportion of Samples")+
  xlim(-3,3)+
  ylim(0,1)+
  theme(panel.background = element_rect(fill = "white", colour = "grey"))

DiffHist20BCG<-  ggplot()+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff>=20,],aes(x=BCGDiff,y=(..count..)/sum(..count..)),binwidth=1,fill=TPColors[1],alpha=0.9)+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BCGDiff,y=(..count..)/sum(..count..)),binwidth=1,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in Taxa Rel Abundance Between Two Samples",y="Proportion of Samples")+
  xlim(-3,3)+
  ylim(0,1)+
  theme(panel.background = element_rect(fill = "white", colour = "grey"))

DiffHist10BCG<-  ggplot()+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff>=10&BCGYrDiff$YrDiff<20,],aes(x=BCGDiff,y=(..count..)/sum(..count..)),binwidth=1,fill=TPColors[2],alpha=0.9)+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BCGDiff,y=(..count..)/sum(..count..)),binwidth=1,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in Taxa Rel Abundance Between Two Samples",y="Proportion of Samples")+
  xlim(-3,3)+
  ylim(0,1)+
  theme(panel.background = element_rect(fill = "white", colour = "grey"))

DiffHist5BCG<- ggplot()+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff>5&BCGYrDiff$YrDiff<10,],aes(x=BCGDiff,y=(..count..)/sum(..count..)),binwidth=1,fill=TPColors[3],alpha=0.9)+
  geom_histogram(data=BCGYrDiff[BCGYrDiff$YrDiff<=5,],aes(x=BCGDiff,y=(..count..)/sum(..count..)),binwidth=1,fill=TPColors[4],alpha=0.4)+
  labs(x="Difference in Taxa Rel Abundance Between Two Samples",y="Proportion of Samples")+
  xlim(-3,3)+
  ylim(0,1)+
  theme(panel.background = element_rect(fill = "white", colour = "grey"))

##histogram plots
ggsave(plot=DiffHistBaseBCG,paste0("bugPlots/DiffHistBaseTierBCG",tier,".jpg"),width=5,height=5,units="in")
ggsave(plot=DiffHist20BCG,paste0("bugPlots/DiffHist20TierBCG",tier,".jpg"),width=5,height=5,units="in")
ggsave(plot=DiffHist10BCG,paste0("bugPlots/DiffHist10TierBCG",tier,".jpg"),width=5,height=5,units="in")
ggsave(plot=DiffHist5BCG,paste0("bugPlots/DiffHist5TierBCG",tier,".jpg"),width=5,height=5,units="in")



Time<-BCGYrDiff$YrDiff<=5
Criteria1<-BCGYrDiff$BCGDiff<(-1)
Criteria2<-BCGYrDiff$BCGDiff>1


dim(BCGYrDiff[Time&Criteria1|Time&Criteria2,])[1]/dim(BCGYrDiff[Time,])[1]



