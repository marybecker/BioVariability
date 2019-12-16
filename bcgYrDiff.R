setwd('P:/Projects/GitHub_Prj/BioVariability')

library(ggplot2)

bcg<-read.csv("data/bug_bcg_10Yrs.csv",header=TRUE)
flow<-read.csv("data/findex_1989_2017.csv",header=TRUE)
flow<-flow[,c(1,14)]
colnames(flow)[1]<-"tripdate"
bcg<-merge(bcg,flow,by.x="tripdate",all.x=TRUE)
bcg<-bcg[bcg$STA_SEQ==14188|
           bcg$STA_SEQ==14314|
           bcg$STA_SEQ==14441|
           bcg$STA_SEQ==14442|
           bcg$STA_SEQ==14444|
           bcg$STA_SEQ==14450|
           bcg$STA_SEQ==14605|
           bcg$STA_SEQ==14706|
           bcg$STA_SEQ==14720,]##'Reference'Sites Only.  Comment out for all sites.
# bcg<-bcg[bcg$STA_SEQ!=14188&
#            bcg$STA_SEQ!=14314&
#            bcg$STA_SEQ!=14441&
#            bcg$STA_SEQ!=14442&
#            bcg$STA_SEQ!=14444&
#            bcg$STA_SEQ!=14450&
#            bcg$STA_SEQ!=14605&
#            bcg$STA_SEQ!=14706&
#            bcg$STA_SEQ!=14720,]##'Non-Reference'Sites Only.  Comment out for all sites.
sites<-unique(bcg$STA_SEQ)
sampleflow<-unique(bcg[c("STA_SEQ","SampYr","index")])



# bcgYrDiff<-data.frame(STA_SEQ=integer(),SampYr1=integer(),
#                       bcgYr1=double(),SampYr2=integer(),
#                       bcgYr2=double())#Create empty dataframe to store combinations
# 
# for (i in 1:length(sites)){
#   bcgSite<-bcg[bcg$STA_SEQ == sites[i],c(2,4:5)]
#   comb<-as.data.frame(t(combn(bcgSite$SampYr,2)))#create a dataframe of all possible combinations
#   bcgSiteDiff<-merge(bcgSite,comb,by.x="SampYr",by.y="V1",all.y=TRUE)
#   bcgSiteDiff<-merge(bcgSite,bcgSiteDiff,by.x="SampYr",by.y="V2",all.y=TRUE)
#   bcgSiteDiff<-bcgSiteDiff[,c(2,1,3,4,6)]
#   colnames(bcgSiteDiff)<-c("STA_SEQ","SampYr1","bcgYr1","SampYr2","bcgYr2")
#   bcgYrDiff<-rbind(bcgYrDiff,bcgSiteDiff)
# }
# 
# 
# 
# bcgYrDiff$bcgDiff<-(bcgYrDiff$bcgYr1-bcgYrDiff$bcgYr2)
# bcgYrDiff$yrDiff<-abs(bcgYrDiff$SampYr1-bcgYrDiff$SampYr2)
# bcgDiff<-bcgYrDiff[,c(1,6)]
# aggregate(bcgDiff$bcgDiff,by=list(STA_SEQ=bcgDiff$STA_SEQ),FUN=mean)
# 
# bcgYrDiff<-merge(bcgYrDiff,sampleflow,by.x=c("STA_SEQ","SampYr1"),by.y=c("STA_SEQ","SampYr"),all.x=TRUE)
# colnames(bcgYrDiff)[8]<-"flow1"
# bcgYrDiff<-merge(bcgYrDiff,sampleflow,by.x=c("STA_SEQ","SampYr2"),by.y=c("STA_SEQ","SampYr"),all.x=TRUE)
# colnames(bcgYrDiff)[9]<-"flow2"
# bcgYrDiff$flowDiff<-bcgYrDiff$flow1-bcgYrDiff$flow2
# 
# dim(bcgYrDiff[bcgYrDiff$bcgDiff>=1,])[1]/dim(bcgYrDiff)[1]
# 
# 
# ggplot(bcgYrDiff[bcgYrDiff$bcgDiff<0,],aes(flowDiff,bcgDiff))+
#   geom_point()
# 
# ggplot(bcg[bcg$STA_SEQ=='14720',],aes(SampYr,NominalTier))+
#   geom_point()
# 
bcg$YrGrp<-ifelse(bcg$SampYr<=1999,"1989-1999",
                  ifelse (bcg$SampYr <=2009,"2000 - 2009","2010-2017"))

bcg<-bcg[,c(2,3,5,11)]
bcgAvgSite<- aggregate(bcg$NominalTier,by=as.list(bcg[,c(1,2,4)]),FUN=function(x) c(mean=mean(x),n=length(x)))

for (i in 1:length(sites)){
  site<-(bcgAvgSite[bcgAvgSite$STA_SEQ==sites[i],])
  plot<- ggplot(site,aes(YrGrp,x.mean))+
            geom_col()+
            labs(title=paste(site$Station_Name,site$STA_SEQ),x=NULL,y="BCG Avg Score")+
            ylim(0,6)
  ggsave(paste0(site$STA_SEQ,"_",site$Station_Name,".jpg"),device="jpeg",width = 4, height = 4,)
}



