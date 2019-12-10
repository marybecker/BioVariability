setwd('P:/Projects/GitHub_Prj/BioVariability')

library(ggplot2)

bcg<-read.csv("data/bug_bcg_10Yrs.csv",header=TRUE)
# bcg<-bcg[bcg$STA_SEQ==14188|
#            bcg$STA_SEQ==14314|
#            bcg$STA_SEQ==14441|
#            bcg$STA_SEQ==14442|
#            bcg$STA_SEQ==14444|
#            bcg$STA_SEQ==14450|
#            bcg$STA_SEQ==14605|
#            bcg$STA_SEQ==14706|
#            bcg$STA_SEQ==14720,]##'Reference'Sites Only.  Comment out for all sites.
bcg<-bcg[bcg$STA_SEQ!=14188&
           bcg$STA_SEQ!=14314&
           bcg$STA_SEQ!=14441&
           bcg$STA_SEQ!=14442&
           bcg$STA_SEQ!=14444&
           bcg$STA_SEQ!=14450&
           bcg$STA_SEQ!=14605&
           bcg$STA_SEQ!=14706&
           bcg$STA_SEQ!=14720,]##'Reference'Sites Only.  Comment out for all sites.
sites<-unique(bcg$STA_SEQ)


bcgYrDiff<-data.frame(STA_SEQ=integer(),SampYr1=integer(),
                      bcgYr1=double(),SampYr2=integer(),
                      bcgYr2=double())#Create empty dataframe to store combinations

for (i in 1:length(sites)){
  bcgSite<-bcg[bcg$STA_SEQ == sites[i],c(1,4:5)]
  comb<-as.data.frame(t(combn(bcgSite$SampYr,2)))#create a dataframe of all possible combinations
  bcgSiteDiff<-merge(bcgSite,comb,by.x="SampYr",by.y="V1",all.y=TRUE)
  bcgSiteDiff<-merge(bcgSite,bcgSiteDiff,by.x="SampYr",by.y="V2",all.y=TRUE)
  bcgSiteDiff<-bcgSiteDiff[,c(2,1,3,4,6)]
  colnames(bcgSiteDiff)<-c("STA_SEQ","SampYr1","bcgYr1","SampYr2","bcgYr2")
  bcgYrDiff<-rbind(bcgYrDiff,bcgSiteDiff)
}

bcgYrDiff$bcgDiff<-abs(bcgYrDiff$bcgYr1-bcgYrDiff$bcgYr2)
bcgYrDiff$yrDiff<-abs(bcgYrDiff$SampYr1-bcgYrDiff$SampYr2)
bcgDiff<-bcgYrDiff[,c(1,6)]
aggregate(bcgDiff$bcgDiff,by=list(STA_SEQ=bcgDiff$STA_SEQ),FUN=mean)

dim(bcgYrDiff[bcgYrDiff$bcgDiff>=1,])[1]/dim(bcgYrDiff)[1]


ggplot(bcgYrDiff,aes(yrDiff,bcgDiff))+
  geom_point()
