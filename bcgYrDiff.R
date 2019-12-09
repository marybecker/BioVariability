setwd('P:/Projects/GitHub_Prj/BioVariability')

bcg<-read.csv("data/bug_bcg_10Yrs.csv",header=TRUE)
sites<-unique(bcg$STA_SEQ)

bcgYrDiff<-data.frame()#Create empty dataframe

for i in length(sites){
  bcgSite<-bcg[bcg$STA_SEQ == sites[i],c(1,4:5)]
  comb<-as.data.frame(t(combn(bcgSite$SampYr,2)))
  bcgSiteDiff<-merge(bcgSite,comb,by.x="SampYr",by.y="V1",all.y=TRUE)
  bcgSiteDiff<-merge(bcgSite,bcgSiteDiff,by.x="SampYr",by.y="V2",all.y=TRUE)
  bcgSiteDiff<-bcgSiteDiff[,c(2,1,3,4,6)]
  colnames(bcgSiteDiff)<-c("STA_SEQ","SampYr1","bcgYr1","SampYr2","bcgYr2")
}
