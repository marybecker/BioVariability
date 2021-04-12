setwd('C:/Users/deepuser/Documents/Projects/ProgramDev/BioVariability/Taxa_Changes_Landscape_Chem')

# Load libraries
library(ggplot2)

TP <- read.csv('data/TotalPhosphorus_MultiYearBCG2_Sites.csv',header=TRUE)
colnames(TP)[1] <-'staSeq'



ggplot(TP_agg,aes(as.factor(Grp),x,fill=as.factor(Grp)))+
  geom_boxplot(colour='#cccccc')+
  scale_y_log10()

LowTP_Grp1 <- data.frame(Pct = dim(TP[TP$Grp==1&TP$Avg_TP<=0.01,])[1]/dim(TP[TP$Grp==1,])[1],
                         YrGrp = '2000s')
LowTP_Grp2<- data.frame(Pct = dim(TP[TP$Grp==2&TP$Avg_TP<=0.01,])[1]/dim(TP[TP$Grp==2,])[1],
                        YrGrp = '2010s')
LowTP <- rbind(LowTP_Grp1,LowTP_Grp2)

LowTPPct<- ggplot(LowTP,aes(as.factor(YrGrp),Pct*100))+
  geom_col(fill="#74c476",alpha=c(0.6,0.9))+
  labs(y='Percent of Sites <= 0.01 mg/L Total Phosphorus')+
  theme(axis.title.x=element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey"))

ggsave(plot=LowTPPct,'plots/lowTPPct.png',width=5,height=5,units="in")

TP_reshp <- reshape(TP, idvar='staSeq',timevar='Grp',direction='wide')
colnames(TP_reshp)[c(1,2,3,5)] <- c('staSeq','locationName','TP_2000s','TP_2010s')
TP_diff <- TP_reshp[,c(1,2,3,5)]
TP_diff$TPDiff <- TP_diff$TP_2000s - TP_diff$TP_2010s

dim(TP_diff[TP_diff$TPDiff<0,])[1]/dim(TP_diff)[1]

mean(TP_diff[TP_diff$TPDiff<0,5]))

TP_diff[TP_diff$TPDiff<=-0.0,5]

wfpath <- 'C:/Users/deepuser/Documents/Projects/GISPrjs/2020/CWA_Pres/data/'
write.csv(TP_diff,paste0(wfpath,'TPDiff.csv'),row.names = FALSE)



 