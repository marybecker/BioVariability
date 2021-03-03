# Set working directory
setwd('C:/Users/deepuser/Documents/Projects/ProgramDev/BioVariability/Landscape_Chem')

# Load libraries
library(ggplot2)

# Load raw data
taxa <- read.csv('data/bcg2Site_AllYears_Taxa.csv',header=TRUE)
taxa$YrGrp <- ifelse(taxa$YearGrp<3,1990,ifelse(taxa$YearGrp>4,2010,2000))
taxa$RegionalTempOptimaClass<- ifelse(taxa$RegionalTempOptimaClass == '','None', 
                                      taxa$RegionalTempOptimaClass)

# Manipulate raw data for boxplots.  Aggregate by Year than by decade.
temp_agg <- aggregate(taxa$RelAbund,
                     by=as.list(taxa[,c('staSeq','Sample_Year','YrGrp','RegionalTempOptimaClass')]),
                     FUN=sum)
temp_agg <- aggregate(temp_agg$x,
                      by=as.list(temp_agg[,c('staSeq','YrGrp','RegionalTempOptimaClass')]),
                      FUN=mean)

# Set parameters for boxplots
temp_class <- unique(temp_agg$RegionalTempOptimaClass)
print(temp_class)
bp_colors <- c("#045a8d","#a6bddb","#d9d9d9","#fcbba1")

# Run boxplot for each temperature class
for(i in 1:length(temp_class)){
  temp_class_agg <- temp_agg[temp_agg$RegionalTempOptimaClass == temp_class[i],]
  temp_class_agg <- temp_class_agg[,c(1,2,4)]
  temp_class_reshape <- reshape(temp_class_agg, idvar='staSeq',timevar='YrGrp',direction='wide')
  temp_class_reshape <- temp_class_reshape[complete.cases(temp_class_reshape),]
  temp_class_reshape <- reshape(temp_class_reshape,direction="long",idvar="staSeq",
                                varying=2:ncol(temp_class_reshape),sep=".")
  
  tempOpt_bp  <-  ggplot(temp_class_reshape,aes(as.factor(time),temp_class_reshape[,3]))+
                  geom_boxplot(fill=bp_colors[i],alpha=c(0.9,0.7,0.5))+
                  ylim(0,round(max(temp_class_agg[,3])+0.05,1))+
                  labs(y=paste("Relative Abundance",temp_class[i],"Water Taxa"))+
                  scale_x_discrete(labels=c("1990s","2000s","2010s"))+
                  theme(axis.title.x=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "grey"))
                  # theme(axis.title.x=element_blank(),
                  #       panel.background = element_rect(fill = "white", colour = "grey"))
  
  tempOpt_bp
  
  file <- paste0("plots/tempOpt_",temp_class[i],".png")
  ggsave(file,width=5,height=5,units="in")
}
