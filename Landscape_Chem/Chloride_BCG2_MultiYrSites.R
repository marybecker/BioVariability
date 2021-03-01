# Set working directory
setwd('C:/Users/deepuser/Documents/Projects/ProgramDev/BioVariability/Landscape_Chem')

# Load ggplot
library(ggplot2)

# Load raw data
chloride_cond <- read.csv('data/Chloride_SpCond_MultiYrBCG2Sites.csv',header=TRUE)
colnames(chloride_cond)[1] <-'staSeq'
chloride_cond$Grp <- ifelse(chloride_cond$YearGrp<5,1,2)

# Set loop parameters
parameter <- c('Chloride','SpCond')
label <- c('Chloride(mg/L)','Conductivity ms/cm')

for (i in 1:length(parameter)){
  # Separate chloride and spCond
  chem_parameter <- chloride_cond[chloride_cond$ChemParameter==parameter[i],]
  
  # Aggregate the data and prep for plotting
  c_agg <- aggregate(chem_parameter$value,by=as.list(chem_parameter[,c('staSeq','Grp')]),FUN=mean)
  c_agg <- reshape(c_agg, idvar='staSeq',timevar='Grp',direction='wide')
  c_agg <- c_agg[complete.cases(c_agg),]
  c_agg$increase <- ifelse(c_agg$x.2>c_agg$x.1,1,0)
  c_agg$diff <- c_agg$x.2 - c_agg$x.1
  
  
  c_multiYr <- merge(chloride,c_agg,by='staSeq')
  c_multiYr <- aggregate(c_multiYr$value,by=as.list(c_multiYr[,c('staSeq','Grp')]),FUN=mean)
  
  
  # Set plot variables
  color <- c("#7bccc4","#a8ddb5")
  
  # Boxplot
  c_group_plot <- ggplot(c_multiYr,aes(as.factor(Grp),x))+
                  geom_boxplot(fill=color[i],alpha=c(0.6,0.9))+
                  labs(y=label[i])+
                  ylim(0,round(max(c_multiYr[,3])+0.05,1))+
                  scale_x_discrete(labels=c("2000s","2010s"))+
                  theme(axis.title.x=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "grey"))
  
  c_group_plot
  
  file <- paste0("plots/bcg2",parameter[i],".png")
  ggsave(plot=c_group_plot,file,width=5,height=5,units="in")
  

  #plot specific sites with long-term records
  sites <- c(14442,14314,14441,14470)
  c_site_plot <-  ggplot(chem_parameter[chem_parameter$staSeq%in%sites,],
                        aes(as.factor(Sample_Year),value,colour=locationName,group=locationName))+
                  geom_point()+
                  geom_line()+
                  labs(x='Year',y=label[i])+
                  theme(axis.title.x=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "grey"))
  
  file2 <- paste0("plots/bcg2Sites",parameter[i],".png")
  ggsave(plot=c_site_plot,file2,width=8,height=5,units="in")
}




                          
                      
                          


