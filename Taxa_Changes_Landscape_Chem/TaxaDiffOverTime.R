# Set working directory
setwd('C:/Users/deepuser/Documents/Projects/ProgramDev/BioVariability/Taxa_Changes_Landscape_Chem')

# Load libraries
library(ggplot2)

bcgTier<- "2"

# Load raw taxa data
taxa <- read.csv(paste0('data/bcg',bcgTier,'Site_AllYears_Taxa.csv'),header=TRUE)
colnames(taxa)[1]<-'staSeq'

master_taxa_list <- read.csv('data/MasterTaxaList.csv',header=TRUE)
bcg_taxa_list <- master_taxa_list[,c('finalID','PHYLUM','CLASS','ORDER','FAMILY',
                                    'GENUS','BCG_Attribute')]
colnames(bcg_taxa_list)[1] <- 'TaxonNameCurrent'

# Load landcover data
bcg2_land <- read.csv('data/bcg2BugSitesCoreForestDevAllYears.csv',header=TRUE)
catchments_CF <- read.csv('data/coreForestAllYears_catchments.csv',header=TRUE)
bcg2_land_catchments <- merge(bcg2_land[,c('STA_SEQ','HydroID')],catchments_CF,by = "HydroID")

# Identify sites with multiple years of data
taxa_multiyr <- unique(taxa[c("staSeq","Sample_Year")])
taxa_multiyr <- aggregate(Sample_Year ~ staSeq, data= taxa_multiyr, FUN = length)
taxa_multiyr <- taxa_multiyr[taxa_multiyr$Sample_Year>1,]
taxa_multiyr <- merge(taxa, taxa_multiyr,by="staSeq")
samples <- unique(taxa_multiyr[c("staSeq","Sample_Year.x")])
colnames(samples) <- c("STA_SEQ","SampYr")

# Function that takes in dataframe as c("STA_SEQ","SampYr","RelAbund" - x)
# Returns all combinations of variable diffs for x across sample years
bug_year_comb <- function(taxa_data){
  sites <- unique(taxa_data$STA_SEQ)
  #Create empty dataframe to store combinations
  bug_year_diff<-data.frame(STA_SEQ=integer(),SampYr1=integer(),BugYr1=double(),
                            SampYr2=integer(),BugYr2=double())
  
  #create a dataframe of all possible combinations of samples across years
  for (i in 1:length(sites)){
    bug_site<-taxa_data[taxa_data$STA_SEQ == sites[i],]
    comb<-as.data.frame(t(combn(bug_site$SampYr,2)))
    bug_siteDiff<-merge(bug_site,comb,by.x="SampYr",by.y="V1",all.y=TRUE)
    bug_siteDiff<-merge(bug_site,bug_siteDiff,by.x="SampYr",by.y="V2",all.y=TRUE)
    bug_siteDiff<-bug_siteDiff[,c(2,1,3,4,6)]
    colnames(bug_siteDiff)<-c("STA_SEQ","SampYr1","BugYr1","SampYr2","BugYr2")
    bug_year_diff<-rbind(bug_year_diff,bug_siteDiff)
  }
  
  bug_year_diff$BugDiff<-bug_year_diff$BugYr1-bug_year_diff$BugYr2
  bug_year_diff$YrDiff<-bug_year_diff$SampYr1-bug_year_diff$SampYr2
  
  return(bug_year_diff)
}

##############################################################################################
#Plots Bug data Over Time
##############################################################################################

######RELATIVE ABUNDANCE######################################################################
##############################################################################################

#Set variables for loop
taxa_attributes <- list(taxa_multiyr[taxa_multiyr$BCG_Attribute<4,],
                   taxa_multiyr[taxa_multiyr$BCG_Attribute==4,],
                   taxa_multiyr[taxa_multiyr$BCG_Attribute>4,])
colors <- c("#008837","#f9711d","#B39DDB")
taxa_attribute_names <- c("Sensitive","Moderate","Tolerant")



#Run for each bcg taxa attribute
for(i in 1:length(taxa_attributes)){
  #Select bcg attribute taxa (sensitive, moderate or tolerant taxa)
  select_taxa <- taxa_attributes[[i]]
  
  # Aggregate selected bcg taxa attribute 
  bugc <- aggregate(select_taxa$RelAbund,
                    by=as.list(select_taxa[,c("staSeq","Sample_Year.x")]),
                    FUN=sum)
  colnames(bugc)<-c("STA_SEQ","SampYr","RelAbund")
  
  #Merge any samples that have 0% rel abund of bcg attribute taxa
  bugc <- merge(samples,bugc,by=c("STA_SEQ","SampYr"),all.x=TRUE)
  bugc[is.na(bugc)]<-0
  
  #call function
  bug_year_diff <- bug_year_comb(bugc)
  
  # Add Year Group
  bug_year_diff$YrGrp<-ifelse(bug_year_diff$YrDiff<=5,1,
                              ifelse(bug_year_diff$YrDiff>=20,4,
                                     ifelse(bug_year_diff$YrDiff>9&bug_year_diff$YrDiff<20,3,2)))
  
  # Run and save plots
  bug_year_diff_bp <- ggplot(bug_year_diff,aes(as.factor(YrGrp),BugDiff))+
    geom_boxplot(fill = colors[i])+
    labs(y=paste("Difference in Relative Abundance of",
                 taxa_attribute_names[i],"Taxa"),
         x="Number of Years Between Two Samples")+
    ylim(-0.5,0.5)+
    scale_x_discrete(labels=c("1-5","5-10","10-20","20-30"))+
    theme(panel.background = element_rect(fill = "white", colour = "grey"))
  
  file1 <- paste0("plots/taxadiffbp_bcg",bcgTier,taxa_attribute_names[i],".png")
  ggsave(plot=bug_year_diff_bp,file1,width=5,height=5,units="in")
  
  # Get the median value for each year group
  by_bug_year_diff <- by(bug_year_diff$BugDiff, bug_year_diff$YrGrp,median,simplify=FALSE)
  median_bug_year_diff <- as.data.frame(cbind(by_bug_year_diff))
  median_bug_year_diff$years <- rownames(median_bug_year_diff)
  
  bug_year_diff_bar <-  ggplot(median_bug_year_diff,aes(as.factor(years),by_bug_year_diff))+
    geom_col(fill = colors[i])+
    labs(y=paste("Median Difference in Relative Abundance of",
                 taxa_attribute_names[i],"Taxa"),
         x="Number of Years Between Two Samples")+
    coord_cartesian(ylim=c(-0.15,0.15))+
    scale_x_discrete(labels=c("1-5","5-10","10-20","20-30"))+
    theme(panel.background = element_rect(fill = "white", colour = "grey"))
  bug_year_diff_bar
  
  file2 <- paste0("plots/taxadiffmedbar_bcg",bcgTier,taxa_attribute_names[i],".png")
  ggsave(plot=bug_year_diff_bar,file2,width=5,height=5,units="in")
}

######COUNT OF TAXA###########################################################################
##############################################################################################

#Set variables for loop
cnt_bcg_taxa <- merge(taxa_multiyr,bcg_taxa_list, by=c('TaxonNameCurrent','BCG_Attribute'))
cnt_bcg_taxa <- unique(cnt_bcg_taxa[c('staSeq','Sample_Year.x','BCG_Attribute','PHYLUM',
                                      'CLASS','ORDER','FAMILY','GENUS')])
cnt_bcg_taxa$RelAbund <- 1
cnt_taxa <- aggregate(cnt_bcg_taxa$RelAbund,
                      by=as.list(cnt_bcg_taxa[,c("staSeq","Sample_Year.x")]),
                      FUN=sum)
pct_bcg_taxa <- merge(cnt_bcg_taxa,cnt_taxa,by=c('staSeq','Sample_Year.x'))
pct_bcg_taxa$RelAbund <- pct_bcg_taxa$RelAbund/pct_bcg_taxa$x

taxa_attributes_cnt <- list(pct_bcg_taxa[which(pct_bcg_taxa$BCG_Attribute<4),],
                            pct_bcg_taxa[which(pct_bcg_taxa$BCG_Attribute==4),],
                            pct_bcg_taxa[which(pct_bcg_taxa$BCG_Attribute>4),])

# taxa_attributes_cnt <- list(cnt_bcg_taxa[cnt_bcg_taxa$BCG_Attribute<4,],
#                         cnt_bcg_taxa[cnt_bcg_taxa$BCG_Attribute==4,],
#                         cnt_bcg_taxa[cnt_bcg_taxa$BCG_Attribute>4,])
colors <- c("#008837","#f9711d","#B39DDB")
taxa_attribute_names <- c("Sensitive","Moderate","Tolerant")


#Run for each bcg taxa attribute
for(i in 1:length(taxa_attributes_cnt)){
  #Select bcg attribute taxa (sensitive, moderate or tolerant taxa)
  select_taxa <- taxa_attributes_cnt[[i]]
  
  # Aggregate selected bcg taxa attribute 
  bugc <- aggregate(select_taxa$RelAbund,
                    by=as.list(select_taxa[,c("staSeq","Sample_Year.x")]),
                    FUN=sum)
  colnames(bugc)<-c("STA_SEQ","SampYr","RelAbund")
  
  #Merge any samples that have 0% rel abund of bcg attribute taxa
  bugc <- merge(samples,bugc,by=c("STA_SEQ","SampYr"),all.x=TRUE)
  bugc[is.na(bugc)]<-0
  
  #call function
  bug_year_diff <- bug_year_comb(bugc)
  
  # Add Year Group
  bug_year_diff$YrGrp<-ifelse(bug_year_diff$YrDiff<=5,1,
                              ifelse(bug_year_diff$YrDiff>=20,4,
                                     ifelse(bug_year_diff$YrDiff>9&bug_year_diff$YrDiff<20,3,2)))
  
  # Run and save plots
  bug_year_diff_bp <- ggplot(bug_year_diff,aes(as.factor(YrGrp),BugDiff))+
    geom_boxplot(fill = colors[i])+
    labs(y=paste("Difference in % Taxa of",
                 taxa_attribute_names[i],"Taxa"),
         x="Number of Years Between Two Samples")+
    #ylim(-0.5,0.5)+
    scale_x_discrete(labels=c("1-5","5-10","10-20","20-30"))+
    theme(panel.background = element_rect(fill = "white", colour = "grey"))
  
  bug_year_diff_bp
  
  file1 <- paste0("plots/pcttaxadiffbp_bcg",bcgTier,taxa_attribute_names[i],".png")
  ggsave(plot=bug_year_diff_bp,file1,width=5,height=5,units="in")
  
  # Get the median value for each year group
  by_bug_year_diff <- by(bug_year_diff$BugDiff, bug_year_diff$YrGrp,median,simplify=FALSE)
  median_bug_year_diff <- as.data.frame(cbind(by_bug_year_diff))
  median_bug_year_diff$years <- rownames(median_bug_year_diff)
  
  bug_year_diff_bar <-  ggplot(median_bug_year_diff,aes(as.factor(years),by_bug_year_diff))+
    geom_col(fill = colors[i])+
    labs(y=paste("Median Difference in % Taxa of",
                 taxa_attribute_names[i],"Taxa"),
         x="Number of Years Between Two Samples")+
    #coord_cartesian(ylim=c(-0.1,0.1))+
    scale_x_discrete(labels=c("1-5","5-10","10-20","20-30"))+
    theme(panel.background = element_rect(fill = "white", colour = "grey"))
  bug_year_diff_bar
  
  file2 <- paste0("plots/pcttaxadiffmedbar_bcg",bcgTier,taxa_attribute_names[i],".png")
  ggsave(plot=bug_year_diff_bar,file2,width=5,height=5,units="in")
}

bug_year_diff$CFYear1 <- ifelse(bug_year_diff$SampYr1 <1995,1990,
                                ifelse(bug_year_diff$SampYr1 < 2000, 1995,
                                       ifelse (bug_year_diff$SampYr1 < 2005, 2002,
                                               ifelse(bug_year_diff$SampYr1 < 2010, 2006,
                                                      ifelse(bug_year_diff$SampYr1 < 2015, 
                                                             2010, 2015)))))
bug_year_diff$CFYear2 <- ifelse(bug_year_diff$SampYr2 <1995,1990,
                                ifelse(bug_year_diff$SampYr2 < 2000, 1995,
                                       ifelse (bug_year_diff$SampYr2 < 2005, 2002,
                                               ifelse(bug_year_diff$SampYr2 < 2010, 2006,
                                                      ifelse(bug_year_diff$SampYr2 < 2015, 
                                                             2010, 2015)))))

bcg2_CF <- reshape(bcg2_land[,c(1,10:16)],direction="long",idvar="STA_SEQ",
                   varying=2:ncol(bcg2_land[,c(1,10:16)]),sep="_")
colnames(bcg2_CF)<-c("STA_SEQ","CFYear1","CFPct1")

taxa_CF <- merge(bug_year_diff,bcg2_CF,by=c("STA_SEQ","CFYear1"))

colnames(bcg2_CF)<-c("STA_SEQ","CFYear2","CFPct2")

taxa_CF <- merge(taxa_CF,bcg2_CF,by=c("STA_SEQ","CFYear2"))
taxa_CF$CFPctDiff <- taxa_CF$CFPct1 - taxa_CF$CFPct2

ggplot(taxa_CF[taxa_CF$YrDiff>15,],aes(CFPctDiff,BugDiff))+
  geom_point()+
  ylim(0,1)+
  xlim(0,1)+
  scale_colour_gradient(low = "white", high = "black")

##############################################################################################
#Merge Bug Data with Landuse
##############################################################################################


#Select bcg attribute taxa (sensitive, moderate or tolerant taxa)
taxa_relAbund <- taxa[,c('staSeq','Sample_Year','RelAbund','TaxonNameCurrent','BCG_Attribute')]
colnames(taxa_relAbund)<-c('staSeq','Date','RelAbund','TaxaName','Attribute')
select_taxa <- taxa_relAbund[which(taxa_relAbund$Attribute<4),]

# Aggregate selected bcg taxa attribute 
bugc <- aggregate(select_taxa$RelAbund,
                  by=as.list(select_taxa[,c("staSeq","Date")]),
                  FUN=sum)
colnames(bugc)<-c("STA_SEQ","SampYr","RelAbund")

#Merge any samples that have 0% rel abund of bcg attribute taxa
bugc <- merge(samples,bugc,by=c("STA_SEQ","SampYr"),all.x=TRUE)
bugc[is.na(bugc)]<-0

bugc$CFYear <- ifelse(bugc$SampYr <1995,1990,
                      ifelse(bugc$SampYr < 2000, 1995,
                             ifelse (bugc$SampYr < 2005, 2002,
                                     ifelse(bugc$SampYr < 2010, 2006,
                                            ifelse(bugc$SampYr < 2015, 2010, 2015)))))

bcg2_CF <- reshape(bcg2_land[,c(1,10:16)],direction="long",idvar="STA_SEQ",
                   varying=2:ncol(bcg2_land[,c(1,10:16)]),sep="_")
colnames(bcg2_CF)<-c("STA_SEQ","CFYear","CFPct")

bcg2_CF_catchments <- reshape(bcg2_land_catchments[,c(2,6:12)],direction = "long", idvar="STA_SEQ",
                              varying=2:ncol(bcg2_land_catchments[,c(2,6:12)]),sep="_")
colnames(bcg2_CF_catchments)<-c("STA_SEQ","CFYear","CF")



taxa_CF <- merge(bugc,bcg2_CF,by=c("STA_SEQ","CFYear"))
taxa_CF <- merge(taxa_CF,bcg2_CF_catchments,by=c("STA_SEQ","CFYear"))

ggplot(taxa_CF,aes(CF,CFPct,colour=RelAbund))+
  geom_point()+
  ylim(0,1)+
  xlim(0,1)+
  scale_colour_gradient(low = "white", high = "black")




