# Set working directory
setwd('C:/Users/deepuser/Documents/Projects/ProgramDev/BioVariability/Taxa_Changes_Landscape_Chem')

# Load libraries
library(ggplot2)

CF_watershed <- read.csv('data/Statewide_CFPct_AllYears.csv',header=TRUE)
colnames(CF_watershed)<-c('HydroID','NextDown','SqMi','CFPct_1985','CFPct_1990','CFPct_1995',
                             'CFPct_2002','CFPct_2006','CFPct_2010','CFPct_2015')
CF_catchment <- read.csv('data/coreForestAllYears_catchments.csv',header=TRUE)
sites <- read.csv('data/BCG_Sites_1987_2019_HydroID.csv',header=TRUE)
BCG <- read.csv('data/BCG_Samples_1989_2019.csv',header=TRUE)
colnames(BCG)[1]<- 'STA_SEQ'
BCG <- BCG[,c(1,4,7)]
BCG <- merge(BCG,sites[,c(1,9)],by='STA_SEQ')

#sites_CF <- merge(sites,CF_watershed,by='HydroID')
#sites_CF_2_sqmi <- sites_CF[sites_CF$SqMi<=2,]

BCG$CFYear <- ifelse(BCG$SampleYear <1995,1990,
                      ifelse(BCG$SampleYear < 2000, 1995,
                             ifelse (BCG$SampleYear < 2005, 2002,
                                     ifelse(BCG$SampleYear < 2010, 2006,
                                            ifelse(BCG$SampleYear < 2015, 2010, 2015)))))

CFPct <- reshape(CF_watershed[,c(1,4:10)],direction="long",idvar="HydroID",
                   varying=2:ncol(CF_watershed[,c(1,4:10)]),sep="_")
colnames(CFPct)<-c("HydroID","CFYear","CFPct")

CFCatch <- reshape(CF_catchment[,c(1,5:11)],direction = "long", idvar="HydroID",
                              varying=2:ncol(CF_catchment[,c(1,5:11)]),sep="_")
colnames(CFCatch)<-c("HydroID","CFYear","CF")

BCG_CF <- merge(BCG,CFPct,by=c('HydroID','CFYear'))
BCG_CF <- merge(BCG_CF,CFCatch,by=c('HydroID','CFYear'))
BCG_CF <- merge(BCG_CF,CF_watershed[,c(1,3)],by='HydroID')
BCG_CF <- merge(BCG_CF,sites[,c(1,2,4,5)],by='STA_SEQ')



ggplot(BCG_CF,aes(as.factor(NominalTier),CFPct*100,fill=as.factor(NominalTier)))+
  geom_boxplot(colour='#cccccc')+
  xlab('BCG Tier')+
  ylab('% Core Forest Drainage Basin')+
  scale_fill_brewer(type = "seq",
                    palette = 4,
                    direction = -1)+ 
  theme(panel.background = element_rect(fill = '#252525', colour = '#969696'),
        plot.background = element_rect(fill = '#252525'),
        panel.grid = element_blank(),
        axis.text = element_text(colour='#cccccc',size=rel(1.1)),
        axis.title = element_text(color='#cccccc',size=rel(1.5)),
        legend.position = 'none')

# ggsave('C:/Users/deepuser/Documents/Projects/liss_legacyTN_biointegrity/images/BCGSites_CF_Drainage.png',
#        width=5,height=5,units="in")

ggplot(BCG_CF,aes(as.factor(NominalTier),CF*100,fill=as.factor(NominalTier)))+
  geom_boxplot(colour='#cccccc')+
  xlab('BCG Tier')+
  ylab('% Core Forest Catchment')+
  scale_fill_brewer(type = "seq",
                    palette = 4,
                    direction = -1)+ 
  theme(panel.background = element_rect(fill = '#252525', colour = '#969696'),
        plot.background = element_rect(fill = '#252525'),
        panel.grid = element_blank(),
        axis.text = element_text(colour='#cccccc',size=rel(1.1)),
        axis.title = element_text(color='#cccccc',size=rel(1.5)),
        legend.position = 'none')

# ggsave('C:/Users/deepuser/Documents/Projects/liss_legacyTN_biointegrity/images/BCGSites_CF_Catchment.png',
#        width=5,height=5,units="in")


ggplot(BCG_CF[BCG_CF$SqMi<6,],aes(CF,CFPct, colour=as.factor(NominalTier)))+
  geom_point(size = 4)+
  ylim(0,1)+
  xlim(0,1)+
  labs(colour = 'BCG Tier')+
  xlab('% Core Forest Catchment')+
  ylab('% Core Forest Drainage Basin')+
  geom_vline(xintercept = 0.5,colour='#969696')+
  geom_hline(yintercept = 0.5,colour='#969696')+
  scale_colour_brewer(type = "seq",
                      palette = 4,
                      direction = -1,
                      aesthetics = "colour") + 
  theme(panel.background = element_rect(fill = '#252525', colour = '#969696'),
            plot.background = element_rect(fill = '#252525'),
            panel.grid = element_blank(),
            axis.text = element_text(colour='#cccccc',size=rel(1.1)),
            axis.title = element_text(color='#cccccc',size=rel(1.1)),
            title = element_text(color = '#cccccc',size=rel(1.5)),
            legend.position = 'top',
            legend.background = element_rect(fill='#252525'),
            legend.key = element_rect(fill='#252525'),
            legend.text = element_text(colour = '#cccccc',size=rel(1.1)),
            legend.title = element_text(colour='#cccccc',size=rel(1.1)))

# ggsave('C:/Users/deepuser/Documents/Projects/liss_legacyTN_biointegrity/images/Small_BCGSites_CF_CFC.png',
#        width=8,height=5,units="in")


bcg_low <- BCG_CF[BCG_CF$CFPct>0.5 & 
                          BCG_CF$CF>0.5 &
                          BCG_CF$NominalTier >2 &
                          BCG_CF$SqMi < 6,]
bcg_low_sites <- unique(bcg_low[c('STA_SEQ','locationNa',
                                  'SqMi','HydroID','ylat','xlong')])
# write.csv(bcg_low_sites,'C:/Users/deepuser/Documents/Projects/GISPrjs/2021/LISS_legacy/data/BCG_Low_Sites_LT6SqMi.csv', row.names=FALSE)

ggplot(bcg_low,aes(CF,CFPct, colour=as.factor(NominalTier)))+
  geom_point(size = 4)+
  ylim(0,1)+
  xlim(0,1)+
  labs(colour = 'BCG Tier')+
  xlab('% Core Forest Catchment')+
  ylab('% Core Forest Drainage Basin')+
  geom_vline(xintercept = 0.5,colour='#969696')+
  geom_hline(yintercept = 0.5,colour='#969696')+
  scale_colour_manual(values = c('3' = '#43a2ca','4' = '#7bccc4', '5' = '#bae4bc','6'='#f0f9e8')) + 
  theme(panel.background = element_rect(fill = '#252525', colour = '#969696'),
        plot.background = element_rect(fill = '#252525'),
        panel.grid = element_blank(),
        axis.text = element_text(colour='#cccccc',size=rel(1.1)),
        axis.title = element_text(color='#cccccc',size=rel(1.1)),
        title = element_text(color = '#cccccc',size=rel(1.5)),
        legend.position = 'top',
        legend.background = element_rect(fill='#252525'),
        legend.key = element_rect(fill='#252525'),
        legend.text = element_text(colour = '#cccccc',size=rel(1.1)),
        legend.title = element_text(colour='#cccccc',size=rel(1.1)))

# ggsave('C:/Users/deepuser/Documents/Projects/liss_legacyTN_biointegrity/images/Small_LowBCGSites_CF_CFC.png',
#        width=8,height=5,units="in")


bcg_high <- merge(BCG_CF[BCG_CF$CFPct<0.5 & BCG_CF$CF<0.5 &BCG_CF$NominalTier == 2,],sites[,c(1:3)],by='STA_SEQ')
bcg_high_sites <- unique(bcg_high[c('STA_SEQ','locationNa','locationDe')])


#Cnt of Samples at Each Site Used in CF Calcs
BCG_Cnt <- aggregate(SampleYear~STA_SEQ,data=BCG_CF,FUN=length)
BCG_Cnt <- merge(BCG_Cnt,sites[,c(1:5,9)],by="STA_SEQ")
write.csv(BCG_Cnt,'C:/Users/deepuser/Documents/Projects/GISPrjs/2020/CWA_Pres/data/BCG_sites_1989_2019.csv', row.names=FALSE)
