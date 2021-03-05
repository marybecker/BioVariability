# Set working directory
setwd('')

# Load libraries
library(ggplot2)
library(reshape2)

# Load data
bcg_land <- read.csv("data/bcg_BugSites_CoreForest_022521.csv",header=TRUE)
bcg_land <- bcg_land[bcg_land$MinOfNomin!=6,]
bcg2_land <- read.csv("data/bcg2BugSitesCoreForestDevAllYears.csv",header=TRUE)

# Manipulate data for boxplots (long to wide)
bcg2_CF <- reshape(bcg2_land[,c(1,10:16)],direction="long",idvar="STA_SEQ",
                   varying=2:ncol(bcg2_land[,c(1,10:16)]),sep="_")
bcg2_Dev <- reshape(bcg2_land[,c(1,17:23)],direction="long",idvar="STA_SEQ",
                   varying=2:ncol(bcg2_land[,c(1,17:23)]),sep="_")
bcg2_CFDiff <- reshape(bcg2_land[,c(1,24:29)],direction="long",idvar="STA_SEQ",
                   varying=2:ncol(bcg2_land[,c(1,24:29)]),sep="_")
bcg2_DevDiff <- reshape(bcg2_land[,c(1,30:35)],direction="long",idvar="STA_SEQ",
                    varying=2:ncol(bcg2_land[,c(1,30:35)]),sep="_")

# Boxplot of CF Across all BCG Tiers
bcg_colors<-c("#a6611a","#404040","#018571","#49006a")

bcg_landplt <-  ggplot(bcg_land,aes(as.factor(MinOfNomin),Statewide_CFPct_AllYears_CFPct15))+
            geom_boxplot(fill=bcg_colors,alpha=0.7)+
            ylim(0,1)+
            labs(y="% Core Forest",x="BCG Tier")+
            theme_bw()

ggsave(plot=bcg_landplt,"plots/bcg_CoreForest2015.png",width=5,height=5,units="in")

# Set the variables for the boxplot.  df org as STA_SEQ,time,x_value - from reshape
i<- 1
years<- c("1985","1995","2006","2015")
cf_colors<-c("#005824","#238b45","#41ae76","#66c2a4","#99d8c9","#ccece6","#edf8fb")
dev_colors<-c("#fee5d9","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#99000d")
list_of_bcg2_df<- list(bcg2_CFDiff,bcg2_DevDiff,bcg2_CF,bcg2_Dev)
bcg2_df <- list_of_bcg2_df[[i]] #choose the df to use
y_label_list<- list("% Core Forest Difference Since 1985",
                    "% Development Difference Since 1985",
                    "% Core Forest",
                    "% Development")
y_label<- y_label_list[[i]] #choose the y label
bplt_colors <- list(cf_colors,dev_colors,cf_colors,dev_colors)
colors <- bplt_colors[[i]]

# boxplot - change from 1985
bcg2_landplt <- ggplot(bcg2_df,aes(as.factor(time),bcg2_df[,3]*100))+
                geom_boxplot(fill=colors[1:length(unique(bcg2_df$time))],alpha=0.8)+
                labs(y=y_label[[1]],x="Year")+
                theme_bw()

file <- paste0("plots/bcg2",names(bcg2_df[3]),".png")
ggsave(plot=bcg2_landplt,file,width=5,height=5,units="in")



