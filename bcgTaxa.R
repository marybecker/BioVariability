setwd('/home/mkozlak/Documents/Projects/GitHub/BioVariability')

library(ggplot2)
library(vegan)
library(rgdal)
library(stringr)

##Read in the taxa, mastertaxa and site data
taxa<-read.csv("data/taxaData_5Yr_Tier2.csv",header=TRUE,stringsAsFactors = FALSE)
taxa$BCG_Attribute<-ifelse(is.na(taxa$BCG_Attribute),0,taxa$BCG_Attribute)
taxa<-taxa[taxa$BCG_Attribute!=0,]##Subset only taxa with a BCG Attribute
mastertaxa<-read.csv("data/masterTaxalist.csv",header=TRUE,stringsAsFactors = FALSE)
sites<-read.csv("data/bcg2_5YrSiteLoc.csv",header=TRUE)
sites<-sites[,c(1,2,4,5)]

##Define Taxa OTU.  Use Genus level. If 'Na' use Family Level
taxa$OTU<-taxa$GENUS
for (i in 1:dim(taxa)[1]){
  if(taxa$GENUS[i]=='Na'&taxa$FAMILY[i]=='Na'){
    taxa$OTU[i]<-taxa$TaxonNameCurrent[i]
  }
  if(taxa$FAMILY[i]!='Na'&taxa$GENUS[i]=='Na'){
    taxa$OTU[i]<-taxa$FAMILY[i]
  }
}

##Merge Taxa List with Master Taxa List
taxaM <- merge(taxa,mastertaxa,by.x="TaxonNameCurrent",by.y="finalID")
taxaM <- taxaM[taxaM$GENUS.x!='Na',] #only get taxa id down to the Genus level with BCG Attribute assigned

##Get and format data for phisViz
PhyloTree<- unique(taxaM[c("PHYLUM","CLASS","ORDER","FAMILY.y","GENUS.y","BCG_Attribute.y")])
colnames(PhyloTree)<- c("PHYLUM","CLASS","ORDER","FAMILY","Taxon_name","Category")
write.csv(PhyloTree,"phylo_tree.csv",row.names=FALSE)

sampData<-aggregate(taxaM$ABUNDANCE,by=as.list(taxaM[,c("STA_SEQ","Sample_Date","GENUS.y")]),FUN=sum)
colnames(sampData)<-c("SID","Collection_Date","Taxon_name","RelAbund")
sampData$dateLen<-str_length(sampData$Collection_Date)
sampData$Year<-str_sub(sampData$Collection_Date,sampData$dateLen-3,sampData$dateLen)
sampData$Month<-0
sampData$Day<-0

for (i in 1:dim(sampData)[1]){
  sampData$Day[i]<-str_sub(sampData$Collection_Date[i],
                             str_locate(sampData$Collection_Date[i],"/")[1]+1,
                             sampData$dateLen[i]-5)
}

for (i in 1:dim(sampData)[1]){
  sampData$Month[i]<-str_sub(sampData$Collection_Date[i],1,
                             str_locate(sampData$Collection_Date[i],"/")[1]-1)
}

sampData$Month<- ifelse(str_length(sampData$Month)==2,sampData$Month,paste0("0",sampData$Month))
sampData$Day<- ifelse(str_length(sampData$Day)==2,sampData$Day,paste0("0",sampData$Day))
sampData$Collection_Date<-paste0(sampData$Year,"/",sampData$Month,"/",sampData$Day)
sampData<-sampData[,1:4]

write.csv(sampData,"sample_data.csv",row.names=FALSE,quote=FALSE)

#Transform coordinates to numeric
sites$YLat  <- as.numeric(sites$YLat)
sites$XLong  <- as.numeric(sites$XLong)
sites.SP  <- SpatialPointsDataFrame(sites[,c(4,3)],
                                      sites[,-c(4,3)])
proj4string(sites.SP) <- CRS("+proj=utm +zone=18 +datum=WGS84") 
#proj4string(dataMap.SP) <- CRS("+init=epsg:4326") #WGS 84

str(sites.SP) # Now is class SpatialPointsDataFrame

#Write as geojson
writeOGR(sites.SP,"sites",layer="sites", driver='GeoJSON')


