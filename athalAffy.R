library("ArrayExpress", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library(affy)
library(dplyr)
library(gcrma)
library(arrayQualityMetrics)
library(readr)
library(ggbiplot)
#query
#ATsets = queryAE(species = "Arabidopsis+thaliana")
#keep ones with raw data
#ATsetsfiltered<-ATsets %>% filter(tolower(Species)=="arabidopsis thaliana") %>% filter(tolower(Raw)=="yes")
#read exps list downloaded from website
AT_AE_Experiments_181227_014344 <- read_delim("ArrayExpress-Experiments-181227-194246.tsv", 
                                                     "\t", escape_double = FALSE, trim_ws = TRUE)
ATsetsfiltered<-AT_AE_Experiments_181227_014344 %>% filter(tolower(Organism)=="arabidopsis thaliana") %>% filter(`Raw Data`!="Data is not available") %>% filter(Assays >4)

smallATset<-head(ATsetsfiltered)

#download all .cel files
failed<-c()
celData<-data.frame()

idlist<-smallATset$Accession
for(id in idlist){
  print(id)
  atSamp = getAE(id, type = "raw")
  thisRawList<-atSamp$rawFiles
  #if there are no .cel files
  celInd<-grep("\\.CEL$",thisRawList)
  if(length(celInd)<1){
    print("Err")
    failed<-c(failed,id)
  }else{
    celData<-bind_rows(celData,data.frame(Exp=rep(id,length(thisRawList)),Cel=thisRawList))
    
    #read data
    #affy.data = ReadAffy(filenames = thisRawList)
    #apply gcrma
    #eset <- gcrma(affy.data)
    #get metadata
    #thisPD<-pData(eset)
    #get expression
    #thisExp<-exprs(eset)
  }
}

#read all cel files downloaded
files = list.files(pattern = "\\.CEL$", full.names = TRUE)
#files = list.files("celfiles2/", full.names = TRUE)
affy.data = ReadAffy(filenames = files)

#apply gcrma
eset <- gcrma(affy.data)
thisPD<-pData(eset)
thisExp<-exprs(eset)
dim(thisExp)
#log 2 normalize
thisExplog<-log2(thisExp)
dim(thisExplog)

#apply PCA to identify batch affects
tData<-t(thisExp) #changes dimention ???
dim(tData)

tDatalog<-t(thisExplog)
dim(tDatalog)
#filter to remove 0 variance cols (genes)
tDatalogF <- tDatalog[,apply(tDatalog, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
esetPCA<-prcomp(tDatalogF, scale. = TRUE)


ggbiplot(esetPCA, obs.scale = 1, var.scale = 1,
         groups = factor(celData$Exp), ellipse = TRUE, circle = TRUE,varname.size=0,var.axes = F) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')




#clean files






data(wine)
wine.pca <- prcomp(wine, scale. = TRUE)
ggbiplot(wine.pca, obs.scale = 1, var.scale = 1,
         groups = wine.class, ellipse = TRUE, circle = TRUE,varname.size=0) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
