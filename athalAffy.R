library("ArrayExpress")
library(affy)
library(dplyr)
library(gcrma)
library(readr)
library(ggbiplot)
library(sva)
library(ggplot2)
library(exploBATCH)
library(arrayQualityMetrics)

#query
#ATsets = queryAE(species = "Arabidopsis+thaliana")
#keep ones with raw data
#ATsetsfiltered<-ATsets %>% filter(tolower(Species)=="arabidopsis thaliana") %>% filter(tolower(Raw)=="yes")
#read exps list downloaded from website
AT_AE_Experiments_181227_014344 <- read_delim("ArrayExpress-Experiments-181227-194246.tsv", 
                                                     "\t", escape_double = FALSE, trim_ws = TRUE)
ATsetsfiltered<-AT_AE_Experiments_181227_014344 %>% filter(tolower(Organism)=="arabidopsis thaliana") %>% filter(`Raw Data`!="Data is not available") %>% filter(Assays >4)

ATsetsfiltered_10<-ATsetsfiltered%>% filter(Assays >9)

ATsetsfiltered_10_100<-ATsetsfiltered%>% filter(Assays >9) %>% filter(Assays <101)
sum(ATsetsfiltered_5_100$Assays)

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
  celInd<-grep("\\.CEL$",thisRawList,ignore.case = T)
  if(length(celInd)<1){
    print("Err")
    failed<-c(failed,id)
  }else{
    celData<-bind_rows(celData,data.frame(Exp=rep(id,length(thisRawList)),Cel=thisRawList))
  }
}

celData<-celData %>% mutate(batchId=as.numeric(factor(celData$Exp)))

#read all cel files downloaded
files = list.files(pattern = "\\.CEL$", full.names = TRUE,ignore.case = T)
affy.data = ReadAffy(filenames = files)

#apply gcrma
eset <- gcrma(affy.data)
thisPD<-pData(eset)
thisExp<-exprs(eset)
save.image("eset.RData")
dim(thisExp)
#log 2 normalize
thisExplog<-log2(thisExp)
dim(thisExplog)


#apply PCA to identify batch affects
tData<-t(thisExp) #changes dimention ???
dim(tData)

tDatalog<-t(thisExplog)
dim(tDatalog)
#


#filter to remove 0 variance cols (genes)
tDatalogF <- tDatalog[,apply(tDatalog, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
#see dist of corr vals
corrmat<-cor(tDatalogF)
dim(corrmat)
corrVec<-data.frame(cors=as.vector(corrmat[upper.tri(corrmat,diag = F)]))
sd<-head(corrVec,1000000)
p<-ggplot(data = corrVec,aes(x=cors))+ geom_histogram(bins = 10,colour="black", fill="white")

pdf("plot_nobatchcorr.pdf")
print(p)
dev.off()


esetPCA<-prcomp(tDatalogF, scale. = TRUE)


ggbiplot(esetPCA, obs.scale = 1, var.scale = 1,
         groups = factor(celData$Exp), ellipse = TRUE, circle = TRUE,varname.size=0,var.axes = F) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

#perform batch correction
batch = celData$Exp

# parametric adjustment
combat_edata1 = ComBat(dat=thisExp, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
dim(combat_edata1)
tDatalog<-t(log2(combat_edata1+1))
dim(tDatalog)
#filter to remove 0 variance cols (genes)
tDatalogF <- tDatalog[,apply(tDatalog, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
esetPCA<-prcomp(tDatalogF, scale. = TRUE)
corrmat<-cor(tDatalog)
# non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=thisExp, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
dim(combat_edata2)
tDatalog<-t(log2(combat_edata2))
dim(tDatalog)
#filter to remove 0 variance cols (genes)
tDatalogF <- tDatalog[,apply(tDatalog, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
esetPCA<-prcomp(tDatalogF, scale. = TRUE)



#clean files






data(wine)
wine.pca <- prcomp(wine, scale. = TRUE)
ggbiplot(wine.pca, obs.scale = 1, var.scale = 1,
         groups = wine.class, ellipse = TRUE, circle = TRUE,varname.size=0) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')




#explobatch
require(exploBATCH)
require(exploBATCHbreast)  # the package with breast cancer dataset
data(Breast)               # load the breast cancer data        
data(batchBreast)          # load the variable defining batches 


expBATCH(
  D=Breast,              # the breast cancer expression data matrix
  batchCL=batchBreast,   # the variable identifying the three batches
  Conf=NA,               # no biological variable of interest in this example
  mindim=2,              # the minimum number of pPCs
  maxdim=9,              # the maximum number of pPCs, we don't want this argument to be large,   
  # otherwise it will slow down exploBATCH.
  method="ppcca",        # use correctBATCH as well as ComBat to correct batch effect
  SDselect=2             # set SD at 2 to reduce computational time.
)

#Batch should be numeric

expBATCH(
  D=tData,              # the breast cancer expression data matrix
  batchCL=as.numeric(factor(celData$Exp)),   # the variable identifying the three batches
  Conf=NA,               # no biological variable of interest in this example
  mindim=2,              # the minimum number of pPCs
  maxdim=9,              # the maximum number of pPCs, we don't want this argument to be large,   
  # otherwise it will slow down exploBATCH.
  method="ppcca",        # use correctBATCH as well as ComBat to correct batch effect
  SDselect=0             # set SD at 2 to reduce computational time.
)




#filter good files
celData <- read_delim("C:/Users/mrbai/Desktop/celData.tsv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
toKeepExp<-ATsetsfiltered_10$Accession

celData_toKeep<-celData[celData$Exp %in% toKeepExp,]
length(unique(celData_toKeep$Exp))
length(unique(toKeepExp))
goodfiles <- read_csv("C:/Users/mrbai/Desktop/goodfiles.csv")
goodfiles$x<-substr(goodfiles$x,3,nchar(goodfiles$x))
celData_toKeep<-celData_toKeep[celData_toKeep$Cel %in% goodfiles$x,]
write.csv(celData_toKeep,"celDataReduced.csv",row.names = F)

celData_toKeep$Cel<-paste("./",celData_toKeep$Cel,sep = "")
write.csv(celData_toKeep$Cel,"goodFilesReduced.csv",row.names = F)

#zip files coressponding to exps
zipList <- read_csv("zipList.txt", col_names = FALSE)
#filter zip list which matches accession in ATsetsfiltered_10_100
zipList$X1<-substr(zipList$X1,3,nchar(zipList$X1))
zipList$X2<-sapply(strsplit(as.character(zipList$X1), ".raw"), "[[", 1)

filterZipList<- zipList[zipList$X2 %in% ATsetsfiltered_10_100$Accession,]
write_tsv(filterZipList,"zipList_10_100.tsv")
length(unique(filterZipList$X1))

#final celdata for 10_100 excluding failed to download files
celData_10_100<- celData[celData$Exp %in% unique(filterZipList$X2),]

