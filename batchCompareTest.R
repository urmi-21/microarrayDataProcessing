library("ArrayExpress")
library(affy)
library(dplyr)
library(gcrma)
library(arrayQualityMetrics)
library(readr)
library(ggbiplot)
library(sva)
library(ggplot2)
library(exploBATCH)


#download E-MTAB-5610
atSamp5610 = getAE("E-MTAB-5610", type = "raw")


affy.data_5610 = ReadAffy(filenames = atSamp5610$rawFiles)
eset_5610 <- gcrma(affy.data_5610)
thisPD_5610<-pData(eset_5610)
thisExp_5610<-exprs(eset_5610)
boxplot(thisExp_5610)
boxplot(thisExp)
dim(thisExp_5610)
#log 2 normalize
thisExp_5610log<-log2(thisExp_5610)
dim(thisExp_5610log)
tDatalog<-t(thisExp_5610log)
dim(tDatalog)
#filter to remove 0 variance cols (genes)
tDatalogF <- tDatalog[,apply(tDatalog, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
PCA<-prcomp(tDatalogF, scale. = TRUE)
batch<-c("C","C","C","I","I","I")
ggbiplot(PCA, obs.scale = 1, var.scale = 1,
         groups = factor(c("C","C","C","I","I","I")), ellipse = TRUE, circle = TRUE,varname.size=0,var.axes = F) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')


#download E-MTAB-5980
atSamp5980 = getAE("E-MTAB-5980", type = "raw")


affy.data_5980 = ReadAffy(filenames = atSamp5980$rawFiles)
eset_5980 <- gcrma(affy.data_5980)
thisPD_5610<-pData(eset_5980)
thisExp_5980<-exprs(eset_5980)
boxplot(thisExp_5980)
boxplot(thisExp)
dim(thisExp_5980)
#log 2 normalize
thisExp_5980log<-log2(thisExp_5980)
dim(thisExp_5980log)
tDatalog<-t(thisExp_5980log)
dim(tDatalog)
#filter to remove 0 variance cols (genes)
tDatalogF <- tDatalog[,apply(tDatalog, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
PCA<-prcomp(tDatalogF, scale. = TRUE)
batch<-c("C","C","C","I","I","I")
ggbiplot(PCA, obs.scale = 1, var.scale = 1,
         groups = factor(c("C","C","C","I","I","I")), ellipse = TRUE, circle = TRUE,varname.size=0,var.axes = F) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')


#download E-MTAB-4782
atSamp4782 = getAE("E-MTAB-4782", type = "raw")


affy.data_4782 = ReadAffy(filenames = atSamp4782$rawFiles)
eset_4782 <- gcrma(affy.data_4782)
thisPD_4782<-pData(eset_4782)
thisExp_4782<-exprs(eset_4782)
boxplot(thisExp_4782)
boxplot(thisExp)
dim(thisExp_4782)
#log 2 normalize
thisExp_4782log<-log2(thisExp_4782)
dim(thisExp_4782log)
tDatalog<-t(thisExp_4782log)
dim(tDatalog)
#filter to remove 0 variance cols (genes)
tDatalogF <- tDatalog[,apply(tDatalog, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
PCA<-prcomp(tDatalogF, scale. = TRUE)
batch<-c("1","1","1","2","2","2","3","3","3","4","4","4","5","5","5","6","6","6","7","7","7","8","8","8")
batch<-c("1","1","1","1","1","1","2","2","2","2","2","2","5","5","5","5","5","5","5","5","5","5","5","5")
ggbiplot(PCA, obs.scale = 1, var.scale = 1,
         groups = factor(batch), ellipse = TRUE, circle = TRUE,varname.size=0,var.axes = F) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')



modcombat = model.matrix(~1, data=thisPD_4782)
combat_edata = ComBat(dat=thisExp_4782, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
cdataLogt<-t(log2(combat_edata))
tDatalogF <- cdataLogt[,apply(cdataLogt, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
PCA<-prcomp(tDatalogF, scale. = TRUE)


#combine datasets
affy.data_comb = ReadAffy(filenames = c(atSamp5980$rawFiles,atSamp5610$rawFiles,atSamp4782$rawFiles))
eset_comb <- gcrma(affy.data_comb,normalize=TRUE)
thisPD_comb<-pData(eset_comb)
thisExp_comb<-exprs(eset_comb)
boxplot(thisExp_comb)
boxplot(thisExp)
dim(thisExp_comb)
#log 2 normalize
thisExp_comblog<-log2(thisExp_comb)
dim(thisExp_comblog)
tDatalog<-t(thisExp_comblog)
dim(tDatalog)
#filter to remove 0 variance cols (genes)
tDatalogF <- tDatalog[,apply(tDatalog, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
PCA<-prcomp(tDatalogF, scale. = TRUE)
batch<-c("C","C","C","I","I","I","C","C","C","I","I","I")
batch<-c("1","1","1","1","1","1","2","2","2","2","2","2","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3")
batch<-c("1C","1C","1C","1B","1B","1B","2C","2C","2C","2B","2B","2B","3a","3a","3a","3b","3b","3b","3c","3c","3c","3d","3d","3d","3e","3e","3e","3f","3f","3f","3g","3g","3g","3h","3h","3h")


ggbiplot(PCA, obs.scale = 1, var.scale = 1,
         groups = factor(batch), ellipse = TRUE, circle = TRUE,varname.size=0,var.axes = F) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')








#batch correction on combined data
# parametric adjustment
rmData<-thisExp_comb[rowSums(thisExp_comb)>100,]
dim(rmData)
rmData<-thisExp_comb+1
combat_edata1 = ComBat(dat=thisExp_comb, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE,mean.only = F)
summary(combat_edata1)
boxplot(combat_edata1)
tcdata<-t(log2(combat_edata1+1))
dim(tcdata)
tDatalogF <- tcdata[,apply(tcdata, 2, var, na.rm=TRUE) != 0]
dim(tDatalogF)
PCA<-prcomp(tDatalogF, scale. = TRUE)

