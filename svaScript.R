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


celData<-celData %>% mutate(batchId=as.numeric(factor(celData$Exp)))

#read all cel files downloaded
files = list.files(pattern = "\\.CEL$", full.names = TRUE)
affy.data = ReadAffy(filenames = files)

#apply gcrma
eset <- gcrma(affy.data)
thisPD<-pData(eset)
thisExp<-exprs(eset)
dim(thisExp)

mod0 = model.matrix(~1,data=celData)

mod = model.matrix(~as.factor(Exp), data=celData)

n.sv = num.sv(thisExp,mod0,method="leek")

svobj = sva(thisExp,mod,mod0,n.sv=n.sv)

pValues = f.pvalue(thisExp,mod,mod0)
qValues = p.adjust(pValues,method="BH")

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)

pValuesSv = f.pvalue(thisExp,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

svaObj0=sva(thisExp,mod)
cleany = cleaningY(thisExp,mod0,svaObj0)


#do pca
#log 2 normalize
thisExplog<-log2(thisExp)
dim(thisExplog)
tData<-t(thisExp) #changes dimention ???
dim(tData)
tDatalog<-t(thisExplog)
dim(tDatalog)
#filter to remove 0 variance cols (genes)
tDatalogF <- tDatalog[,apply(tDatalog, 2, var, na.rm=TRUE) != 0]
esetPCA<-prcomp(tDatalogF, scale. = TRUE)
ggbiplot(esetPCA, obs.scale = 1, var.scale = 1,
         groups = factor(celData$Exp), ellipse = TRUE, circle = TRUE,varname.size=0,var.axes = F) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')


cleanylog<-log2(cleany+10)
dim(cleanylog)
tcleany<-t(cleany) #changes dimention ???
dim(cleany)
tcleanylog<-t(cleanylog)
dim(tcleanylog)
#filter to remove 0 variance cols (genes)
tDatalogF <- tcleanylog[,apply(tcleanylog, 2, var, na.rm=TRUE) != 0]
esetPCA<-prcomp(tDatalogF, scale. = TRUE)
ggbiplot(esetPCA, obs.scale = 1, var.scale = 1,
         groups = factor(celData$Exp), ellipse = TRUE, circle = TRUE,varname.size=0,var.axes = F) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')








#ref: https://support.bioconductor.org/p/47350/
cleaningY = function(y, mod, svaobj) {
  X = cbind(mod, svaobj$sv) 
  Hat = solve(t(X) %*% X) %*% t(X) 
  beta = (Hat %*% t(y)) 
  P = ncol(mod) 
  cleany = y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P), ]) 
  return(cleany)
}



