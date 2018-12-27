BiocManager::install("ArrayExpress", version = "3.8")
BiocManager::install("affy", version = "3.8")
library("ArrayExpress")
library(affy)
library(affydata)
library(ggbiplot)
memory.limit(size=56000)

#code Ref: http://math.usu.edu/jrstevens/stat5570/2.4.Visualization.pdf
data(Dilution)
eset<-log2(exprs(Dilution))
X <-eset[, 1]
Y <-eset[, 2]
M <- Y-X
A <- (Y + X) / 2
plot(A,M,main= 'default M-A plot',pch= 16,cex = 0.1);abline(h = 0)

pc <- prcomp(t(eset))

names(pc)
group <- as.numeric(colnames(eset))
plot(pc$x[, 1], pc$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2",col=c(1,2,3,4))
g <- ggbiplot(pc, scale = 1,obs.scale = 1, varname.abbrev = FALSE,var.axes = FALSE,pc.biplot =TRUE,circle = TRUE)
g
##using NOISeq
myPCA<-PCA.GENES(t(eset))
myPCA$var.exp
explo.plot(myPCA)
explo.plot(myPCA)

myPCA = dat(eset, type = "PCA")
