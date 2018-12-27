BiocManager::install("ArrayExpress", version = "3.8")
BiocManager::install("affy", version = "3.8")
library("ArrayExpress")
library(affy)
library(affydata)

#code Ref: http://math.usu.edu/jrstevens/stat5570/2.4.Visualization.pdf
data(Dilution)
eset<-log2(exprs(Dilution))
X <-eset[, 1]
Y <-eset[, 2]
M <- Y-X
A <- (Y + X) / 2
plot(A,M,main= 'default M-A plot',pch= 16,cex = 0.1);abline(h = 0)
