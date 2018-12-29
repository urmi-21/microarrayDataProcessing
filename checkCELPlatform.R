#checkPlatform in CEL files
library(affyio)
files = list.files(pattern = "\\.CEL$", full.names = TRUE,ignore.case=T)
platforms<- table(sapply(files, function(x) read.celfile.header(x)$cdfName))
filePlatform<-sapply(files, function(x) read.celfile.header(x)$cdfName)
ATH1files<-names(filePlatform[which(filePlatform=="ATH1-121501")])
