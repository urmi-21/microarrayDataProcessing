#checkPlatform in CEL files
library(affyio)
files = list.files(pattern = "\\.CEL$", full.names = TRUE,ignore.case=T)
platforms<- table(sapply(files, function(x) read.celfile.header(x)$cdfName))
filePlatform<-sapply(files, function(x) read.celfile.header(x)$cdfName)
ATH1files<-names(filePlatform[which(filePlatform=="ATH1-121501")])
#affy.data = ReadAffy(filenames = ATH1files)

corrupted.files <- c()
goodFiles<-c()
for(i in 1:length(ATH1files)) {
    err <<- F
    x <- tryCatch(ReadAffy(filenames=ATH1files[i]), error = function(e){print(e); err <<- T}, finally=print(paste("finished with", ATH1files[i], "at", Sys.time())))
    if(err){
        corrupted.files <- c(corrupted.files,ATH1files[i])
        print(paste("curroptedFile", ATH1files[i]))
    }
    else{
        goodFiles<-c(goodFiles,ATH1files[i])
        }
}

affy.data = ReadAffy(filenames = goodFiles)
