library(readr)
library(magrittr)
library(dplyr)
library(purrr)
library(stringr)
makeasChar<-function(listOfDF){
  t<-length(listOfDF)
  
  for(i in 1:t){
    listOfDF[[i]]<- as.data.frame(listOfDF[[i]] %>% mutate_all(as.character), stringsAsFactors=FALSE)
    tfile<- listOfDF[[i]]
    tfile[is.na(tfile)]<-"NA"
    tfile<-tfile[,1:2]
    ttfile<-data.frame(t(tfile))
    colnames(ttfile)<- as.character(unlist(ttfile[1,]))
    
    
   
    
    listOfDF[[i]]<-ttfile
    
  }
  
  return(listOfDF)
}


#megre sdrf files downloaded from arrayexpress
df_idf<-list.files(full.names = TRUE,path = "idfDir/", pattern = "*.txt") %>% lapply(read_tsv,col_names = F)%>% makeasChar %>% bind_rows()

