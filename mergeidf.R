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
    ttfile<-ttfile[-1,]
    listOfDF[[i]]<-ttfile
    
  }
  
  return(listOfDF)
}


#megre sdrf files downloaded from arrayexpress
df_idf<-list.files(full.names = TRUE,path = "idfDir/", pattern = "*.txt") %>% lapply(read_tsv,col_names = F)%>% makeasChar %>% bind_rows()

write_tsv(df_idf,"idf_summary.tsv")

#read and join idf and sdrf
combined_srdf_Summary <- read_delim("~/Downloads/Immuno-Navigator_datasets/human/combined_srdf_Summary.tsv", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
idf_summary <- read_delim("~/Downloads/Immuno-Navigator_datasets/human/idf_summary.tsv", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

combinedMD<- inner_join(idf_summary,combined_srdf_Summary)

write_tsv(combinedMD,"combined_metadata.tsv")
