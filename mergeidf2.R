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

cleanSdrf<-function(df){
  
  df<-as.data.frame(df %>% mutate_all(as.character), stringsAsFactors=FALSE)
  tfile<- df
  tfile[is.na(tfile)]<-"NA"
  tfile<-tfile[,1:2]
  ttfile<-data.frame(t(tfile))
  colnames(ttfile)<- as.character(unlist(ttfile[1,]))
  ttfile<-ttfile[-1,]
  return(ttfile)
}

#megre sdrf files downloaded from arrayexpress
df_idf<-list.files(full.names = TRUE,path = "idf/", pattern = "*.txt") %>% lapply(read_tsv,col_names = F)%>% makeasChar %>% bind_rows()

write_tsv(df_idf,"idf_summary.tsv")

#add accession

#read all files
flist<-list.files(full.names = TRUE,path = "idf/", pattern = "*.txt")

listFiles<-list()
k<-1
for(f in flist){
  
  print(f)
  #read f and store in list
  thisF<-read_tsv(f,col_names = F)
  #clean this F (transpose)
  thisF<-cleanSdrf(thisF)
  thisFname<-tools::file_path_sans_ext(tools::file_path_sans_ext(basename(f)))
  #remove duplicate columns
  colnames(thisF)<-make.unique(colnames(thisF))
  #remove .idf
  thisF <- thisF %>% mutate(Accession=thisFname)
  listFiles[[k]]<-thisF
  k<-k+1
  
}

#bind rows
df_idf<-bind_rows(listFiles)




#read and join idf and sdrf
combined_srdf_Summary <- read_delim("sdrf_summary.tsv", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
idf_summary <- read_delim("idf_summary.tsv", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)


colnames(combined_srdf_Summary)[ncol(combined_srdf_Summary)] <- "Cel"

combined_srdf_Summary<-left_join(Acc_to_Cel,combined_srdf_Summary)

#column 7 has accessions
colnames(idf_summary)[7]<-"Exp"

combinedMD<- inner_join(idf_summary,combined_srdf_Summary)

#remove NA cols
combinedMD2 <- combinedMD[, colSums(is.na(combinedMD)) < nrow(combinedMD) ]

#remove repeated rows

combinedMD3<- data.frame(lapply(combinedMD2,as.character),stringsAsFactors = F)

combinedMD3[is.na(combinedMD3)] <- "NA"

combinedMD3 <- combinedMD3 %>% distinct()



#duplicated runids
combinedMD3$Combined[duplicated(combinedMD3$Combined)]
















write_tsv(combinedMD,"combined_metadata.tsv")
