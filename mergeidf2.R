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
#df_idf<-list.files(full.names = TRUE,path = "idf/", pattern = "*.txt") %>% lapply(read_tsv,col_names = F)%>% makeasChar %>% bind_rows()



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
write_tsv(df_idf,"idf_summary.tsv")

##read sdrf
#read all files
# sdrf from E-MTAB-1264 has duplicated rows for the CEL files it was edited and duplicates were removed
flistsdrf<-list.files(full.names = TRUE,path = "sdrf/", pattern = "*.txt")

listFiles<-list()
k<-1
for(f in flistsdrf){
  
  print(f)
  #read f and store in list
  thisF<-read_tsv(f,col_names = T)
  #thisFname<-tools::file_path_sans_ext(tools::file_path_sans_ext(basename(f)))
  thisFname<-unlist(strsplit(basename(f), "." ,fixed = T))[1]
  #remove duplicate columns
  colnames(thisF)<-make.unique(colnames(thisF))
  #remove .idf
  thisF <- thisF %>% mutate(Accession=thisFname)
  thisF <- thisF %>% mutate_all(as.character)
  listFiles[[k]]<-thisF
  k<-k+1
}

#bind rows
df_sdrf<-bind_rows(listFiles)
#rename col
colnames(df_sdrf)[which(colnames(df_sdrf)=="Array Data File")] <- "Array_Data_File"
df_sdrf[is.na(df_sdrf)]<-"NA"
#check all rows have Array Data File
df_sdrf$`Source Name`[which(df_sdrf$Array_Data_File=="NA")]

#extract useful info
keyWords<-c("organism","cell","disease","treatment","arrayexpress","Array_Data_File","Material Type","Protocol","Sample","Description","Characteristics")

colsToKeep<-grepl(paste(keyWords,collapse = "|"),colnames(df_sdrf),ignore.case = T)

colsToKeep[which(colnames(df_sdrf)=="Accession")]=TRUE

df_sdrf_final<-df_sdrf[,colsToKeep]

#order by col name
df_sdrf_final<-df_sdrf_final[,colnames(df_sdrf_final)[order(colnames(df_sdrf_final))]]
#remove all cols with NA
df_sdrf_final[df_sdrf_final=="NA"]<-NA
df_sdrf_final[df_sdrf_final==""]<-NA
df_sdrf_final2 <- df_sdrf_final[, colSums(is.na(df_sdrf_final)) < nrow(df_sdrf_final) ]

#merge columns by keywords
#init new df
test<-data.frame()
test<-as.data.frame(df_sdrf_final[,"Accession"])


df_sdrf_final[is.na(df_sdrf_final)]<-""

for(keyw in keyWords){
  print(keyw)
  thisColName<-paste(keyw,"Info",sep = "_")
  thisCols<-colnames(df_sdrf_final)[grepl(keyw,colnames(df_sdrf_final),ignore.case = T)]
  tempdf<-as.data.frame(df_sdrf_final[ , thisCols ])
  test[,thisColName] <- apply( tempdf , 1 , function(row) paste(row[nzchar(row)], collapse = ";") )
}

#remove first column
#test<-test[,c(2:ncol(test))]

df_sdrf_final <- test %>% distinct

df_sdrf_final<- df_sdrf_final %>% mutate(DataColum = paste(Accession,Array_Data_File_Info,sep = "_"))

#duplicates
df_sdrf_final$DataColum[duplicated(df_sdrf_final$DataColum)]

#write to file
write_tsv(df_sdrf_final,"sdrf_summary.tsv")





combinedMD<- left_join(df_idf,df_sdrf_final)
combinedMD <- combinedMD %>% distinct

#datacolumns in data
Acc_to_Cel <- read_delim("Acc_to_Cel.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

combinedMD<-combinedMD %>% filter(DataColum %in% Acc_to_Cel$Combined)

combinedMD <- data.frame(lapply(combinedMD, as.character), stringsAsFactors=FALSE)

combinedMD[is.na(combinedMD)]<-"NA"
combinedMD[combinedMD==""]<-"NA"
combinedMD2 <- combinedMD %>% distinct




repeatedDC <- as.data.frame(combinedMD$DataColum[duplicated(combinedMD$DataColum)])


write_tsv(repeatedDC,"repeatedDC.tsv")

write_tsv(combinedMD,"combined_metadata.tsv")




