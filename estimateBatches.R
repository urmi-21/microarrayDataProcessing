library(dplyr)
library(purrr)
library(text2vec)
library(stringr)
#script to estimate batches from metadata

#Read all idf and sdrf files
files_idf= list.files(pattern = "\\.idf.txt$", full.names = TRUE,ignore.case = T)
files_sdrf= list.files(pattern = "\\.sdrf.txt$", full.names = TRUE,ignore.case = T)

#read all sdrf in data frame
sdrf<- files_sdrf %>% lapply(read_tsv) %>% bind_rows

celFileColumn<-"Array Data File"






