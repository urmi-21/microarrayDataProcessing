library(dplyr)
library(purrr)
library(text2vec)
library(stringr)
library('lsa')
library(cluster)
#script to estimate batches from metadata

#Read all idf and sdrf files
files_idf= list.files(pattern = "\\.idf.txt$", full.names = TRUE,ignore.case = T)
files_sdrf= list.files(pattern = "\\.sdrf.txt$", full.names = TRUE,ignore.case = T)

#read all sdrf in data frame
sdrf<- files_sdrf %>% lapply(read_tsv) %>% bind_rows
#donot combine files
sdrf<- files_sdrf %>% lapply(read_tsv)


celFileColumn<-"Array Data File"






#text2vec ex
data("movie_review")
# select 500 rows for faster running times
movie_review = movie_review[1:500, ]
prep_fun = function(x) {
  x %>% 
    # make text lower case
    str_to_lower %>% 
    # remove non-alphanumeric symbols
    str_replace_all("[^[:alnum:]]", " ") %>% 
    # collapse multiple spaces
    str_replace_all("\\s+", " ")
}
movie_review$review_clean = prep_fun(movie_review$review)

doc_set_1 = movie_review[1:300, ]
it1 = itoken(doc_set_1$review_clean, progressbar = FALSE)
# specially take different number of docs in second set
doc_set_2 = movie_review[301:500, ]
it2 = itoken(doc_set_2$review_clean, progressbar = FALSE)


it = itoken(movie_review$review_clean, progressbar = FALSE)
v = create_vocabulary(it) %>% prune_vocabulary(doc_proportion_max = 0.1, term_count_min = 5)
vectorizer = vocab_vectorizer(v)


# they will be in the same space because we use same vectorizer
# hash_vectorizer will also work fine
dtm1 = create_dtm(it1, vectorizer)
dim(dtm1)

dtm2 = create_dtm(it2, vectorizer)
dim(dtm2)

d1_d2_jac_sim = sim2(dtm1,method = "jaccard", norm = "none")


d1_d2_cos_sim = sim2(dtm1, dtm2, method = "cosine", norm = "l2")

d1_d2_cos_sim[1:2, 1:5]











#########compute similarity of rows################


f1<-sdrf[[2]]

#1 remove columns with no information
f1<-Filter(function(x)(length(unique(x))>1), f1)
#compute similarity of each row
td = tempfile()
dir.create(td)
for(i in 1:nrow(f1)){
  write( as.character(f1[i,]), file=paste(td, paste("R",i,sep=""), sep="/"))
}
myMatrix = textmatrix(td, minWordLength=2)
sim<-cosine(myMatrix)
boxplot(sim)



HC=hclust(as.dist(1-sim), method="single")
cutree(HC,k=3)

plot(HC, cex = 0.6)
rect.hclust(HC, k = 4, border = 2:5)

kmeanclust <- kmeans(sim,3)

kmeanclust$cluster

fviz_cluster(kmeanclust, data = sim, geom = "point",
             stand = FALSE, frame.type = "norm")
#elbow method
fviz_nbclust(sim, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)

#gap stat

gap_stat <- clusGap(sim, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method="Tibs2001SEmax")

fviz_gap_stat(gap_stat)




library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(NbClust)
