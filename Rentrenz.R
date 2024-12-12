# install.packages("devtools")
# library(devtools)
# install_github("ropensci/rentrez")
library("rentrez")

query_Q6 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynam* AND \"nonlinear incidence\""
query_Q5 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynam* AND nonlinear AND incidence"
query_Q4 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynamics AND nonlinear"
query_Q3 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynamics"
query_Q2 <- "epidem* AND infect* AND (heterogeneity OR structure OR network)"
query_Q1 <- "epidem* AND (heterogeneity OR structure OR network)"

query_Q7 <- "(epidem* OR disease) AND (hetero* OR structure OR network OR )"

query_Author <- "AND ((Gomes, Mgm[Author]) OR (Dwyer, G[Author]) OR (Dushoff, J[Author]) OR (Elkinton, Js[Author]) OR (Keeling, Mj[Author]) OR (Grenfell, Bt[Author]) OR (Granich, Rm[Author]) OR (Wilson, EB[Author]) OR (Worcester, J[Author]) OR (Hethcote, HW[Author]) OR (Levin, SA[Author]) OR (Liu, W[Author]))"

query_Q1A <- paste(query_Q1, query_Author, sep=" ")
query_Q2A <- paste(query_Q2, query_Author, sep=" ")

TargetList <- c("1a"="32511451","2a"="18811331","2b"="10856195","3b"="10343409","4a"="19038438","5a"="16588678","5b"="3668394","5c"="3958634")
names(TargetList)

Q1_result <- entrez_search(db="pubmed", term=query_Q1, retmax=10000)
length(Q1_result$ids)
Q1_ids<-Q1_result$ids

Q1A_result <- entrez_search(db="pubmed", term=query_Q1A, retmax=10000)
length(Q1A_result$ids)
Q1A_ids<-Q1A_result$ids

Q2_result <- entrez_search(db="pubmed", term=query_Q2, retmax=10000)
length(Q2_result$ids)
Q2_ids<-Q2_result$ids

Q2A_result <- entrez_search(db="pubmed", term=query_Q2A, retmax=10000)
length(Q2A_result$ids)
Q2A_ids<-Q2A_result$ids

Q3_result <- entrez_search(db="pubmed", term=query_Q3, retmax=10000)
length(Q3_result$ids)
Q3_ids<-Q3_result$ids

Q4_result <- entrez_search(db="pubmed", term=query_Q4, retmax=10000)
length(Q4_result$ids)
Q4_ids<-Q4_result$ids

Q5_result <- entrez_search(db="pubmed", term=query_Q5, retmax=10000)
length(Q5_result$ids)
Q5_ids<-Q5_result$ids

Q6_result <- entrez_search(db="pubmed", term=query_Q6, retmax=10000)
length(Q6_result$ids)
Q6_ids<-Q6_result$ids


TargetList[TargetList %in% Q1_ids]
TargetList[TargetList %in% Q2_ids]
TargetList[TargetList %in% Q3_ids]
TargetList[TargetList %in% Q4_ids]
TargetList[TargetList %in% Q5_ids]
TargetList[TargetList %in% Q6_ids]

TargetList[TargetList %in% Q1A_ids]
TargetList[TargetList %in% Q2A_ids]


query_Q7 <- "(epidem* OR disease) AND (hetero* OR structure OR network) "
query_Q7A <- paste(query_Q7, query_Author, sep=" ")

Q7_result <- entrez_search(db="pubmed", term=query_Q7, retmax=10000)
length(Q7_result$ids)
Q7_result$count

Q7A_result <- entrez_search(db="pubmed", term=query_Q7A, retmax=10000)
length(Q7A_result$ids)
Q7A_ids<-Q7A_result$ids
TargetList[TargetList %in% Q7A_ids]

query_Q8 <- "(epidem* OR disease) AND (hetero* OR structure OR network) AND dynamic"
query_Q8A <- paste(query_Q8, query_Author, sep=" ")

Q8_result <- entrez_search(db="pubmed", term=query_Q8, retmax=10000)
length(Q8_result$ids)
Q8_result$count
# Q8_ids<-Q8_result$ids
# TargetList[TargetList %in% Q8_ids]

Q8A_result <- entrez_search(db="pubmed", term=query_Q8A, retmax=10000)
length(Q8A_result$ids)
Q8A_ids<-Q8A_result$ids
TargetList[TargetList %in% Q8A_ids]

query_Q9 <- "(epidem* OR disease) AND (hetero* OR structure OR network) AND nonlinear"
query_Q9A <- paste(query_Q9, query_Author, sep=" ")

Q9_result <- entrez_search(db="pubmed", term=query_Q9, retmax=10000)
length(Q9_result$ids)
Q9_result$count
Q9_ids<-Q9_result$ids
TargetList[TargetList %in% Q9_ids]

