# install.packages("devtools")
# library(devtools)
## install_github("ropensci/rentrez")
## BMB: OK to use CRAN version instead?

## Richard: Tried several times on Dec 11th, but not working at the moment
## (It said rentrez is not available for current version of R (V4.3.3, V4.4.1,
## V4.4.2.) 
## Tried on Dec 15th on the same device and worked. Seems like an update based
## on V4.4.2 has just been updated to CRAN.
## install.packages("rentrez")

library("rentrez")
library("stringr")
library("easyPubMed")
## https://support.nlm.nih.gov/kbArticle/?pn=KA-05317
## https://account.ncbi.nlm.nih.gov/
if (!nzchar(ncbi_api_key <- Sys.getenv("NCBI_API_KEY"))) {
    stop("please get an API key from https://account.ncbi.nlm.nih.gov/ ",
         "(see https://support.nlm.nih.gov/kbArticle/?pn=KA-05317) ",
         "and set the environmental variable NCBI_API_KEY to it")
}

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

## https://academia.stackexchange.com/questions/191088/how-can-i-get-around-the-10000-search-result-limit-in-pubmed
Q1_result <- entrez_search(db="pubmed", term=query_Q1, retmax=10000)
length(Q1_result$ids)
Q1_ids<-Q1_result$ids

## https://academia.stackexchange.com/questions/191088/how-can-i-get-around-the-10000-search-result-limit-in-pubmed

## need to set the 'retstart' parameter, loop over batches ...

## Richard: I get the idea of 'retstart' and I'll give it a try if possible(not 
## sure if I can still connect to NCBI from China use VPNs, I'll see when I arri
## ve). I think it is reasonable to just get PMIDs for queries with >10k results
## using entrez_search(with a loop machine), since for these queries we care more 
## about if they overlap with our target papers.

## https://www.nlm.nih.gov/dataguide/eutilities/utilities.html
## https://github.com/ropensci/rentrez/issues/180
if (FALSE) {
    
    system.time(
        Q2_results_EPM <- easyPubMed::batch_pubmed_download(
                                          pubmed_query_string = query_Q2,
                                          api_key = ncbi_api_key,
                                          batch_size = 1000,
                                          ## restart counter
                                          res_cn = 1, 
                                          dest_file_prefix = "Test1_Q2", 
                                          encoding = "ASCII")
    )
    ## 72 batches, ~25M each, 1.8G total
    ## started at approx 5000/minute, but seems throttled after 10 batches?

    ## would like to be able to get *just* PMIDs (rather than fetching PMIDs + metadata), ... ??

}

## we can get just PMIDs by looping over this:
xx <- readLines("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=cancer&reldate=60&datetype=edat&retmax=100&retstart=101")

## we can definitely make this into a machine, either by brute force (grep "<Id>..." and extract) or by using an XML parser

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

## BMB: extract info to data frame (suitable for write.csv)

rentrez_df <- function(x) {
    csv_df <- t(extract_from_esummary(summs, c("uid", "pubdate", "title", "fulljournalname")))
    collapse_auths <- function(x) paste(x$authors$name[x$authors$authtype == "Author"], collapse = ";")
    auths <- (extract_from_esummary(summs, "authors", simplify = FALSE)
        |> sapply(collapse_auths)
    )
    csv_df <- data.frame(csv_df, auths)
    csv_df$pubdate <- stringr::str_extract(csv_df$pubdate, "[0-9]+")
    return(csv_df)
}

Q6_df <- rentrez_df(Q6_result)

## Richard: Got it! Thanks a lot.