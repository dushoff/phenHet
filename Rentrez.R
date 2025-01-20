library("rentrez")
library("stringr")
library("easyPubMed")
## https://support.nlm.nih.gov/kbArticle/?pn=KA-05317
## https://account.ncbi.nlm.nih.gov/
## please DO NOT PUT API KEYS into public code ...

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

query_Author <- "AND ((Gomes, Mgm[Author]) OR (Dwyer, G[Author]) OR (Dushoff, J[Author]) OR (Elkinton, Js[Author]) OR (Keeling, Mj[Author]) OR (Grenfell, Bt[Author]) OR (Granich, Rm[Author]) OR (Wilson, EB[Author]) OR (Worcester, J[Author]) OR (Hethcote, HW[Author]) OR (Levin, SA[Author]) OR (Liu, W[Author]) OR (Berestycki, H[Author]) OR (Rose, C[Author]) OR (Korobeinikov, A[Author]) OR (Novozhilov, AS[Author]))"

query_Q1A <- paste(query_Q1, query_Author, sep=" ")
query_Q2A <- paste(query_Q2, query_Author, sep=" ")

TargetList <- c("1a"="35189135","2a"="18811331","2b"="10856195","3a"="10343409","4a"="19038438","5a"="16588678","5b"="3668394","5c"="3958634","6a"="34314731","7a"="36964799","8a"="16794947","8b"="17443392","9a"="18722386")
names(TargetList)

## https://academia.stackexchange.com/questions/191088/how-can-i-get-around-the-10000-search-result-limit-in-pubmed
Q1_result <- entrez_search(db="pubmed", term=query_Q1, retmax=5000)
length(Q1_result$ids)
Q1_result$count
Q1_ids<-Q1_result$ids


## https://academia.stackexchange.com/questions/191088/how-can-i-get-around-the-10000-search-result-limit-in-pubmed

## need to set the 'retstart' parameter, loop over batches ...

## https://www.nlm.nih.gov/dataguide/eutilities/utilities.html
## https://github.com/ropensci/rentrez/issues/180

## Richard: I get the idea of 'retstart' and I'll give it a try if possible(not 
## sure if I can still connect to NCBI from China use VPNs, I'll see when I arri
## ve). I think it is reasonable to just get PMIDs for queries with >10k results
## using entrez_search(with a loop machine), since for these queries we care more 
## about if they overlap with our target papers.

# entrez_search(db="pubmed", term=query_Q1, retmax=5000,retstart=10001)
## Not working with retstart>10k: Error in ans[[1]] : subscript out of bounds 
## Loop with smaller retstart is not working as well

## Try loop over batches, batchsize=5000
# batchsize <- 5000
# loop_n_Q1 <- ceiling(Q1_result$count/batchsize)
# 
# Q1_ids <- c()
# for (i in c(1:loop_n_Q1)){
#   temp_result <- entrez_search(db="pubmed"
#                                , term=query_Q1
#                                , retmax=batchsize
#                                , retstart=batchsize*(i-1))
#   #print(temp_id$count)
#   temp_ids <- temp_result$ids
#   Q1_ids <- append(Q1_ids,temp_ids)
# }
# length(Q1_ids)


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
## https://www.nlm.nih.gov/dataguide/eutilities/utilities.html#esearch
my_fetch0 <- function(query, retmax, retstart, rettype = "uilist",
                      max_retry = 5,
                      retry_pause = 5) {
    qq <- gsub(" ", "+", query)
    base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed"
    url <- sprintf("%s&term=%s&retmax=%d&retstart=%d&rettype=%s",
                   base_url, qq, retmax, retstart, rettype)
    tries <- 0
    repeat {
        out <- readLines(url)
        msg <- grep("<ERROR>", out, value = TRUE)
        if (length(msg) == 0 || tries >= max_retry) break
        cat(msg, "\n")
        cat(sprintf("retstart = %d, retmax = %d\n", retstart, retmax))
        cat(sprintf("trying again, pausing for %d seconds ...\n", retry_pause))
        Sys.sleep(retry_pause)
    }
    ## FIXME: return character(0) + warning if we fail?
    if (tries > max_retry) {
        stop("query failure: ", msg)
    }
    tag <- switch(rettype,
                  count = "</?Count>",
                  uilist = "</?Id>")
    vals <- (grep(tag, out, value = TRUE)
        |> gsub(pattern = sprintf("^.*%s(.*)%s.*$", tag, tag), replacement = "\\1")
    )
    if (rettype == "count") vals <- as.integer(vals)
    ## could convert PMIDs to numeric but ... ?
    return(vals)
}

my_fetch_all <- function(query, batchsize = 1000, verbose = FALSE, pause = 5) {
    ## retmax doesn't seem to matter for 'count'
    cc <- my_fetch0(query, 100, 1, rettype = "count")
    retstart_vec <- seq(1, cc, by = batchsize)
    if (verbose) pb <- txtProgressBar(max = cc, style = 3)
    ## FIXME: break if we start to fail? return partial answers?
    ids <- lapply(retstart_vec,
                  function(rs) {
                      if (verbose) setTxtProgressBar(pb, rs)
                      my_fetch0(retstart = rs,
                                retmax = batchsize,
                                query = query)
                      Sys.sleep(pause)
                  })
    browser()
}

# x <- my_fetch0("cancer", 100, 1)
# x <- my_fetch0(query, 100, 1, "count")
# x <- my_fetch0(query_Q1, 100, 1)
# x <- my_fetch0(query_Q1, 100, 1, "count")

## still fails once we get to retstart = 10001??
if (FALSE) {
    my_fetch_all(query_Q1, verbose = TRUE)
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


query_Q7 <- "(epidem* OR disease OR infect*) AND (dynam* OR equilibri*) AND (inciden* OR susceptib* OR spread OR transmission) AND (hetero* OR network OR variation OR nonlinear OR non-linear)"
query_Q7A <- paste(query_Q7, query_Author, sep=" ")

Q7_result <- entrez_search(db="pubmed", term=query_Q7, retmax=10000)
length(Q7_result$ids)
Q7_result$count

Q7A_result <- entrez_search(db="pubmed", term=query_Q7A, retmax=10000)
length(Q7A_result$ids)
Q7A_ids<-Q7A_result$ids
TargetList[TargetList %in% Q7A_ids]

query_Q8 <- "(epidem* OR disease OR infect*) AND (dynam* OR equilibri*) AND (inciden* OR susceptib* OR spread OR transmission) AND (hetero* OR variation OR nonlinear OR non-linear)"
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

query_Q9 <- "(epidem* OR disease OR infect*) AND (dynam* OR equilibri*) AND (inciden* OR susceptib* OR spread OR transmission) AND (hetero* OR (individual variation) OR nonlinear OR non-linear)"
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
