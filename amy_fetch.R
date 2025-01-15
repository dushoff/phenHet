if (packageVersion("easyPubMed") < "3.0.0") {
    stop('please install a later version of easyPubMed ',
         '(e.g. remotes::install_github("dami82/easyPubMed") )')
}
library(easyPubMed)
library(tidyverse)
do_slow <- FALSE

if (!nzchar(ncbi_api_key <- Sys.getenv("NCBI_API_KEY"))) {
    warning("you may want to get an API key from https://account.ncbi.nlm.nih.gov/ ",
         "(see https://support.nlm.nih.gov/kbArticle/?pn=KA-05317) ",
         "and set the environmental variable NCBI_API_KEY to it")
    ncbi_api_key <- NULL
}



query_Q6 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynam* AND \"nonlinear incidence\""
query_Q5 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynam* AND nonlinear AND incidence"
query_Q4 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynamics AND nonlinear"
query_Q3 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynamics"
query_Q2 <- "epidem* AND infect* AND (heterogeneity OR structure OR network)"
query_Q1 <- "epidem* AND (heterogeneity OR structure OR network)"

epm <- epm_query(query_Q4)
str(epm@uilist)  ## empty list, no fetch yet
epm_f <- epm_fetch(epm
                 , format = 'medline'
                 , write_to_file = TRUE
                 , outfile_prefix = "testPubMed_Q4"
                 , store_contents = TRUE)
unlist(epm_f@uilist[1:5])

## uilist ONLY
epm_f2 <- epm_fetch(epm, format = 'uilist')

epm_Q1 <- epm_query(query_Q1)
if (do_slow) {
    epm_Q1_f2 <- epm_fetch(epm_Q1, format = 'uilist',
                           api_key = ncbi_api_key)
} else warning("skipping slow step!")

