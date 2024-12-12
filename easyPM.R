## following https://cran.r-project.org/web/packages/easyPubMed/vignettes/getting_started_with_easyPubMed.html

library(easyPubMed)
## example from notes.Rmd: should retrieve 276 records?
query_Q6 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynam* AND \"nonlinear incidence\""
query_Q5 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynam* AND nonlinear AND incidence"
query_Q4 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynamics AND nonlinear"
query_Q3 <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynamics"
query_Q2 <- "epidem* AND infect* AND (heterogeneity OR structure OR network)"
query_Q1 <- "epidem* AND (heterogeneity OR structure OR network)"

my_entrez_id <- get_pubmed_ids(query_Q2)
as.numeric(my_entrez_id$Count)

# batch_pubmed_download(query_Q2, dest_dir = NULL,
#                       dest_file_prefix = "Q2_data_",
#                       format = "xml", api_key = NULL,
#                       batch_size = 2000, res_cn = 1,
#                       encoding = "UTF8")

id_Q2 <- fetch_all_pubmed_ids(my_entrez_id)
length(id_Q2)

# Both easyPM and Rentrenz cannot download record or fetch ids more than 10,000
# from pubmed/NCBI.
# 

# my_abstracts_xml <- fetch_pubmed_data(pubmed_id_list = my_entrez_id)
# ## my_df <- article_to_df(my_abstracts_xml)
# my_PM_list <- articles_to_list(pubmed_data = my_abstracts_xml)
# 
# length(my_PM_list)

## etc.

