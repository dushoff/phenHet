## following https://cran.r-project.org/web/packages/easyPubMed/vignettes/getting_started_with_easyPubMed.html

library(easyPubMed)
## example from notes.Rmd: should retrieve 276 records?
my_query <- "epidem* AND infect* AND (heterogeneity OR structure OR network) AND dynamics AND nonlinear"
my_entrez_id <- get_pubmed_ids(my_query)

as.numeric(my_entrez_id$Count)

my_abstracts_xml <- fetch_pubmed_data(pubmed_id_list = my_entrez_id)
## my_df <- article_to_df(my_abstracts_xml)
my_PM_list <- articles_to_list(pubmed_data = my_abstracts_xml)

length(my_PM_list)

## etc.
