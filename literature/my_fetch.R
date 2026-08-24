query_Q1 <- "epidem* AND (heterogeneity OR structure OR network)"
query <- "cancer"

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

print(x <- my_fetch0("cancer", 100, 1))
print(x <- my_fetch0(query, 100, 1, "count"))
print(x <- my_fetch0(query_Q1, 100, 1))
print(x <- my_fetch0(query_Q1, 100, 1, "count"))

## still fails once we get to retstart = 10001??
if (TRUE) {
    my_fetch_all(query_Q1, verbose = TRUE)
}
