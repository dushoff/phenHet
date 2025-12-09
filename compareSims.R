### Package Part
library(igraph)
library(Rcpp)

library(shellpipes)
source("scaleFuns.R")
sourceFiles()

pars0 <- list(beta = 1, gamma = 1, r = 1, lambda = 5)
Rcpp::sourceCpp("edgelist.cpp")
simEdge <- GilAlgoCpp
Rcpp::sourceCpp("NetSimulator.cpp")
simVertex <- GilAlgoCpp

runSim <- function(N = 1e4, params, simtype = c("vertex", "edge"), seed = NULL, MaxTime = 100) {
  attach(params)
  if (!is.null(seed)) set.seed(seed)

  start_time <- proc.time()
  nseq <- rnbinom(N,r,mu=lambda)
  while(!CheckSeq(nseq)){
    nseq <- rnbinom(N,r,mu=lambda)
  }
  ## generating graph
  G <- sample_degseq(  nseq
                   , method = "fast.heur.simple"
                     )
  
  stopifnot(!any(sort(degree(G))-sort(nseq)!=0))

  ## Translate igraph network object into adjacency matrix
  Adj_list <- as_adj_list(  G
                        , mode = "all"
                        , loops = "once"
                        , multiple = TRUE
                          )
  network_setup_time <- proc.time() -  start_time


  simFun <- switch(simtype,
                   vertex = simVertex,
                   edge = simEdge)
  
  sim_time <- system.time(sim_result <- simFun(Adj_list, N, beta, gamma, MaxTime = MaxTime))
  return(tibble::lst(param, N, simtype, network_setup_time, sim_time, sim_result))
}
  
