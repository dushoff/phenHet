## short/self-contained tests for GIllespie algorithm code development
library(Rcpp)

source("NetworkSimulator.R")

lambda <- 5
kappa <- 2
r <- 1/kappa

(p<-1/(1+kappa*lambda))
(v<-lambda/p)

#### Disease Parameter
beta <- 0.25
# gamma <- 0.75
gamma <- 0.2

N <- 250000
set.seed(2853)
seq <- rnbinom(N,r,mu=lambda)
while(!CheckSeq(seq)){
  seq <- rnbinom(N,r,mu=lambda)
}
## generating graph
G <- sample_degseq(  seq
                   , method = "fast.heur.simple"
                   #, method = "configuration.simple"
                   )


# check realization is successful
stopifnot(!any(sort(degree(G))-sort(seq)!=0))

# Translate igraph network object into adjacency matrix
Adj_list <- as_adj_list(  G
                        , mode = "all"
                        , loops = "once"
                        , multiple = TRUE
)

### Rcpp Version
sourceCpp('Full_Copilot.cpp')
args1 <- list(N, beta, gamma, MaxTime = 2.7, TrackDyn = FALSE, debug = TRUE,
              debug_freq=100)

set.seed(101)
Cpp_result <- do.call("GilAlgoCpp", c(list(Adj_list), args1))
Cpp_result$FinalStat


set.seed(101)
system.time(R_result <- do.call("GilAlgo", c(list(G), args1)))
R_result$FinalStat

## problem is somewhere in TrackDyn ... 

if (FALSE) {
  source("gillespie.R")
  Rprof("tmp.Rprof")
  set.seed(101)
  system.time(R_result <- GilAlgo(G, N, beta, gamma, MaxTime = 2))
  Rprof(NULL)
  source("https://raw.githubusercontent.com/noamross/noamtools/refs/heads/master/R/proftable.R")
  proftable("tmp.Rprof")
}
