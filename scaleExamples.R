### Package Part
library(igraph)
library(cbinom)
library(Rcpp)
## library(RcppClock)

library(shellpipes)
sourceFiles()

#### Disease Parameter
beta <- 1
gamma <- 1
N <- 1e6
r <- 1
lambda <- 5

# Seed
set.seed(2639)

nseq <- rnbinom(N,r,mu=lambda)
while(!CheckSeq(nseq)){
  nseq <- rnbinom(N,r,mu=lambda)
}

CheckSeq(nseq)
# generating graph
G <- sample_degseq(  nseq
                   , method = "fast.heur.simple"
                   #, method = "configuration.simple"
                   )


# check realization is successful
stopifnot(!any(sort(degree(G))-sort(nseq)!=0))

#igraph::simple_cycles(G)


# Translate igraph network object into adjacency matrix
Adj_list <- as_adj_list(  G
                        , mode = "all"
                        , loops = "once"
                        , multiple = TRUE
)

## rds(Save Adj_list)
## quit()

### Rcpp Version
sourceCpp(matchFile(exts=c("cpp", "Cpp")))
set.seed(2941)
system.time(Cpp_result <- GilAlgoCpp(Adj_list, N, beta, gamma, MaxTime = 100))
# print(profile_cpp)
Cpp_result$FinalStat
