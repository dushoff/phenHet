### Package Part
library(igraph)
library(cbinom)
library(Rcpp)
## library(RcppClock)

library(shellpipes)
manageConflicts()
sourceFiles()

#### Disease Parameter
gamma <- 0.05
beta <- gamma
MaxTime <- 80
N <- 1e5
r <- 1
lambda <- 5

# Seed
set.seed(2637)

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
system.time(result <- GilAlgoCpp(Adj_list, N, beta, gamma, MaxTime))
# print(profile_cpp)
result$FinalStat

saveEnvironment()
