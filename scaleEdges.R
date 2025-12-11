### Package Part
library(igraph)
library(cbinom)
library(Rcpp)
## library(RcppClock)

library(shellpipes)
loadEnvironments()
sourceFiles()

#### Disease Parameter
beta <- 0.2
gamma <- 0.25
N <- 1e6
r <- 1
lambda <- 5
MaxTime <- 1000

# Seed
set.seed(2639)

seq <- rnbinom(N,r,mu=lambda)
while(!CheckSeq(seq)){
  seq <- rnbinom(N,r,mu=lambda)
}

CheckSeq(seq)
# generating graph
G <- sample_degseq(  seq
                   , method = "fast.heur.simple"
                   #, method = "configuration.simple"
                   )


# check realization is successful
# should be True
!any(sort(degree(G))-sort(seq)!=0)

#igraph::simple_cycles(G)


# Translate igraph network object into adjacency matrix
Adj_list <- as_adj_list(  G
                        , mode = "all"
                        , loops = "once"
                        , multiple = TRUE
)

sourceCpp(matchFile(exts=c("cpp", "Cpp")))

system.time(result <- simFun(Adj_list, N, beta, gamma, MaxTime))

print(result$FinalStat)

saveEnvironment()

