### Package Part
library(igraph)
library(cbinom)
library(Rcpp)
## library(RcppClock)

library(shellpipes)
sourceFiles()

#### Disease Parameter
beta <- 0.25
gamma <- 0.2
N <- 40000
r <- 1
lambda <- 5

# Seed
set.seed(2639)

r
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

### Rcpp Version
sourceCpp(matchFile(exts=c("cpp", "Cpp")))
set.seed(2941)
system.time(Cpp_result <- GilAlgoCpp(Adj_list, N, beta, gamma, MaxTime = 100))
# print(profile_cpp)
Cpp_result$FinalStat
