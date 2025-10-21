## short/self-contained tests for Gillespie algorithm code development
library(Rcpp)

## Profiling
#install.packages("RcppClock")
library(RcppClock)

source("NetworkSimulator.R")

TrackDyn <- TRUE
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
sourceCpp('NetSimulator.cpp')
args1 <- list(N, beta, gamma, MaxTime = 3.0
              , TrackDyn = TrackDyn
              , debug = TRUE
              , debug_freq=1
              , debug_low=500, debug_up=550
              #, debug_low=930, debug_up=940
              )

# args2 <- list(50000, beta, gamma, MaxTime = 60.0
#               , TrackDyn = TrackDyn
#               #, debug = TRUE
#               #, debug_freq=1
#               #, debug_low=900, debug_up=5000
#               #, debug_low=930, debug_up=940
# )


set.seed(101)
#set.seed(201)
system.time(Cpp_result <- do.call("GilAlgoCpp", c(list(Adj_list), args1)))
#Cpp_result$FinalStat
print(profile_cpp)
#plot(profile_cpp)
#Cpp_result

set.seed(101)
#set.seed(201)
system.time(R_result <- do.call("GilAlgo", c(list(G), args1)))
#R_result$FinalStat


### R and Cpp sample not consistent for some case: Seed 101
### Solved? Now use a Rcpp:sample wrapper in R version, seems fixed.
# #Cpp_result$Details[526:530,]
# #R_result$Details[526:530,]

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


