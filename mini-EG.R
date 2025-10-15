library(Rcpp)
library(igraph)

sourceCpp("mini-EG.cpp")
#source("NetworkSimulator.R")

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
# CheckSeq(seq)
# while(!CheckSeq(seq)){
#   seq <- rnbinom(N,r,mu=lambda)
# }
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

prob <- degree(G)
x <- c(29701,54212,177425)
#x <- c(2,3,4)

### Baseline
set.seed(101)
sample(c(1:N),1,prob=prob)
B_vec <- runif(1200,0,1)
B_vec[1074:1080]

### R::sample
set.seed(101)
sample(c(1:N),1,prob=prob)
R_vec <- runif(1074,0,1)
R_vec[1074]
sample(x,1)
runif(2,0,1)

### Cpp::sample
set.seed(101)
sample(c(1:N),1,prob=prob)
C_vec <- Rrunif_eg(1074)
C_vec[1074]
CppSample_eg(x,1)
Rrunif_eg(2)

### Some weird numerical issue just happen at this 
### specific value!!!!!!!!!!!!!!!!!!
### Depends on seeds 101 and event 527
### Independent on x vector value

### A more consistent sample function in R??

#set.seed(101)
#sample(c(1:N),1,prob=prob)
#Cpp_vec <- Rrunif_eg(1500)
#Cpp_vec-R_vec
