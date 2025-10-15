library(Rcpp)

sourceCpp("mini-EG.cpp")
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
prob <- degree(G)

x <- c(6,3,5,7)
y <- c(1,2)


### Rcpp pure unif draw
set.seed(101)
Rrunif_eg(4)

### R pure unif draw
set.seed(101)
runif(4,0,1)

### Rcpp unif+sample x
set.seed(101)
Rrunif_eg(1)
CppSample_eg(x)
CppSample_eg(x)
CppSample_eg(x)
Rrunif_eg(2)

### R unif+sample x
set.seed(101)
runif(1,0,1)
sample(x,1)
sample(x,1)
sample(x,1)
runif(2,0,1)

