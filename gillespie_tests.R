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

N <- 50000
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
args1 <- list(N, beta, gamma, MaxTime = 50.0
              , TrackDyn = TrackDyn
              #, debug = TRUE
              #, debug_freq=1
              #, debug_low=900, debug_up=5000
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

### 
prob = degree(G)
set.seed(101)
sample(c(1:N),1,prob=prob)
R_vec <- runif(2000,0,1)
R_vec[1064:1083]
#R_vec[1930:1950]

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


# #Cpp_result$Details[526:530,]
# #R_result$Details[526:530,]
# 
# Cpp_result$Details[500:600,1]-R_result$Details[500:600,1]
# 
# ### RZ test code
# #Cpp_result$Reff[1,]
# #R_result$Reff[1,]
# 
# Cpp_out<-as.data.frame(Cpp_result$Reff)
# R_out <- as.data.frame(R_result$Reff)
# 
# Cpp_sim<- Cpp_out[!is.na(Cpp_out$Infect_time),]
# R_sim<- R_out[!is.na(R_out$Infect_time),]
# 
# Cpp_sim<-Cpp_sim[order(Cpp_sim$Infect_time),]
# R_sim<-R_sim[order(R_sim$Infect_time),]
# 
# Cpp_sim[499:504,2:7]
# R_sim[499:504,2:7]
# 
# idx_cpp <- which(Cpp_sim$Node==40592)
# R_sim[idx_cpp,2:7]
# Cpp_sim[idx_cpp,2:7]
# 
# idx_R <- which(R_sim$ind==56862)
# R_sim[idx_R,2:7]
# Cpp_sim[idx_R,2:7]
# 
# degree(G)[235867]
# 
# Cpp_sim$Degree-R_sim$Deg_vec
# 
