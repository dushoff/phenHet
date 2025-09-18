#install.packages("igraph")
library(igraph)

# N <- 10000
# delta <- 10
# # Poisson network
# set.seed(15812)
# seq <- rpois(N,delta)

### Check the Even total degree
# sum(seq)%%2

### Checking Realization by Erdos-Gallai Thm
EG_check <- function(DegreeDist){
  check_vec <- c(0)
  DegreeDist <- sort(DegreeDist,decreasing = T)
  N <- length(DegreeDist)
  Cum <- cumsum(DegreeDist)
  
  #Mark3: corrected Durfee number m 
  Dlist <- DegreeDist- c(0:(N-1)) >=0
  m <- length(which(Dlist==TRUE))
  
  for (k in c(1:m)) {
    RHS <- k*(k-1) + k*(length(DegreeDist[which(which(DegreeDist>=k)>k)]))+sum(DegreeDist[which(DegreeDist<k)])
    LHS <- Cum[k]    
    if (LHS<=RHS) {
      check_vec[k] <- 1
    } else {
      check_vec[k] <- 0
    }
  }
  # m<-min(check_vec)
  # return(m)
  return(!any(check_vec==0))
}

# EG_check(seq)

### A bundle check of both even (fast) and EG(formal)
CheckSeq <- function(Seq){
  if((sum(Seq)%%2==0) & (EG_check(Seq))){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# CheckSeq(seq)

# G <- sample_degseq(  seq
#                    , method = "fast.heur.simple"
#                    )

### Perhaps similar to Blitzstein Diaconiz or BKM algorithm???
# The “fast.heur.simple” method generates simple graphs. 
# It is similar to “configuration” but tries to avoid multiple and loop edges 
# and restarts the generation from scratch if it gets stuck. 
# It can generate all simple realizations of a degree sequence,
# but it is not guaranteed to sample them uniformly. 
# This method is relatively fast and it will eventually succeed if the provided 
# degree sequence is graphical, but there is no upper bound 
# on the number of iterations.


# E(G)
# plot(G)
# Check degree
# any(sort(degree(G))-sort(seq)!=0)

# Adj_list <- as_adj_list(  G
#             , mode = "all"
#             , loops = "once"
#             , multiple = TRUE
# )

# as.vector(Adj_list[[1]])
# degree(G)
# sample(c(1:N),5,prob=degree(G))

######### SSA Algorithm for transmission
GilAlgo <- function(  Network
                    , size
                    , beta
                    , gamma
                    , MaxTime
                    , InitInfSize=1
                    , TrackDyn=T
                    ){
  Net <- Network
  G <- as_adj_list(  Net
                   , mode = "all"
                   , loops = "once"
                   , multiple = TRUE
  )
  Deg_vec <- degree(Net)
  N <- size
  g <- gamma
  b <- beta
  
  ind <- c(1:N)
  
  # random initial infection with size i_0
  i_0 <- InitInfSize
  ###Radomly chose s vertices to be infected
  # InitIndex <- c(sample.int(N,i_0))
  
  # Chose infected node with weight of degree
  InitIndex <- sample(c(1:N),i_0,prob=Deg_vec)
  
  #Status: S=0, I=1, R=2
  
  #Initialize status and rate
  t <- 0
  Status <- rep(0,N)
  Status[InitIndex] <- 1
  Istep <- length(which(Status==1))
  
  if (TrackDyn==T){
    NumStep <- 1
    t_vec <- c(t)
    S_vec <- c(length(which(Status==0))/N)
    I_vec <- c(length(which(Status==1))/N)
    R_vec <- c(0)
    
    Infect_time <- rep(NA,N)
    Infect_time[InitIndex] <- 0
    Infect_num <- rep(0,N)
  }
  
  Rate <- rep(0,N)
  Rate[InitIndex] <- g
  
  for (i in c(1:i_0)) {
    x <- InitIndex[i]
    # Network neighbor
    Neighbor <- as.vector(G[[x]])
    # Susceptible neighbor: update their rate
    # Contact <- Neighbor[which(Status[Neighbor]==0)] 
    Contact <- Neighbor[!Status[Neighbor]]
    Rate[Contact] <- b
  }
  cat("Init Sum", sum(Rate),"\n")
  cat("Init index", InitIndex,"\n")
  # while loop: keep looping if t<tmax & Istep != 0 
  # i.e. there is still active infection
  while(t<MaxTime & Istep != 0){
    
    ## SSA Calculation
    Sum <- sum(Rate)
    Cum <- cumsum(Rate)
    # the vertex index of event:
    # cat("Sum is", Sum, ",")
    
    r <- runif(2, min = 0, max = 1)
    Event <- min(which(Cum>r[1]*Sum))
    # Infection: status 0 to 1
    # Recovery: status 1 to 2
    Status[Event] <- Status[Event]+1
    
    # Network neighbor of event index
    Neighbor <- as.vector(G[[Event]])
    
    # Susceptible neighbor: update their rate
    Contact <- Neighbor[!Status[Neighbor]]
    
    # cat("contact: ", Contact,"\n")
    # Infected neighbor: Potential infectors
    Infector <- Neighbor[which(Status[Neighbor]==1)]
    
    ## Time spent for event happen
    Tstep <- -log(r[2])/Sum  
    t <- t+Tstep
    
    # cat("Event index is", Event,",")
    # cat("Status is", Status[Event],"\n")
    # if (Sum<0){
    #   return(Rate)
    #   break
    # }
    ## Update status
    if (Status[Event]==2){               ## Recovery
      Rate[Event] <- 0
      Rate[Contact] <- Rate[Contact]-b
    } else if (Status[Event]==1){        ## Infection
      Rate[Event] <- g
      Rate[Contact] <- Rate[Contact]+b   ## Independence: linear
      
      if (TrackDyn==T){
        # vector of infection time of vertices
        # NA if not being infected eventually
        Infect_time[Event] <- t
        
        # For each infection event in SSA, we might not be able
        # to figure out the exactly one infector as the event is
        # determined by the rate of infectee i.e. number of its
        # actively infected neighbor.
        
        # But since exponential distribution of infection time
        # have the Memorylessness property, and we are assuming all
        # neighbor are iid and considering an expectation, we can average
        # out the new infection event to all active infected neighbor
        # at the moment of event.
        # Infect_num[Infector] <- Infect_num[Infector]+1/(length(Infector))
        
        # As suggested by Ben, we now randomly chose one infector (if more than 
        # one) instead of do the average
        samp_inf <- sample(Infector,1)
        Infect_num[samp_inf] <- Infect_num[samp_inf]+1
      }
    } else {
    }
    
    ## Active number of infections of the whole network
    Istep <- length(which(Status==1))
    
    ## Update proportion
    if (TrackDyn==T){
      NumStep <- NumStep+1
      t_vec[NumStep] <- t
      S_vec[NumStep] <- length(which(Status==0))/N
      I_vec[NumStep] <- length(which(Status==1))/N
      R_vec[NumStep] <- length(which(Status==2))/N
    }
  }
  
  ## Final sizes
  FinishTime <- t
  Ssize <- length(which(Status==0))/N
  Isize <- length(which(Status==1))/N
  Rsize <- length(which(Status==2))/N
  FinalStat <- data.frame(FinishTime,Ssize,Isize,Rsize)
  
  if (TrackDyn==T){
    Track <- cbind(t_vec,S_vec,I_vec,R_vec)
    Infect <- cbind(ind,Infect_time,Infect_num)
    return(list(FinalStat=FinalStat,Details=Track,Reff=Infect))
  } else {
    return(FinalStat) 
  }
}

# beta <- 0.25
# gamma <- 0.2

# result <- GilAlgo(G, N, beta, gamma, MaxTime = 150)
# 
# result$FinalStat
# result$Details
# result$Reff
#max(result$Reff[,3])
