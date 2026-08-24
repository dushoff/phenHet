
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

CheckSeq <- function(Seq){
  return((sum(Seq)%%2==0) && (EG_check(Seq)))
}
