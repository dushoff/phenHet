library(shellpipes); manageConflicts()
library(dplyr)

## This function can generate negative numbers;
## should correspond to impossible values of S
sub_zhao1 <- function(x, delta, kappa, kThresh=1e-4){
	kfun <- function(S, k){
		l <- log(S)
		if (k < kThresh)
			return(l + kappa*l^2/2)
		return((exp(kappa*l)-1)/kappa)
	}
   out <- (x^kappa+kfun(x, kappa)/delta)
	if (out <= 0) return(NA)
	return(out)
}

zhao1 <- Vectorize(sub_zhao1)
saveEnvironment()
