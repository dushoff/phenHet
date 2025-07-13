zhao1 <- function(x, delta, kappa, kThresh=1e-4){
	kfun <- function(S, k){
		l <- log(S)
		if (k < kThresh)
			return(l + kappa*l^2/2)
		return((exp(kappa*l)-1)/kappa)
	}
   return(x^kappa+kfun(x, kappa)/delta)
}
