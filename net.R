library(shellpipes)

manageConflicts()
library(igraph)

loadEnvironments()

set.seed(seed)

seq <- rnbinom(N,r,mu=lambda)
while(!CheckSeq(seq)){
  seq <- rnbinom(N,r,mu=lambda)
}
G <- sample_degseq(seq, method = "fast.heur.simple")
stopifnot(!any(sort(degree(G))-sort(seq)!=0))

Adj_list <- as_adj_list(G
	, mode = "all" , loops = "once" , multiple = TRUE
)

rdsSave(Adj_list)
