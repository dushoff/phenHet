library(shellpipes)

manageConflicts()
library(igraph)
library(Rcpp)

loadEnvironments()

set.seed(seed)

system.time(seq <- rnbinom(N,r,mu=lambda))
system.time(while(!CheckSeq(seq)){
  seq <- rnbinom(N,r,mu=lambda)
})
system.time(G <- sample_degseq(seq, method = "fast.heur.simple"))
stopifnot(!any(sort(degree(G))-sort(seq)!=0))

system.time(Adj_list <- as_adj_list(G
	, mode = "all" , loops = "once" , multiple = TRUE
))

system.time(sourceCpp(matchFile(exts=c("cpp", "Cpp"))))

system.time(result <- simFun(Adj_list, N, beta, gamma, MaxTime))

print(result$FinalStat)

system.time(rdsSave(result))

