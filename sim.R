library(shellpipes)

manageConflicts()
library(igraph)
library(Rcpp)
## library(RcppClock)

loadEnvironments()
Adj_list <- rdsRead()

# Seed
set.seed(seed)

sourceCpp(matchFile(exts=c("cpp", "Cpp")))

system.time(result <- simFun(Adj_list, N, beta, gamma, MaxTime = 100))

print(result$FinalStat)

rdsSave(result)

