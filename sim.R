system.time(library(shellpipes))

system.time(manageConflicts())
system.time(library(igraph))
system.time(library(Rcpp))

system.time(loadEnvironments())
system.time(Adj_list <- rdsRead())

system.time(set.seed(seed))

system.time(sourceCpp(matchFile(exts=c("cpp", "Cpp"))))

system.time(result <- simFun(Adj_list, N, beta, gamma, MaxTime))

system.time(print(result$FinalStat))

system.time(rdsSave(result))

