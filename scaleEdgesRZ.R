### Package Part
library(igraph)
library(cbinom)
library(Rcpp)
library(zoo)
## library(RcppClock)

# library(shellpipes)
# loadEnvironments()
# sourceFiles()

source("scaleFuns.R")
source("NetworkODE.R")
sourceCpp("edgelist.cpp")

#### Disease Parameter
beta <- 0.2
gamma <- 0.25
N <- 1e6
r <- 1
lambda <- 5

# Seed
set.seed(2639)

nseq <- rnbinom(N,r,mu=lambda)
while(!CheckSeq(nseq)){
  nseq <- rnbinom(N,r,mu=lambda)
}

CheckSeq(nseq)
# generating graph
G <- sample_degseq(  nseq
                   , method = "fast.heur.simple"
                   #, method = "configuration.simple"
                   )


# check realization is successful
# should be True
stopifnot(!any(sort(degree(G))-sort(nseq)!=0))

#igraph::simple_cycles(G)


# Translate igraph network object into adjacency matrix
Adj_list <- as_adj_list(  G
                        , mode = "all"
                        , loops = "once"
                        , multiple = TRUE
)

# sourceCpp(matchFile(exts=c("cpp", "Cpp")))

system.time(result <- simFun(Adj_list, N, beta, gamma, MaxTime = 100))

print(result$FinalStat)


# ### Order R_c for infect time
dat_sim_out<-result$Infector
dat_Rsim<- dat_sim_out[!is.na(dat_sim_out$InfectTime),]
dat_Rsim<-dat_Rsim[order(dat_Rsim$InfectTime),]

## rolling mean
rn <- 5
edge <- (rn-1)/2
roll_mean <- rep(NA,length(dat_Rsim$NumInfected))
(roll_mean[c((edge+1):(length(dat_Rsim$NumInfected)-edge))]
  <-rollmean(dat_Rsim$NumInfected,rn)
  #<- rollmean(Sim_RcS,rn)
)
dat_Rsim <- cbind(dat_Rsim,roll_mean)


### R_0 calculation

kvalue <- seq(0,400)
#Pk <- dpois(kvalue,lambda)
Pk<- dnbinom(kvalue,r,mu=lambda)
# DDist
DDist <- data.frame(kvalue,Pk)
R_c0 <- beta/(beta+gamma)*PGFd2G0(1,DDist)/lambda
R_imax <- beta/gamma*(PGFd2G0(1,DDist)/lambda-1)


ggplot()+theme_bw()+
  #geom_point(data=dat_Rsim, aes(x=Infect_time, y=Infect_num_rnd,color="Sim RND"),size=0.2)+
  #geom_point(data=dat_Rsim, aes(x=Infect_time, y=Sim_RcS,color="Sim R_c star"),size=0.2)+
  #geom_line(data=dat_Rsim, aes(x=InfectTime, y=NumInfected,color="Sim"))+
  geom_line(data=dat_Rsim, aes(x=InfectTime, y=roll_mean,color="Roll mean n=5"))+
  #geom_line(aes(x=time, y=R_i,color="R_i"))+
  #geom_line(data=dat_reff, aes(x=time, y=R_cstar,color="R_c star"))+
  #geom_line(data=Rvs_df, aes(x=time, y=R_c,color="R_c"))+
  #geom_line(aes(x=time, y=R_cc,color="R_c"))+
  #geom_line(aes(x=time, y=est, color="Estimation"))+
  #geom_hline(yintercept=beta/(beta+gamma)*lambda,color="purple")+
  #geom_hline(yintercept=peak,color="black")+
  geom_hline(yintercept=R_c0,color="orange")+
  ylim(0,10)+
  xlim(0,20)+
  #scale_color_manual(values=c("red", "black","brown"))
  labs(y = "R_eff")

# 
# 
# #betasaveEnvironment()
# 
