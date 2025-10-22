## rm(list=ls())
## wd <- getwd()
# install.packages("zoo")

### Package Part
library(pracma)
library(zoo)
library(gsl)
library(deSolve)
library(igraph)
library(ggplot2)
library(cbinom)
library(Rcpp)
library(RcppClock)

source("NetworkSimulator.R")
source("NetworkODE.R")

######################### Parameter Setting
#### Network Parameters
lambda <- 5
kappa <- 2
r <- 1/kappa

(p<-1/(1+kappa*lambda))
(v<-lambda/p)

N <- 50000

kvalue <- seq(0,400)
#Pk <- dpois(kvalue,lambda)
Pk<- dnbinom(kvalue,r,mu=lambda)

# DDist
DDist <- data.frame(kvalue,Pk)

#### Initial Condition
# Initial Condition Solver based on I0=1-S0-(R0=0)
(it_omega <- Init_omega_func(1))
(S0Count <- PGFG0(1-it_omega,DDist)*N)

##### Disease Parameter
beta <- 0.25
gamma <- 0.75
#gamma <- 0.2

#### Eigen Direction R(0)
(Eigen_R <- EigenR(DDist, beta, gamma, lambda, init_omega = it_omega))
######################### Parameter End

######################### MSV ODE part
CM_Opt<- ModProc_CM(  DDist,beta,gamma
                    , ODEmaxTime = 50
                    , ODEstep = 5e-2
                    , init_omega = it_omega
                    , TrackDyn = TRUE
                    , init_R = Eigen_R
                    )
CM_Opt$R0
theta_inf <- CM_Opt$ThetaInfinity
CM_out <- CM_Opt$Dynamic

# Just quickly pull the ODE result for later calculation
time <- CM_out[,1]
theta <- CM_out[,2]
CM_R <- CM_out[,3]
CM_S <- CM_out[,4]
CM_I <- CM_out[,5]

#### Reverse ODE for Todd's p(t)
# manually picking point
theta_inf
tps <- 451
print(CM_out[tps,2])

(P_inf <- beta/(beta+gamma)*PGFd1G0(theta_inf,DDist)/lambda)

CM_out[tps,]
Rvs_vec1<-as.numeric(CM_out[tps,])

{ t_Rvs<-Rvs_vec1[1]
  theta_Rvs <- Rvs_vec1[2]
  R_Rvs <- Rvs_vec1[3]
}

RVS_args <- list(  DDist
                 , beta
                 , gamma
                 , theta_Rvs
                 , R_Rvs
                 , ODEmaxTime = t_Rvs
                 , ODEstep = 5e-2)

Rvs_out <- do.call("Rvs_ODE", c(RVS_args))

Rvs_out[,1]<- max(Rvs_out[,1])-Rvs_out[,1]
Rvs_out <- Rvs_out[order(Rvs_out[,1]),]

Rvs_df<-as.data.frame(Rvs_out)
Rvs_df$I_out[1:30]

#### Check RVS: check with CM ODE for S and I
CM_df <- as.data.frame(CM_out)

ggplot()+theme_bw()+
  geom_point(data=Rvs_df, aes(x=time, y=S_out, color="S_Rvs"), alpha=0.1)+
  geom_line(data=CM_df, aes(x=time, y=S_out,color="S"))+
  geom_point(data=Rvs_df, aes(x=time, y=I_out, color="I_Rvs"), alpha=0.1)+
  geom_line(data=CM_df, aes(x=time, y=I_out,color="I"))+
  xlim(0,10)+
  labs(y = "Proportion")

##################################

### Derivatives based on Negative Binomial
#theta_dot <- -beta*theta+gamma*(1-theta)+beta/lambda*(lambda*CM_S)/(1+kappa*lambda-theta*kappa*lambda)
#S_dot <- theta_dot*(lambda*CM_S)/(1+kappa*lambda-theta*kappa*lambda)

theta_dot <- -beta*theta+gamma*(1-theta)+beta/lambda*PGFd1G0(theta,DDist)
S_dot <- theta_dot*PGFd1G0(theta,DDist)

#def_reff<- -lambda*CM_S*(-(beta+gamma)*(1+log(CM_S)/lambda)+beta*CM_S+gamma)/(CM_I*gamma)
#est_reff<- beta/(gamma)*(lambda-1)*CM_S
#cal_reff<- beta/(beta+gamma)*lambda*CM_S*(1+log(CM_S)/lambda)

R_c0 <- beta/(beta+gamma)*PGFd2G0(1,DDist)/lambda
#R_c0 <- beta/(beta+gamma)*lambda*(kappa+1)

R_imax <- beta/gamma*(PGFd2G0(1,DDist)/lambda-1)
#R_imax <- beta/gamma*(lambda*(kappa+1)-1)

R_i <- -S_dot/(CM_I*gamma)

R_cstar <- beta/(beta+gamma)*PGFd2G0(theta,DDist)/lambda
#R_cstar <- R_c0*CM_S^(1+2*kappa)

PLen<-length(Rvs_df$P)
### E(K_I)
#R_c_plus1 <- Rvs_df$P*(theta[1:PLen]*(PGFd2G0(theta[1:PLen],DDist)/PGFd1G0(theta[1:PLen],DDist))+1)
### E(K_I-1)
R_c <- Rvs_df$P*(PGFd2G0(theta[1:PLen],DDist)/PGFd1G0(theta[1:PLen],DDist))

Rvs_df$R_c <- R_c
#Rvs_df$R_c1<- R_c_plus1

dat_reff <- cbind(time
                  ,R_i
                  #,cal_reff
                  ,R_cstar
                  #,est
                  #,theta
                  #,R_c
                  #,R_c1
                  )


ggplot(data=dat_reff)+theme_bw()+
  geom_line(aes(x=time, y=R_i,color="non-Eigen R_i"))+
  geom_line(aes(x=time, y=R_cstar,color="R_c*"))+
  geom_line(data=Rvs_df,aes(x=time, y=R_c,color="R_c"))+
  #geom_line(data=Rvs_df,aes(x=time, y=R_c1,color="R_c(K+1)"))+
  #geom_line(data=Rvs_df, aes(x=time, y=P*(lambda^2*(kappa+1))/lambda, color="Rev P"))+
  #geom_line(aes(x=time, y=theta, color="theta x10"))+
  geom_hline(yintercept=R_imax,color="black")+
  geom_hline(yintercept=R_c0,color="orange")+
  ylim(0,5)+
  xlim(0,15)+
  #scale_color_manual(values=c("red", "black","brown"))
  labs(y = "R_eff") 


############# Simulation
# Seed
set.seed(2272)

seq <- rnbinom(N,r,mu=lambda)
while(!CheckSeq(seq)){
  seq <- rnbinom(N,r,mu=lambda)
}

CheckSeq(seq)
# generating graph
G <- sample_degseq(  seq
                   , method = "fast.heur.simple"
                   #, method = "configuration.simple"
                   )


# check realization is successful
# should be True
!any(sort(degree(G))-sort(seq)!=0)

#igraph::simple_cycles(G)


# Translate igraph network object into adjacency matrix
Adj_list <- as_adj_list(  G
                        , mode = "all"
                        , loops = "once"
                        , multiple = TRUE
)

### Rcpp Version
sourceCpp('NetSimulator.cpp')
set.seed(2351)
system.time(Cpp_result <- GilAlgoCpp(Adj_list, N, beta, gamma, MaxTime = 100))

# print(profile_cpp)
result <- Cpp_result

## simulation final size
result$FinalStat
## MSV final size
CM_Opt$RInfinity

## simulation dynamic
dat_sim <- result$Details

## comparing I dynamic
ggplot(data = CM_out)+theme_bw()+
  geom_point(data=dat_sim,aes(x=t_vec, y=I_vec,color="Simulation"),alpha=0.1)+
  geom_line(aes(x=time, y=CM_I,color="MSV"))+
  #geom_point(aes(x=time, y=Mod_I,color="Modified"),alpha=0.1)+
  scale_color_manual(values=c("black", "red","blue"))+
  xlim(0,20)+
  labs(y = "I(t)")


## peak values
# random Real infection
max(result$Reff[,6])
# susceptible neighbor at time of infection
max(result$Reff[,5])

### verify calculation: Sum Reff=R_end-1
# random
sum(result$Reff[,6])
result$FinalStat[4]*N
CM_Opt$RInfinity*N
### verified

# visualization
dat_sim_out<-as.data.frame(result$Reff)
dat_Rsim<- dat_sim_out[!is.na(dat_sim_out$Infect_time
),]
#dat_Rsim[1,]
dat_Rsim<-dat_Rsim[order(dat_Rsim$Infect_time),]

## Simulated R_c star: # of Sus Nbr * beta/(beta+gamma)
Sim_RcS <- round(dat_Rsim$S_NbrDeg*(beta/(beta+gamma)),2)

## rolling mean
rn <- 5
edge <- (rn-1)/2
roll_mean <- rep(NA,length(dat_Rsim$Infect_num_rnd))
(roll_mean[c((edge+1):(length(dat_Rsim$Infect_num_rnd)-edge))]
  <-rollmean(dat_Rsim$Infect_num_rnd,rn)
  #<- rollmean(Sim_RcS,rn)
)
dat_Rsim <- cbind(dat_Rsim,Sim_RcS,roll_mean)

#### Visualize
ggplot()+theme_bw()+
  #geom_point(data=dat_Rsim, aes(x=Infect_time, y=Infect_num_rnd,color="Sim RND"),size=0.2)+
  #geom_point(data=dat_Rsim, aes(x=Infect_time, y=Sim_RcS,color="Sim R_c star"),size=0.2)+
  #geom_smooth(data=dat_Rsim, aes(x=Infect_time, y=Infect_num,color="Smooth"))+
  geom_line(data=dat_Rsim, aes(x=Infect_time, y=roll_mean,color="Roll mean n=5"))+
  #geom_line(aes(x=time, y=R_i,color="R_i"))+
  geom_line(data=dat_reff, aes(x=time, y=R_cstar,color="R_c star"))+
  geom_line(data=Rvs_df, aes(x=time, y=R_c,color="R_c"))+
  #geom_line(aes(x=time, y=R_cc,color="R_c"))+
  #geom_line(aes(x=time, y=est, color="Estimation"))+
  #geom_hline(yintercept=beta/(beta+gamma)*lambda,color="purple")+
  #geom_hline(yintercept=peak,color="black")+
  geom_hline(yintercept=R_c0,color="orange")+
  ylim(0,10)+
  xlim(0,10)+
  #scale_color_manual(values=c("red", "black","brown"))
  labs(y = "R_eff") 

dat_Rsim[1:5,]
# Sim_RcS[1:30]
# dat_Rsim$Active_NbrDeg[1:10]
# beta/(beta+gamma)
# dat_Rsim[1:20,c(2:6,8)]
# mean(dat_Rsim$Recovery_time-dat_Rsim$Infect_time)
# mean((dat_Rsim$Infect_num/dat_Rsim$Deg_vec)[1:50])

# all_simple_paths(G,from = 1)

### 20 runs
#set.seed(32025)
# {
#   seq <- rnbinom(N,r,mu=lambda)
#   while(!CheckSeq(seq)){
#     seq <- rnbinom(N,r,mu=lambda)
#   }
#   CheckSeq(seq)
#   G <- sample_degseq(  seq
#                      , method = "fast.heur.simple"
#                      )
#   Adj_list <- as_adj_list(  G
#                             , mode = "all"
#                             , loops = "once"
#                             , multiple = TRUE
#   )
#   system.time(result <- GilAlgoCpp(Adj_list, N, beta, gamma, MaxTime = 100))
#   }
# result$FinalStat
# # CM_Opt$RInfinity
# {
#   dat_sim_out<-as.data.frame(result$Reff)
#   dat_Rsim<- dat_sim_out[!is.na(dat_sim_out$Infect_time),]
#   dat_Rsim<-dat_Rsim[order(dat_Rsim$Infect_time),]
# }
# 
# write.csv2(  dat_Rsim
#           , file="./SimData/sim_50k_g02_round20.csv")

#### readback
{
df1<-read.csv2("./SimData/sim_50k_g02_round01.csv")
#df1$norm_time <- df1$Infect_time-df1$Infect_time[2]
df1$X <- 1

df2<-read.csv2("./SimData/sim_50k_g02_round02.csv")
#df2$norm_time <- df2$Infect_time-df2$Infect_time[2]
df2$X <- 2

df3<-read.csv2("./SimData/sim_50k_g02_round03.csv")
#df3$norm_time <- df3$Infect_time-df3$Infect_time[2]
df3$X <- 3

df4<-read.csv2("./SimData/sim_50k_g02_round04.csv")
#df4$norm_time <- df4$Infect_time-df4$Infect_time[2]
df4$X <- 4

df5<-read.csv2("./SimData/sim_50k_g02_round05.csv")
#df5$norm_time <- df5$Infect_time-df5$Infect_time[2]
df5$X <- 5

df6<-read.csv2("./SimData/sim_50k_g02_round06.csv")
#df6$norm_time <- df6$Infect_time-df6$Infect_time[2]
df6$X <- 6

df7<-read.csv2("./SimData/sim_50k_g02_round07.csv")
#df7$norm_time <- df7$Infect_time-df7$Infect_time[2]
df7$X <- 7

df8<-read.csv2("./SimData/sim_50k_g02_round08.csv")
#df8$norm_time <- df8$Infect_time-df8$Infect_time[2]
df8$X <- 8

df9<-read.csv2("./SimData/sim_50k_g02_round09.csv")
#df9$norm_time <- df9$Infect_time-df9$Infect_time[2]
df9$X <- 9

df10<-read.csv2("./SimData/sim_50k_g02_round10.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df10$X <- 10

df11<-read.csv2("./SimData/sim_50k_g02_round11.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df11$X <- 11

df12<-read.csv2("./SimData/sim_50k_g02_round12.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df12$X <- 12

df13<-read.csv2("./SimData/sim_50k_g02_round13.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df13$X <- 13

df14<-read.csv2("./SimData/sim_50k_g02_round14.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df14$X <- 14

df15<-read.csv2("./SimData/sim_50k_g02_round15.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df15$X <- 15

df16<-read.csv2("./SimData/sim_50k_g02_round16.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df16$X <- 16

df17<-read.csv2("./SimData/sim_50k_g02_round17.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df17$X <- 17

df18<-read.csv2("./SimData/sim_50k_g02_round18.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df18$X <- 18

df19<-read.csv2("./SimData/sim_50k_g02_round19.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df19$X <- 19

df20<-read.csv2("./SimData/sim_50k_g02_round20.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df20$X <- 20

}

### Combine
df_all<-rbind( df1,df2,df3,df4,df5
              ,df6,df7,df8,df9
              #,df10
              #,df11
              ,df12,df13
              #,df14
              ,df15
              ,df16,df17,df18
              #,df19
              ,df20
              )
# For seed 32025, some realization have a long
# warm-up period
df10[1:2,]
df11[1:2,]
df14[1:2,]
df19[1:2,]
# These cause a time phase issue for the sample size 20,
# Creating a second peak for infect numbers due to the phase shifting
# Currently they are withdrawn
# An idea would be "normalize" the time from the second infection.

dat_all<-df_all[order(df_all$Infect_time),]

### Compare the initial condition
mean(dat_all$Infect_num_rnd[1:20])
#R_c0
#dat_all[21:40,c(1:6,8)]

## RcStar sim
Sim_RcS <- round(dat_all$S_NbrDeg*(beta/(beta+gamma)),2)

## P sim?
Psim <- dat_all$Infect_num_rnd/dat_all$Degree
### rolling mean

l <- length(dat_all$Infect_num_rnd)
rn <- 21
edge <- (rn-1)/2
rd<-16

roll_mean <- rep(NA,l)
(roll_mean[c((edge+rd+1):(l-edge))]
  <- rollmean(dat_all$Infect_num_rnd[(rd+1):l],rn)
  #<- rollmean(Sim_RcS[(rd+1):l],rn)
  )
roll_mean[1:rd]<-mean(dat_all$Infect_num_rnd[1:rd])
dat_all <- cbind(dat_all, Sim_RcS, roll_mean, Psim)

### alternative "roll mean"
time<-seq(0,10,0.01)
inf_exp<-time
inf_time<-dat_all$Infect_time
inf_num <-dat_all$Infect_num_rnd
#inf_num <- dat_all$Sim_RcS

len<-0.005

for (i in c(1:length(time))) {
  t <- time[i]
  inf_exp[i]<-mean(inf_num[inf_time>(t-len) & inf_time<(t+len)])
}
dat_exp <- as.data.frame(cbind(time,inf_exp))

ggplot(data=dat_all)+theme_bw()+
  #geom_point(aes(x=Infect_time, y=Infect_num_rnd,color="Sim RND"),size=0.2)+
  #geom_point(aes(x=Infect_time, y=Infect_num_avg,color="Sim AVG"),size=0.2)+
  #geom_point(data=dat_Rsim, aes(x=Infect_time, y=Infect_num_cf,color="Sim CF"),size=0.2)+
  #geom_smooth(data=dat_Rsim, aes(x=Infect_time, y=Infect_num,color="Smooth"))+
  geom_line(aes(x=Infect_time, y=roll_mean,color="Roll mean n=21"),size=0.2,alpha=0.6)+
  #geom_point(aes(x=Infect_time, y=Psim, color="Psim"),size=0.2,alpha=0.6)+
  #geom_point(aes(x=Infect_time, y=roll_mean,color="Roll mean n=5"),size=0.2,alpha=0.2)+
  #geom_line(data=dat_reff,aes(x=time, y=R_i,color="R_i"))+
  #geom_line(aes(x=time, y=cal_reff,color="Zhao1"))+
  geom_line(data=dat_reff,aes(x=time, y=R_cstar,color="R_c star"))+
  geom_line(data=Rvs_df,aes(x=time, y=R_c,color="R_c"))+
  #geom_line(data=dat_reff,aes(x=time, y=R_cc,color="R_c corre"))+
  #geom_line(data=dat_exp,aes(x=time, y=inf_exp,color="exp"))+
  #geom_line(aes(x=time, y=est, color="Estimation"))+
  #geom_hline(yintercept=R_c0*lambda,color="purple")+
  #geom_hline(yintercept=R_imax,color="black")+
  geom_hline(yintercept=R_c0,color="orange")+
  ylim(0,10)+
  xlim(0,20)+
  #scale_color_manual(values=c("red", "black","brown"))
  labs(y = "R_eff") 



# dat_all[1:50,c(1:5,6,10)]
# df20[1:5,c(1:5,6)]

