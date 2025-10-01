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

source("NetworkSimulator.R")


################################### Network Functions
{
  #Degree Distribution data frame
  DDistPK <- function(df){
    m <- max(df[,2])
    kvalue <- c(0:m)
    pk <- rep(0,m)
    for (k in kvalue) {
      pk[k+1] <- length(c(which(df[,2]==k)))/nrow(df)
    }
    Pk <- data.frame(kvalue,pk)
    return(Pk)
  }
  # Creat degree distribution frame from raw data df
  
  #PGFs and Derivatives
  
  #G0
  PGFG0 <- function(x,Pk){
    G0 <- 0
    for (k in Pk[,1]) {
      G0 <- G0+(x^k)*(Pk[k+1,2])
    }
    return(G0)
  }
  
  #G'0
  PGFd1G0 <- function(x,Pk){
    d1G0 <- 0
    for (k in c(1:max(Pk[,1]))) {
      d1G0 <- d1G0+(x^(k-1))*k*(Pk[k+1,2])
    }
    return(d1G0)
  }
  
  #G''0
  PGFd2G0 <- function(x,Pk){
    d2G0 <- 0
    m <- max(Pk[,1])
    for (k in c(2:m)) {
      d2G0 <- d2G0+(x^(k-2))*k*(k-1)*(Pk[k+1,2])
    }
    return(d2G0)
  }
  
  #<K^n>
  Kn <- function(Pk,n){
    Knvalue <- 0
    for (k in Pk[,1]) {
      Knvalue <- Knvalue+(k^n)*(Pk[k+1,2])
    }
    return(Knvalue)
  }
  
  #G1
  PGFG1 <- function(x,Pk){
    G1 <- PGFd1G0(x,Pk)/Kn(Pk,1)
    return(G1)
  }
  
  #G'1
  PGFd1G1 <- function(x,Pk){
    G1 <- PGFd2G0(x,Pk)/Kn(Pk,1)
    return(G1)
  }
  
  
  #u_T=G_q(u_T) self contain equation
  ueqn <- function(x) {
    PGFG1(1+(x-1)*Tvalue,Pk_value)-x
  }
  
  
  ##Changing input from beta, gamma to T
  ##Two different type of constant T assumption
  #Newman's concentration assumption
  Tconst_Newman <- function(beta, gamma){
    Tvalue <-1-exp(-beta/gamma)
    return(Tvalue)
  }
  
  #SIR exponential assumption
  Tconst_exp <- function(beta,gamma){
    Tvalue <-beta/(beta+gamma)
    return(Tvalue)
  }
  
  # Typical Percolation Process
  TypProc <- function(Pk,Tvalue,tol=1e-3){
    Tc_value <- 1/(PGFd1G1(1,Pk))
    OBType <- ''
    s <- 0
    Rinfty <- 0
    u <- 0
    v <- 0
    ueqn <- function(x) {
      PGFG1(1+(x-1)*Tvalue,Pk)-x
    }
    if (Tvalue<Tc_value){
      OBType <- 'Limited'
      s <- 1+Tvalue*PGFd1G0(1,Pk)/(1-Tvalue*PGFd1G1(1,Pk))
      Rinfty <- 0
      u <- 1
      v <- 1
    }
    else if(Tvalue==Tc_value){
      OBType <- 'Undefined'
      s <- 0
      Rinfty <- 0
    }
    else if(Tvalue>Tc_value){
      OBType <- 'Epidemic'
      s <- 0
      #usol <- uniroot(ueqn,c(0+tol,1-tol),tol = 1e-11)
      
      LB_u <- 0
      UB_u <- 1
      u_vec <- seq(from=LB_u,to=UB_u,by=tol)
      u_mat <- matrix(0,nrow = length(u_vec),ncol = 3)
      u_mat[,1] <- u_vec
      u0 <- ueqn(0) 
      
      for (i in c(1:length(u_vec))) {
        u_mat[i,2] <- ueqn(u_vec[i])
        u_mat[i,3] <- u_mat[i,2]/u0
      }
      
      if (length(which(u_mat[,3]<0)) == 0){
        u <- 1
      }
      else{
        UB_u <- u_mat[min(which(u_mat[,3]<0)),1]
        usol <- uniroot(ueqn,c(LB_u,UB_u),tol = 1e-10)
        u <- usol$root
      }
      v <- 1-Tvalue+Tvalue*u
      Rinfty <- 1-PGFG0(v,Pk)
    }
    
    eta <- 0
    m <- max(Pk[,1])
    for (j in c(0:m)) {
      eta <- eta+Pk[,2][which(Pk[,1]==j)]*(v^j)
    }
    
    OutputDF <- data.frame(Tvalue,Tc_value,OBType,s,Rinfty,u,v,eta)
    return(OutputDF)
  }
  
  
  
  
  #Miller Slim and Voltz
  #Configuration model
  ModProc_CM <- function(Pk, beta, gamma, init_theta=1e-3, ODEmaxTime=50, ODEstep=1e-2,ThetaTol=1e-9, TrackDyn=TRUE,init_R=0, s_theta=1e-2){
    if (TrackDyn==TRUE){
      #p0 <- beta/(beta+gamma)
      Sys <- function(t, y, parms){
        with(as.list(c(parms,y)),{
          dtheta <- (-b)*theta+b*PGFd1G0(theta,Pk)/PGFd1G0(1,Pk)+g*(1-theta)
          dR <- g*(1-PGFG0(theta,Pk)-R)
          dP <- (b+g)*(-P+PGFd1G0(theta,Pk)/PGFd1G0(1,Pk))
          #dP <- (-b)*(PGFd1G0(theta,Pk)/PGFd1G0(1,Pk))+(b+g)*P
          return(list(c(dtheta,dR,dP))) 
        }) 
      }
      parms <- c(b=beta,g=gamma)
      times <- seq(0,ODEmaxTime,by=ODEstep)
      y <- c(theta=1-init_theta,R=init_R,P=1)
      
      Sys_out <- ode(y,times,Sys,parms)
      S_out <- PGFG0(Sys_out[,2],Pk)
      I_out <- 1-S_out-Sys_out[,3]
      P_out <- Sys_out[,4]
      Sys_out <- as.matrix(cbind(Sys_out,S_out,I_out))
    }
    
    g <- gamma
    b <- beta
    
    thetaEqn<- function(x) {
      g/(b+g)+b/(b+g)*PGFd1G0(x,Pk)/PGFd1G0(1,Pk)-x
    }
    
    LB_theta <- 0
    UB_theta <- 1
    step_theta <- s_theta
    
    Btheta_vec <- seq(from=LB_theta,to=UB_theta,by=step_theta)
    Btheta_mat <- matrix(0,nrow = length(Btheta_vec),ncol = 3)
    Btheta_mat[,1] <- Btheta_vec
    theta0 <- thetaEqn(0) 
    
    for (i in c(1:length(Btheta_vec))) {
      Btheta_mat[i,2] <- thetaEqn(Btheta_vec[i])
      Btheta_mat[i,3] <- Btheta_mat[i,2]/theta0
    }
    
    if (length(which(Btheta_mat[,3]<0)) == 0){
      thetaInf <- 1
    }else{
      UB_theta <- Btheta_mat[min(which(Btheta_mat[,3]<0)),1]
      theta_sol <- uniroot(thetaEqn,c(LB_theta,UB_theta),tol = 1e-9)
      thetaInf <- theta_sol$root
    }
    
    R0 <- b/(b+g)*PGFd2G0(1,Pk)/PGFd1G0(1,Pk)
    RInf <- 1-PGFG0(thetaInf,Pk)
    
    if (TrackDyn==TRUE){
      return(list(R0=R0,RInfinity=RInf, ThetaInfinity=thetaInf, Dynamic=Sys_out))
    } else {
      return(list(R0=R0,RInfinity=RInf, ThetaInfinity=thetaInf))
    }
  }
  
}
#################### Network Fuctions Part END ###########################################

################################ Distribution part####################################################
#### Distribution Parameter
lambda <- 5
kappa <- 2
r <- 1/kappa

(p<-1/(1+kappa*lambda))
(v<-lambda/p)

#### Disease Parameter
beta <- 0.25
gamma <- 0.75
# gamma <- 0.2

N <- 50000

# dnbinom(10,r,mu=lambda)

kvalue <- seq(0,400)
#Pk <- dpois(kvalue,lambda)
Pk<- dnbinom(kvalue,r,mu=lambda)

DDist <- data.frame(kvalue,Pk)

PGFd2G0(1,DDist)/lambda*beta/(beta+gamma)
# DDist
################################ Distribution part END ##############

################################ Initial Condition
# Initial Condition Solver based on I0=1-S0-(R0=0)
Init_theta_func <- function(I0_val){
  S0_val <- N-I0_val
  Init_eqn <- function(theta){
    PGFG0(theta,DDist)*N-S0_val
  }
  
  Init_sol <- uniroot(Init_eqn,c(0,1),tol = 1e-9)
  init_theta_val <- 1-Init_sol$root
  return(init_theta_val)
}
it_theta <- Init_theta_func(1)
S0Count <- PGFG0(1-it_theta,DDist)*N


#### Fully mixed/Mass Action SIR Model
# MASIR_Proc <- function(b,g,lambda,init_S=1e-3, ODEmaxTime=50, ODEstep=1e-2,TrackDyn=TRUE){
#   if (TrackDyn==TRUE){
#     Sys <- function(t, y, parms){
#       with(as.list(c(parms,y)),{
#         dS <- (-b*lambda)*X*S
#         dX <- (b*lambda)*X*S-(g+b)*X
#         dI <- (b*lambda)*X*S-(g)*I
#         dR <- g*I
#         return(list(c(dS,dI,dR,dX)))
#       })
#     }
#     parms <- c(b=b,g=g,lambda=lambda)
#     times <- seq(0,ODEmaxTime,by=ODEstep)
#     y <- c(S=init_S,I=1-init_S,R=0,X=(1-init_S))
#     
#     Sys_out <- ode(y,times,Sys,parms)
#   }
#   
#   S_0 <- init_S
#   R_0 <- 0
#   
#   R0 <- b/g
#   RInf <- Sys_out[length(Sys_out[,4]),4]
#   
#   if (TrackDyn==TRUE){
#     return(list(R0=R0,RInfinity=RInf, Dynamic=Sys_out))
#   } else {
#     return(list(R0=R0,RInfinity=RInf))
#   }
# }

#### Network Model
CM_Opt<- ModProc_CM(DDist,beta,gamma,ODEmaxTime = 200, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = TRUE)
#MA_Opt<- MASIR_Proc(beta, gamma, lambda, init_S = (N-1)/N, ODEmaxTime=100, ODEstep=1e-1,TrackDyn = TRUE)
#Mod_Opt<- MAmod_Proc(beta, gamma, lambda, init_S = (N-1)/N, ODEmaxTime=100, ODEstep=1e-1,TrackDyn = TRUE)
CM_Opt$R0


#1+log((N-1)/N)/lambda

CM_out <- CM_Opt$Dynamic
# CM_out[10,]
#MA_out <- MA_Opt$Dynamic
#Mod_out <- Mod_Opt$Dynamic


time <- CM_out[,1]
CM_I <- CM_out[,6]
#MA_I <- MA_out[,3]
#Mod_I <- Mod_out[,3]

CM_R <- CM_out[,3]
#MA_R <- MA_out[,4]
#Mod_R <- Mod_out[,4]

CM_P <- CM_out[,4]
theta <- CM_out[,2]

CM_S <- CM_out[,5]
#MA_S <- MA_out[,2]
#Mod_S <- Mod_out[,2]

dat_S <- cbind(time,CM_S
               #, MA_S
               #, Mod_S
               )
dat_R <- cbind(time,CM_R
               #, MA_R
               #, Mod_R
               )
dat <- cbind(time,CM_I
             #, MA_I
             #, Mod_I
             )

 
#CM_P

#Mod_S
#Mod_I

### Derivatives based on Negative Binomial
theta_dot <- -beta*theta+gamma*(1-theta)+beta/lambda*(lambda*CM_S)/(1+kappa*lambda-theta*kappa*lambda)
S_dot <- theta_dot*(lambda*CM_S)/(1+kappa*lambda-theta*kappa*lambda)

#def_reff<- -lambda*CM_S*(-(beta+gamma)*(1+log(CM_S)/lambda)+beta*CM_S+gamma)/(CM_I*gamma)
#est_reff<- beta/(gamma)*(lambda-1)*CM_S
#cal_reff<- beta/(beta+gamma)*lambda*CM_S*(1+log(CM_S)/lambda)
# new_reff<- beta/(beta+gamma)*lambda*CM_S
#new_reff<- beta/(gamma)*(lambda-1)*CM_S

R_c0 <- beta/(beta+gamma)*(lambda^2*(kappa+1))/lambda
peak <- beta/gamma*(lambda*(kappa+1)-1)

R_i <- -S_dot/(CM_I*gamma)
R_c <- R_c0*CM_S^(1+2*kappa)
est <- peak*CM_S^(1+2*kappa)
R_cc <- R_c*CM_P
#R_c1 <- beta/(beta+gamma)*CM_P*PGFd2G0(theta,DDist)/PGFd1G0(theta,DDist)


#new <- beta/gamma*(lambda*(kappa+1)*(2*CM_S^(2*kappa+1)-CM_S^(kappa)*(theta-gamma/beta*(1-theta)))-1)

dat_reff <- cbind(time,R_i
                  #,cal_reff
                  ,R_c
                  ,est
                  #,new
                  ,theta
                  ,R_cc
                  ,R_c1
                  )

ggplot(data=dat_reff)+theme_bw()+
  #geom_line(aes(x=time, y=R_i,color="Instantaneous"))+
  #geom_line(aes(x=time, y=cal_reff,color="Zhao1"))+
  geom_line(aes(x=time, y=R_c,color="Case with no correction"))+
  geom_line(aes(x=time, y=R_cc, color="Corrected Rc"))+
  geom_line(aes(x=time, y=theta, color="theta x10"))+
  #geom_hline(yintercept=beta/(beta+gamma)*lambda,color="purple")+
  #geom_hline(yintercept=peak,color="black")+
  geom_hline(yintercept=R_c0,color="orange")+
  ylim(0,5)+
  xlim(0,20)+
  #scale_color_manual(values=c("red", "black","brown"))
  labs(y = "R_eff") 

CM_P[1:20]
time[1:20]
############# Simulation

# Seed
set.seed(3528)

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
# Adj_list <- as_adj_list(  G
#                         , mode = "all"
#                         , loops = "once"
#                         , multiple = TRUE
# )

result <- GilAlgo(G, N, beta, gamma, MaxTime = 100)

## simulation final size
result$FinalStat
## MSV final size
CM_Opt$RInfinity

## simulation dynamic
dat_sim <- result$Details
## MSV dynamic
# CM_out

## comparing I dynamic
ggplot(data = dat)+theme_bw()+
  geom_point(data=dat_sim,aes(x=t_vec, y=I_vec,color="Simulation"),alpha=0.1)+
  geom_line(aes(x=time, y=CM_I,color="MSV"))+
  #geom_point(aes(x=time, y=Mod_I,color="Modified"),alpha=0.1)+
  scale_color_manual(values=c("black", "red","blue"))+
  xlim(0,60)+
  labs(y = "I(t)")


### Reff
result$Reff[1,]
## peak values
# random Real infection
max(result$Reff[,6])
# susceptible neighbor at time of infection
max(result$Reff[,5])


## avg
# max(result$Reff[,7])
## counter factural
# max(result$Reff[,8])
#peak

### verify calculation: Sum Reff=R_end-1
# random
sum(result$Reff[,6])
# avg
#sum(result$Reff[,7])

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
rn <- 9
edge <- (rn-1)/2
roll_mean <- rep(NA,length(dat_Rsim$Infect_num_rnd))
(roll_mean[c((edge+1):(length(dat_Rsim$Infect_num_rnd)-edge))]
  <-rollmean(dat_Rsim$Infect_num_rnd,rn)
  #<- rollmean(Sim_RcS,rn)
)
dat_Rsim <- cbind(dat_Rsim,Sim_RcS,roll_mean)

#### Visualize
ggplot(data=dat_reff)+theme_bw()+
  geom_point(data=dat_Rsim, aes(x=Infect_time, y=Infect_num_rnd,color="Sim RND"),size=0.2)+
  #geom_point(data=dat_Rsim, aes(x=Infect_time, y=Sim_RcS,color="Sim R_c star"),size=0.2)+
  #geom_point(data=dat_Rsim, aes(x=Infect_time, y=Infect_num_avg,color="Sim AVG"),size=0.2)+
  #geom_point(data=dat_Rsim, aes(x=Infect_time, y=Infect_num_cf,color="Sim CF"),size=0.2)+
  #geom_smooth(data=dat_Rsim, aes(x=Infect_time, y=Infect_num,color="Smooth"))+
  geom_line(data=dat_Rsim, aes(x=Infect_time, y=roll_mean,color="Roll mean n=5"))+
  #geom_line(aes(x=time, y=R_i,color="R_i"))+
  #geom_line(aes(x=time, y=cal_reff,color="Zhao1"))+
  geom_line(aes(x=time, y=R_c,color="R_c star"))+
  #geom_line(aes(x=time, y=R_cc,color="R_c"))+
  #geom_line(aes(x=time, y=est, color="Estimation"))+
  #geom_hline(yintercept=beta/(beta+gamma)*lambda,color="purple")+
  #geom_hline(yintercept=peak,color="black")+
  geom_hline(yintercept=R_c0,color="orange")+
  ylim(0,10)+
  xlim(0,10)+
  #scale_color_manual(values=c("red", "black","brown"))
  labs(y = "R_eff") 

# Sim_RcS[1:30]
# dat_Rsim$Active_NbrDeg[1:10]
# beta/(beta+gamma)
dat_Rsim[11:20,c(2:6,8)]
mean(dat_Rsim$Recovery_time-dat_Rsim$Infect_time)
#mean((dat_Rsim$Infect_num/dat_Rsim$Deg_vec)[1:50])

# all_simple_paths(G,from = 1)

### 20 runs
#set.seed(32025)

{
  seq <- rnbinom(N,r,mu=lambda)
  while(!CheckSeq(seq)){
    seq <- rnbinom(N,r,mu=lambda)
  }
  CheckSeq(seq)
  G <- sample_degseq(  seq
                     , method = "fast.heur.simple"
                     )
  result <- GilAlgo(G, N, beta, gamma, MaxTime = 100)
  dat_sim_out<-as.data.frame(result$Reff)
  dat_Rsim<- dat_sim_out[!is.na(dat_sim_out$Infect_time),]
  dat_Rsim<-dat_Rsim[order(dat_Rsim$Infect_time),]
}
result$FinalStat

write.csv2(  dat_Rsim
             , file="./SimData/sim_50k_md5_round20.csv")

#### readback
{
df1<-read.csv2("./SimData/sim_50k_md5_round01.csv")
#df1$norm_time <- df1$Infect_time-df1$Infect_time[2]
df1$X <- 1

df2<-read.csv2("./SimData/sim_50k_md5_round02.csv")
#df2$norm_time <- df2$Infect_time-df2$Infect_time[2]
df2$X <- 2

df3<-read.csv2("./SimData/sim_50k_md5_round03.csv")
#df3$norm_time <- df3$Infect_time-df3$Infect_time[2]
df3$X <- 3

df4<-read.csv2("./SimData/sim_50k_md5_round04.csv")
#df4$norm_time <- df4$Infect_time-df4$Infect_time[2]
df4$X <- 4

df5<-read.csv2("./SimData/sim_50k_md5_round05.csv")
#df5$norm_time <- df5$Infect_time-df5$Infect_time[2]
df5$X <- 5

df6<-read.csv2("./SimData/sim_50k_md5_round06.csv")
#df6$norm_time <- df6$Infect_time-df6$Infect_time[2]
df6$X <- 6

df7<-read.csv2("./SimData/sim_50k_md5_round07.csv")
#df7$norm_time <- df7$Infect_time-df7$Infect_time[2]
df7$X <- 7

df8<-read.csv2("./SimData/sim_50k_md5_round08.csv")
#df8$norm_time <- df8$Infect_time-df8$Infect_time[2]
df8$X <- 8

df9<-read.csv2("./SimData/sim_50k_md5_round09.csv")
#df9$norm_time <- df9$Infect_time-df9$Infect_time[2]
df9$X <- 9

df10<-read.csv2("./SimData/sim_50k_md5_round10.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df10$X <- 10

df11<-read.csv2("./SimData/sim_50k_md5_round11.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df11$X <- 11

df12<-read.csv2("./SimData/sim_50k_md5_round12.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df12$X <- 12

df13<-read.csv2("./SimData/sim_50k_md5_round13.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df13$X <- 13

df14<-read.csv2("./SimData/sim_50k_md5_round14.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df14$X <- 14

df15<-read.csv2("./SimData/sim_50k_md5_round15.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df15$X <- 15

df16<-read.csv2("./SimData/sim_50k_md5_round16.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df16$X <- 16

df17<-read.csv2("./SimData/sim_50k_md5_round17.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df17$X <- 17

df18<-read.csv2("./SimData/sim_50k_md5_round18.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df18$X <- 18

df19<-read.csv2("./SimData/sim_50k_md5_round19.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df19$X <- 19

df20<-read.csv2("./SimData/sim_50k_md5_round20.csv")
#df10$norm_time <- df10$Infect_time-df10$Infect_time[2]
df20$X <- 20

}

### Combine
df_all<-rbind( df1,df2,df3,df4,df5
              ,df6,df7,df8,df9,df10
              ,df11,df12,df13,df14,df15
              ,df16,df17,df18
              ,df19
              ,df20
              )
dat_all<-df_all[order(df_all$Infect_time),]

### Compare the initial condition
mean(dat_all$Infect_num_rnd[1:20])
R_c0
dat_all[1:20,c(1:6,8)]

### rolling mean
Sim_RcS <- round(dat_all$S_NbrDeg*(beta/(beta+gamma)),2)

l <- length(dat_all$Infect_num_rnd)
rn <- 9
edge <- (rn-1)/2
rd<-20

roll_mean <- rep(NA,l)
(roll_mean[c((edge+rd+1):(l-edge))]
  #<- rollmean(dat_all$Infect_num_rnd[(rd+1):l],rn)
  <- rollmean(Sim_RcS[(rd+1):l],rn)
  )
roll_mean[1:rd]<-mean(dat_all$Infect_num_rnd[1:rd])
dat_all <- cbind(dat_all, Sim_RcS, roll_mean)

### alternative "roll mean"
time<-seq(0,10,0.01)
inf_exp<-time
inf_time<-dat_all$Infect_time
#inf_num <-dat_all$Infect_num_rnd
inf_num <- dat_all$Sim_RcS

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
  #geom_line(aes(x=Infect_time, y=roll_mean,color="Roll mean n=5"),size=0.2,alpha=0.6)+
  #geom_point(aes(x=Infect_time, y=roll_mean,color="Roll mean n=5"),size=0.2,alpha=0.2)+
  #geom_line(data=dat_reff,aes(x=time, y=R_i,color="R_i"))+
  #geom_line(aes(x=time, y=cal_reff,color="Zhao1"))+
  geom_line(data=dat_reff,aes(x=time, y=R_c,color="R_c star"))+
  #geom_line(data=dat_reff,aes(x=time, y=R_cc,color="R_c corre"))+
  geom_line(data=dat_exp,aes(x=time, y=inf_exp,color="exp"))+
  #geom_line(aes(x=time, y=est, color="Estimation"))+
  #geom_hline(yintercept=beta/(beta+gamma)*lambda,color="purple")+
  #geom_hline(yintercept=peak,color="black")+
  geom_hline(yintercept=R_c0,color="orange")+
  ylim(0,20)+
  xlim(0,10)+
  #scale_color_manual(values=c("red", "black","brown"))
  labs(y = "R_eff") 



# dat_all[1:50,c(1:5,6,10)]
# df20[1:5,c(1:5,6)]

