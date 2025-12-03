## rm(list=ls())
## wd <- getwd()

### Package Part
library(pracma)
#library(zoo)
library(gsl)
library(deSolve)
#library(igraph)
library(ggplot2)
library(cbinom)
#library(Rcpp)
#library(RcppClock)

source("NetworkSimulator.R")

################################### Network Functions

#### Degree Distribution data frame
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
  
#### PGFs and Derivatives
  
#G_p
PGFG0 <- function(x,Pk){
  G0 <- 0
  for (k in Pk[,1]) {
    G0 <- G0+(x^k)*(Pk[k+1,2])
  }
  return(G0)
}
  
#G'_p
PGFd1G0 <- function(x,Pk){
  d1G0 <- 0
  for (k in c(1:max(Pk[,1]))) {
    d1G0 <- d1G0+(x^(k-1))*k*(Pk[k+1,2])
  }
  return(d1G0)
}
  
#G''_p
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

#G_q
PGFG1 <- function(x,Pk){
  G1 <- PGFd1G0(x,Pk)/Kn(Pk,1)
  return(G1)
}

#G'_q
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
  
### Typical Percolation Process
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
  
  
  
  
#### Miller Slim and Volz
#Configuration model
ModProc_CM <- function(  Pk, beta, gamma
                       , init_omega=1e-3
                       , ODEmaxTime=50
                       , ODEstep=1e-2
                       , ThetaTol=1e-9
                       , TrackDyn=TRUE
                       , init_R=0
                       , s_theta=1e-2){
  if (TrackDyn==TRUE){
    #p0 <- beta/(beta+gamma)
    Sys <- function(t, y, parms){
      with(as.list(c(parms,y)),{
        dtheta <- (-b)*theta+b*PGFd1G0(theta,Pk)/PGFd1G0(1,Pk)+g*(1-theta)
        dR <- g*(1-PGFG0(theta,Pk)-R)
        #dP <- (b+g)*(-P+PGFd1G0(theta,Pk)/PGFd1G0(1,Pk))
        #dP <- (-b)*(PGFd1G0(theta,Pk)/PGFd1G0(1,Pk))+(b+g)*P
        return(list(c(dtheta,dR
                      #,dP
                      ))) 
      }) 
    }
    parms <- c(b=beta,g=gamma)
    times <- seq(0,ODEmaxTime,by=ODEstep)
    y <- c(theta=1-init_omega,R=init_R
           #,P=1
           )
    
    Sys_out <- ode(y,times,Sys,parms)
    S_out <- PGFG0(Sys_out[,2],Pk)
    I_out <- 1-S_out-Sys_out[,3]
    #P_out <- Sys_out[,4]
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
######## Network Functions Part END ###########################################


################################ Initial Condition
# Initial Condition Solver based on I0=1-S0-(R0=0)
Init_omega_func <- function(I0_val,N){
  S0_val <- N-I0_val
  Init_eqn <- function(theta){
    PGFG0(theta,DDist)*N-S0_val
  }
  
  Init_sol <- uniroot(Init_eqn,c(0,1),tol = 1e-9)
  init_theta_val <- 1-Init_sol$root
  return(init_theta_val)
}

##### Eigen Direction R(0)
EigenR <- function(Pk, beta, gamma, lambda, init_omega){
  R_c0 <- beta/(beta+gamma)*PGFd2G0(1,DDist)/lambda
  Eigen_R <- lambda*gamma/(gamma+(beta+gamma)*(R_c0-1))*init_omega
  return(Eigen_R)
}

EigenP <- function(Pk, beta, gamma, lambda, init_omega){
  X <- beta*PGFd2G0(1,Pk)/lambda
  Eigen_P <- X/(X-2*(beta+gamma))*init_omega
  return(Eigen_P)
}

#### Reverse ODE for P
# we want P_inf satisfy the ODE=0 with theta_inf
#sigma_inf <- PGFd1G0(theta_inf,DDist)/lambda
#P_inf <- beta/(beta+gamma)*sigma_inf

#(gamma*(1-PGFG0(theta_inf,DDist)-R_inf))

### Reverse ODE Function
Rvs_ODE <- function(  Pk, beta, gamma
                    , theta_inf
                    , R_inf
                    , ODEmaxTime = 100
                    , ODEstep = 1e-1
                    , disturb = 1e-6){

  P_inf <- beta/(beta+gamma)*PGFd1G0(theta_inf,Pk)/lambda
  lambda <- PGFd1G0(1,Pk)

  Sys <- function(t, y, parms){
    with(as.list(c(parms,y)),{
      dtheta <- -((-b)*theta+b*PGFd1G0(theta,Pk)/l+g*(1-theta))
      dR <- -(g*(1-PGFG0(theta,Pk)-R))
      dP <- +b*PGFd1G0(theta,Pk)/l-(b+g)*P
      return(list(c(dtheta,dR,dP)))
    })
  }
  parms <- c(b=beta,g=gamma,l=lambda)
  times <- seq(0,ODEmaxTime,by=ODEstep)
  y <- c(theta=theta_inf,R=R_inf,P=P_inf)

  Sys_out <- ode(y,times,Sys,parms)
  S_out <- PGFG0(Sys_out[,2],Pk)
  I_out <- 1-S_out-Sys_out[,3]
  P_out <- Sys_out[,4]
  Sys_out <- as.matrix(cbind(Sys_out,S_out,I_out))
  return(Sys_out)
}
