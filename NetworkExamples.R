rm(list=ls())

wd <- getwd()
### Package Part
library(pracma)
library(gsl)
library(deSolve)
library(igraph)
library(ggplot2)
library(rriskDistributions)
library(tidyr)
library(cbinom)
################################### Network Part
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
  ModProc_CM <- function(Pk, beta, gamma, init_theta=1e-3, ODEmaxTime=50, ODEstep=1e-2,ThetaTol=1e-9, TrackDyn=T,init_R=0, s_theta=1e-2){
    if (TrackDyn==T){
      Sys <- function(t, y, parms){
        with(as.list(c(parms,y)),{
          dtheta <- (-b)*theta+b*PGFd1G0(theta,Pk)/PGFd1G0(1,Pk)+g*(1-theta)
          dR <- g*(1-PGFG0(theta,Pk)-R)
          return(list(c(dtheta,dR))) 
        }) 
      }
      parms <- c(b=beta,g=gamma)
      times <- seq(0,ODEmaxTime,by=ODEstep)
      y <- c(theta=1-init_theta,R=init_R)
      
      Sys_out <- ode(y,times,Sys,parms)
      S_out <- PGFG0(Sys_out[,2],Pk)
      I_out <- 1-S_out-Sys_out[,3]
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
    
    if (TrackDyn==T){
      return(list(R0=R0,RInfinity=RInf, ThetaInfinity=thetaInf, Dynamic=Sys_out))
    } else {
      return(list(R0=R0,RInfinity=RInf, ThetaInfinity=thetaInf))
    }
  }
  
}
#################### Network Part END ###########################################

################################ Distribution part####################################################
lambda <- 10
kvalue <- seq(0,500)
Pk <- dpois(kvalue,lambda)
DDist <- data.frame(kvalue,Pk)
# DDist
################################ Distribution part END ##############


################################ Initial Condition
N <- 500000

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
MAmod_Proc <- function(beta, gamma,lambda, init_S=1-(1e-3), ODEmaxTime=50, ODEstep=1e-2,TrackDyn=T){
  R_net <- beta*gamma*lambda/(beta+gamma)
  if (TrackDyn==T){
    Sys <- function(t, y, parms){
      with(as.list(c(parms,y)),{
        dS <- lambda*S*(-(b+g)*(1+log(S)/lambda)+b*S+g)
        dI <- -lambda*S*(-(b+g)*(1+log(S)/lambda)+b*S+g)-g*I
        dR <- g*I
        return(list(c(dS,dI,dR)))
      })
    }
    parms <- c(b=beta,g=gamma,R_net=R_net, lambda=lambda)
    times <- seq(0,ODEmaxTime,by=ODEstep)
    y <- c(S=init_S,I=1-init_S,R=0)
    
    Sys_out <- ode(y,times,Sys,parms)
  }
  
  g <- gamma
  b <- beta
  S_0 <- init_S
  R_0 <- 0
  
  R0 <- b/g
  RInf <- Sys_out[length(Sys_out[,4]),4]
  
  if (TrackDyn==T){
    return(list(R0=R0,RInfinity=RInf,Rnet=R_net, Dynamic=Sys_out))
  } else {
    return(list(R0=R0,RInfinity=RInf))
  }
}

MASIR_Proc <- function(b,g,lambda,init_S=1e-3, ODEmaxTime=50, ODEstep=1e-2,TrackDyn=T){
  if (TrackDyn==T){
    Sys <- function(t, y, parms){
      with(as.list(c(parms,y)),{
        dS <- (-b*lambda)*X*S
        dX <- (b*lambda)*X*S-(g+b)*X
        dI <- (b*lambda)*X*S-(g)*I
        dR <- g*I
        return(list(c(dS,dI,dR,dX)))
      })
    }
    parms <- c(b=b,g=g,lambda=lambda)
    times <- seq(0,ODEmaxTime,by=ODEstep)
    y <- c(S=init_S,I=1-init_S,R=0,X=(1-init_S))
    
    Sys_out <- ode(y,times,Sys,parms)
  }
  
  S_0 <- init_S
  R_0 <- 0
  
  R0 <- b/g
  RInf <- Sys_out[length(Sys_out[,4]),4]
  
  if (TrackDyn==T){
    return(list(R0=R0,RInfinity=RInf, Dynamic=Sys_out))
  } else {
    return(list(R0=R0,RInfinity=RInf))
  }
}

beta <- 0.45
gamma <- 0.2
CM_Opt<- ModProc_CM(DDist,beta,gamma,ODEmaxTime = 100, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
MA_Opt<- MASIR_Proc(beta, gamma, lambda, init_S = (N-1)/N, ODEmaxTime=100, ODEstep=1e-1,TrackDyn = T)
Mod_Opt<- MAmod_Proc(beta, gamma, lambda, init_S = (N-1)/N, ODEmaxTime=100, ODEstep=1e-1,TrackDyn = T)

CM_Opt$R0
beta/(beta+gamma)*lambda
Mod_Opt$Rnet/gamma


#1+log((N-1)/N)/lambda

CM_out <- CM_Opt$Dynamic
MA_out <- MA_Opt$Dynamic
Mod_out <- Mod_Opt$Dynamic

time <- CM_out[,1]
CM_I <- CM_out[,5]
MA_I <- MA_out[,3]
Mod_I <- Mod_out[,3]

CM_R <- CM_out[,3]
MA_R <- MA_out[,4]
Mod_R <- Mod_out[,4]

theta <- CM_out[,2]

St<-CM_out[,4]

CM_S <- CM_out[,4]
MA_S <- MA_out[,2]
Mod_S <- Mod_out[,2]

dat_S <- cbind(time,CM_S, MA_S, Mod_S)
dat_R <- cbind(time,CM_R, MA_R, Mod_R)
dat <- cbind(time,CM_I, MA_I, Mod_I)

ggplot(data = dat)+theme_bw()+
  geom_line(aes(x=time, y=CM_I,color="Network"))+
  geom_point(aes(x=time, y=MA_I,color="MASIR"),alpha=0.1)+
  geom_point(aes(x=time, y=Mod_I,color="Modified"),alpha=0.1)+
  scale_color_manual(values=c("black", "red","blue"))+
  xlim(0,10)+
  labs(y = "I(t)") 

ggplot(data=dat_S)+theme_bw()+
  geom_line(aes(x=time, y=CM_S,color="Network"))+
  geom_line(aes(x=time, y=MA_S,color="MASIR"))+
  geom_point(aes(x=time, y=Mod_S,color="Modified"),alpha=0.1)+
  scale_color_manual(values=c("black", "red","blue"))+
  xlim(0,10)+
  labs(y = "S(t)") 

ggplot(data = dat_R)+theme_bw()+
  geom_line(aes(x=time, y=CM_R,color="Network"))+
  geom_line(aes(x=time, y=MA_R,color="MASIR"))+
  geom_point(aes(x=time, y=Mod_R,color="Modified"),alpha=0.1)+
  scale_color_manual(values=c("black", "red","blue"))+
  xlim(0,10)+
  labs(y = "R(t)") 


Mod_S
Mod_I
beta
gamma
lambda
theta <- (1-log(Mod_S)/lambda)

def_reff<- -lambda*CM_S*(-(beta+gamma)*(1+log(CM_S)/lambda)+beta*CM_S+gamma)/(CM_I*(gamma+beta))
cal_reff<- beta/(beta+gamma)*lambda*CM_S*(1+log(CM_S)/lambda)
dat_reff <- cbind(time,def_reff,cal_reff,Mod_I,Mod_S,theta)
ggplot(data=dat_reff)+theme_bw()+
  geom_line(aes(x=time, y=def_reff,color="Def"))+
  geom_line(aes(x=time, y=cal_reff,color="cal_eff"))+
  #geom_line(aes(x=time, y=disc_reff,color="disc_eff"))+
  geom_line(aes(x=time, y=theta*5,color="Theta"))+
  geom_line(aes(x=time, y=CM_I*5, color="I"))+
  geom_hline(yintercept=beta/(beta+gamma)*lambda,color="blue")+
  geom_hline(yintercept=1,color="green")+
  xlim(0,10)+
  #scale_color_manual(values=c("red", "black","brown"))
  labs(y = "R_eff") 


which.max(def_reff)
which.max(Mod_I)
max(which(cal_reff>1))
cal_reff[329]
cal_reff[330]
