library(shellpipes)
rpcall("plotsMSV.Rout plotsMSV.R scaleFuns.rda big.params.R slow/big.post.rds NetworkODE.R")
rpcall("simMSV.Rout simMSV.R scaleFuns.rda big.params.rda NetworkODE.R")
manageConflicts()
rpcall("scaleFuns.R edgelist.cpp big.params.R slow/big.post.rds NetworkODE.R")

loadEnvironments()
sourceFiles()

gamma <- 1/D
beta <- rho*gamma
r <- if(kappa==0) Inf else(1/kappa)
omega <- lambda*(1+kappa)+1
Rc0 <- rho*(omega-1)/(rho+1)
Ri0 <- rho*(omega-2)

kvalue <- seq(0,400)
Pk<- dnbinom(kvalue,r,mu=lambda)
DDist <- data.frame(kvalue,Pk)

# source("NetworkODE.R")
it_omega <- Init_omega_func(1,N)
Eigen_R <- EigenR(DDist, beta, gamma, lambda, init_omega = it_omega)

### Run MSV forward starting from eigen-direction
CM_Opt<- ModProc_CM(  DDist,beta,gamma
                      , ODEmaxTime = 100
                      , ODEstep = 1e-1
                      , init_omega = it_omega
                      , TrackDyn = TRUE
                      , init_R = Eigen_R
)
theta_inf <- CM_Opt$ThetaInfinity
CM_out <- CM_Opt$Dynamic

#### Reverse ODE for p(t) and R_c
theta_inf
tps <- length(CM_out[,1])
# print(CM_out[tps,])

## Pull out endpoint of previous simulation
Rvs_vec<-as.numeric(CM_out[tps,])
t_Rvs<-Rvs_vec[1]
theta_Rvs <- Rvs_vec[2]
R_Rvs <- Rvs_vec[3]

RVS_args <- list(  DDist
                   , beta
                   , gamma
                   , theta_Rvs
                   , R_Rvs
                   , ODEmaxTime = t_Rvs
                   , ODEstep = 1e-1)

Rvs_out <- do.call("Rvs_ODE", c(RVS_args))

Rvs_out[,1]<- t_Rvs-Rvs_out[,1]
Rvs_out <- Rvs_out[order(Rvs_out[,1]),]

### Two Data Frame for MSV
Rvs_df<-as.data.frame(Rvs_out)
CM_df <- as.data.frame(CM_out)

### R_c*
theta<-CM_df$theta
S <- CM_df$S_out

### GF version
# R_cstar <- beta/(beta+gamma)*PGFd2G0(theta,DDist)/lambda

### Analytically +S
R_cstar <- beta/(beta+gamma)*(omega-1)*S^(2*kappa+1)

### R_c
PLen<-length(Rvs_df$P)
### GF version
# R_c <- Rvs_df$P*(PGFd2G0(theta[1:PLen],DDist)/PGFd1G0(theta[1:PLen],DDist))

### Analytically +S
R_c <- Rvs_df$P*(lambda*(omega-1)*S^(2*kappa+1)/(S*lambda/(1+kappa*lambda-theta*kappa*lambda)))
# Rvs_df$R_c <- R_c

### R_i
### GF version
# theta_dot <- -beta*theta+gamma*(1-theta)+beta/lambda*PGFd1G0(theta,DDist)
# S_dot <- theta_dot*PGFd1G0(theta,DDist)

### Analytically +S
theta_dot <- -beta*theta+gamma*(1-theta)+beta/lambda*(S*lambda/(1+kappa*lambda-theta*kappa*lambda))
S_dot <- theta_dot*(S*lambda/(1+kappa*lambda-theta*kappa*lambda))

R_i <- -S_dot/(CM_df$I_out*gamma)

time <- Rvs_df$time
MSV<-as.data.frame(cbind(time, R_cstar, R_c, R_i))
