library(shellpipes)
manageConflicts()
rpcall("scaleFuns.R edgelist.cpp big.params.R slow/big.post.rds NetworkODE.R")
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
S0Count <- PGFG0(1-it_omega,DDist)*N
Eigen_R <- EigenR(DDist, beta, gamma, lambda, init_omega = it_omega)

### MSV with eigen direction
CM_Opt<- ModProc_CM(  DDist,beta,gamma
                      , ODEmaxTime = 100
                      , ODEstep = 2e-2
                      , init_omega = it_omega
                      , TrackDyn = TRUE
                      , init_R = Eigen_R
)
theta_inf <- CM_Opt$ThetaInfinity
CM_out <- CM_Opt$Dynamic


#### Reverse ODE for p(t) and R_c
theta_inf
tps <- 5001
# print(CM_out[tps,])

Rvs_vec<-as.numeric(CM_out[tps,])

{ t_Rvs<-Rvs_vec[1]
  theta_Rvs <- Rvs_vec[2]
  R_Rvs <- Rvs_vec[3]
}

RVS_args <- list(  DDist
                   , beta
                   , gamma
                   , theta_Rvs
                   , R_Rvs
                   , ODEmaxTime = t_Rvs
                   , ODEstep = 2e-2)

Rvs_out <- do.call("Rvs_ODE", c(RVS_args))

Rvs_out[,1]<- t_Rvs-Rvs_out[,1]
Rvs_out <- Rvs_out[order(Rvs_out[,1]),]

### Two Data Frame for MSV

Rvs_df<-as.data.frame(Rvs_out)
CM_df <- as.data.frame(CM_out)

### R_c*
theta<-CM_df$theta
R_cstar <- beta/(beta+gamma)*PGFd2G0(theta,DDist)/lambda

### R_c
PLen<-length(Rvs_df$P)
R_c <- Rvs_df$P*(PGFd2G0(theta[1:PLen],DDist)/PGFd1G0(theta[1:PLen],DDist))
# Rvs_df$R_c <- R_c

### R_i
theta_dot <- -beta*theta+gamma*(1-theta)+beta/lambda*PGFd1G0(theta,DDist)
S_dot <- theta_dot*PGFd1G0(theta,DDist)
R_i <- -S_dot/(CM_df$I_out*gamma)

time <- Rvs_df$time
MSV<-as.data.frame(cbind(time, R_cstar, R_c, R_i))

library(ggplot2); theme_set(theme_bw())


loadEnvironments()
Rf <- rdsRead("slow/big.post.rds")

pal <- okabe_ito <- c("#E69F00", "#56B4E9", "#009E73") 

cplot <- (ggplot(Rf)
	+ aes(t+1, value, color=name, linetype="Sim")
	+ geom_line()
	+ geom_point(aes(size=obs),alpha=0.2)
	+ geom_line(data=MSV,aes(time,R_c,color = "Rc", linetype="MSV"))
	+ geom_line(data=MSV,aes(time,R_cstar,color = "Rcstar", linetype="MSV"))
	+ geom_line(data=MSV,aes(time,R_i,color = "Ri", linetype="MSV"))
	+ coord_cartesian(xlim=c(0, 80))
	+ scale_color_manual(values = pal)
	+ scale_size_area()
	+ geom_hline(yintercept=Ri0, color=pal[[3]])
	+ geom_hline(yintercept=Rc0, color=pal[[2]])
)
print(cplot)
# print(cplot + scale_y_log10())

