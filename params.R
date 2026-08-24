library(shellpipes)

sourceFiles()

gamma <- 1/D
beta <- rho*gamma
r <- if(kappa==0) Inf else(1/kappa)
omega <- lambda*(1+kappa)+1
Rc0 <- rho*(omega-1)/(rho+1)
Ri0 <- rho*(omega-2)

saveEnvironment()
