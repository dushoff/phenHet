library(shellpipes)
manageConflicts()

library(dplyr)
library(tidyr)

loadEnvironments()
result <- rdsRead()

attach(result)
stopifnot(FinalStat$Isize == 0)

State <- (State
	|> mutate(Ri = rho*VE/I)
	|> filter(I>0)
)
summary(State)

Infector <- (Infector
	|> filter(!is.na(RecoveryTime))
)
summary(Infector)

Cohort <- (Infector
	|> mutate(c = round(InfectTime))
	|> group_by(c)
	|> summarize(obs = n()/N, Rc = mean(NumInfected)
		, Rcstar = rho*mean(Targets)/(rho+1)
	)
	|> pivot_longer(c(Rc, Rcstar))
	|> mutate(t=c/D)
)

summary(Cohort)

Obs <- (State
	|> transmute(t=t/D, obs=I, Ri)
	|> pivot_longer(Ri)
)

summary(Obs)

Rf <- bind_rows(Obs, Cohort)

rdsSave(Rf)
