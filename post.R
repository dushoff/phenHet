library(shellpipes)
manageConflicts()

library(dplyr)
library(tidyr)

loadEnvironments()
result <- rdsRead()

attach(result)
stopifnot(FinalStat$Isize == 0)

## This code used to calculate Rcstar -- but not correctly
## Need to know about _new_ infectors, which means going back to the sim.
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
	|> mutate(t = round(InfectTime))
	|> group_by(t)
	|> summarize(obs = n()/N, Rc = mean(NumInfected)
		, Rcstar = rho*mean(Targets)/(rho+1)
	)
	|> pivot_longer(c(Rc, Rcstar))
)

summary(Cohort)

## TODO: get I values and use them as sizes
Obs <- (State
	|> select(t, obs=I, Ri)
	|> pivot_longer(Ri)
)

summary(Obs)

Rf <- bind_rows(Obs, Cohort)

rdsSave(Rf)
