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
	|> mutate(Ri = beta*VE/(N*I)/gamma)
	|> filter(I>0)
)

Infector <- (Infector
	|> filter(!is.na(RecoveryTime))
)

## summary(State)
## summary(Infector)

Cohort <- (Infector
	|> mutate(t = round(InfectTime))
	|> group_by(t)
	|> summarize(I = n()/N, value = mean(NumInfected))
	|> mutate(name="Rc")
)

summary(Cohort)

## TODO: get I values and use them as sizes
Obs <- (State
	|> select(t, Ri)
	|> pivot_longer(-t)
)

summary(Obs)

Rf <- bind_rows(Obs, Cohort)

rdsSave(Rf)
