library(shellpipes)
manageConflicts()

library(ggplot2); theme_set(theme_bw())
library(dplyr)

loadEnvironments()

attach(result)

print(FinalStat)

State <- (State
	|> mutate(
		, Ri = beta*VE/(N*I)/gamma
		, Rcstar = Ri*gamma/(beta+gamma)
	)
)

Infector <- (Infector
	|> filter(!is.na(InfectTime))
)

summary(State)
summary(Infector)

print(ggplot(State)
	+ aes(t, VE)
	+ geom_point()
	+ geom_smooth()
)

cplot <- (ggplot(State)
	+ aes(t, Ri)
	+ geom_line()
	+ geom_smooth(color="black")
	+ geom_line(aes(y=Rcstar), color="red")
	+ geom_smooth(aes(y=Rcstar), color="red")
	+ geom_smooth(data=Infector
		, aes(InfectTime, NumInfected)
	)
	+ coord_cartesian(xlim=c(0, 50))
)
print(cplot)
print(cplot + scale_y_log10())
