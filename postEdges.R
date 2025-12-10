library(shellpipes)
manageConflicts()

library(ggplot2); theme_set(theme_bw())
library(dplyr)

loadEnvironments()

names(result)
attach(result)

summary(State)
summary(Infector)

State <- (State
	|> mutate(Ri = VE/(N*I)/gamma)
)

print(ggplot(State)
	+ aes(t, VE)
	+ geom_line()
	+ geom_smooth()
)

cplot <- (ggplot(State)
	+ aes(t, Ri)
	+ geom_line()
	+ geom_smooth(color="black")
	+ geom_smooth(data=Infector
		, aes(InfectTime, NumInfected)
	)
	+ coord_cartesian(xlim=c(0, 50))
)
print(cplot)
print(cplot + scale_y_log10())
