library(shellpipes); manageConflicts()
library(ggplot2); theme_set(theme_bw())
library(dplyr)

loadEnvironments()

ssteps <- 100
d <- (expand.grid(
	x=(1:ssteps)/ssteps
	, delta=c(5, 20, Inf)
	, kappa=c(0, 0.5, 1, 2)
))

d <- (d
	|> mutate (sig=zhao1(x, delta, kappa))
)

print(ggplot(d |> mutate(delta = as.factor(delta)))
	+ aes(x, sig, color=delta)
	+ geom_line()
	+ facet_wrap(~kappa, labeller = label_both)

)

