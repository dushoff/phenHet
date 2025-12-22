library(shellpipes)
rpcall("viewMSV.Rout viewMSV.R simMSV.rda")

library(dplyr)

loadEnvironments()

dat <- as.data.frame(Rvs_out)

print(dat
	|> mutate(slope = (0.5-P)/(1-theta))
)

library(ggplot2)
p <- (ggplot(dat)
	+ aes(theta, P)
	+ geom_point()
)

print(p + dat |> filter(theta>0.99))

