
## Resources 

notes_NovoANDNetwork.md has the netSusc result; try to pull it out and contextualize a bit more.

## Draft

In early stages of spread, the reproductive number is given by R0. This can be thought of as a product of the average infectiousness of infectors, and the average susceptibility of susceptible individuals. As the disease spreads through the population, the effective reproductive number, Reff, declines typically declines from R0. Setting aside Behavioral or policy changes, there are still three main reasons for this decreasing infectiousness of the infected population, decreasing susceptibility of the susceptible population, and decreasing proportions susceptible in the population.

We are interested in the rate of this decrease, and for this purpose, we like to think about the ratio σ, which is defined as the ratio between Reff, and a naive value of Reff given by S Reff. In a homogeneous, well mixed model σ is one.

We are interested in functional forms sensible functional forms for σ in a variety of scenarios, in particular with or without accounting for: changes infection in the infected population, the discrete nature of contacts, and over-dispersion. For now, we should set aside the local nature of mixing, but it would be a fun thing to come back to. 

The Dwyer-Parsons approach accounts for “over-dispersion” (continuous heterogeneity). It can be derived either with a moment approximation, or via a gamma assumption, and is equivalent to σ=S^{-κ}. The Novozhilov approach accounts for discrete contacts, but none of the other factors. Recently, Richard and Jonathan have come up with an approach that generalizes these two, maybe call it the network-susceptibility approach. So it is still focused only on changing susceptibility, and it's not addressing locality, which is very much in our future. 

We are particularly interested in mining more from the Miller approach, and seeing if we can come up with simple approaches that account for changing changes in both infectiousness and susceptibility. There are some dangers there are some difficulties here in terms of accounting for recovery.

## What?

:%s/Arnold\c/R0/g
:%s/Arthur\c/Reff/g
:%s/Sigma\c/σ/g
:%s/Nova\c/Novozhilov/g

