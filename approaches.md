
## Resources 

[This notes document](notes_NovoANDNetwork.md) has the netSusc result; try to pull it out and contextualize a bit more.

## Draft

In early stages of spread, the reproductive number is given by $R_0$. This can be thought of as a product of the average infectiousness of infectors, and the average susceptibility of susceptible individuals. As the disease spreads through the population, the effective reproductive number, $R_\text{eff}$, declines typically declines from $R_0$. Setting aside Behavioral or policy changes, there are still three main reasons for this: 
- Decreasing infectiousness of the infected population
- Decreasing susceptibility of the susceptible population
- And decreasing proportions susceptible in the population.

We are interested in the rate of this decrease, and for this purpose, we like to think about the ratio $\sigma$, which we define as the ratio between $R_\text{eff}$, and the naive value of $R_\text{eff}$ given by $S\times R_0$. In a homogeneous, well mixed model $\sigma=1$, meaning that $R_\text{eff}$ is exactly proportional to $S$; in a heterogeneous model, we generally expect $\sigma<1$.

We are interested in sensible functional forms for $\sigma$ in a variety of scenarios, in particular with or without accounting for: 1) changes infection in the infected population, 2) the discrete nature of contacts, and 3) over-dispersion. For now, we should set aside the local nature of mixing, but it would be a fun thing to come back to. 

The Dwyer-Parsons approach accounts for “over-dispersion” (continuous heterogeneity). It can be derived either with a moment approximation, or via a gamma assumption, and is equivalent to $σ=S^{κ}$. The Novozhilov approach accounts for discrete contacts, but none of the other factors. 

Recently, Richard and Jonathan have come up with an approach that generalizes these two for a negative binomial degree-distribution, maybe call it the network-susceptibility approach. A straightforward argument gives:
$$σ=S^κ+\frac{S^κ-1}{κ \delta}, $$
where $\kappa$ is the Dushoff-style squared CV of an underlying continuous distribution (in this case, the Gamma distribution that underlies the negative binomial), so that $\kappa=0$ corresponds to the Poisson, with  $\sigma=1+\frac{log(S)}{\delta}$ (the Novozhilov result), $\delta\to\infty$ corresponds to the Dwyer result $\sigma=S^{\kappa}$, and $\kappa=1$ corresponds to the geometric case, where $\sigma=S+\frac{S-1}{\delta}$ (is this a known result?)

This approach still focuses only on changing susceptibility (not changes in average infectiousness of infectors); this should be the next thing to think about. We're also not addressing locality, which will be more difficult, but we should not shy away from that, because there is too much shying away going on.

We are particularly interested in mining more from the Miller approach, and seeing if we can come up with simple approaches that account for changing changes in both infectiousness and susceptibility. We anticipate some difficulties here in terms of accounting for recovery.

<!-- Dictation code

Ignore this.

:%s/Arnold\c/R0/g
:%s/Arthur\c/Reff/g
:%s/Sigma\c/σ/g
:%s/Nova\c/Novozhilov/g

-->

