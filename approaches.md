## Resources 

[This notes document](notes_NovoANDNetwork.md) has the netSusceptibility result; try to pull it out and contextualize a bit more.

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
where $\kappa$ is the Dushoff-style squared CV of an underlying continuous distribution (in this case, the Gamma distribution that underlies the negative binomial), so that $\kappa=0$ corresponds to the Poisson, with  $\sigma=1+\frac{log(S)}{\delta}$ (the Novozhilov Poisson result), $\delta\to\infty$ corresponds to the Gamma, with $\sigma=S^{\kappa}$(the Dwyer result), and $\kappa=1$ corresponds to the geometric case, where $\sigma=S+\frac{S-1}{\delta}$ (an equivalent result has been derived by [RomanescuEtAL(2023)](https://doi.org/10.1016/j.epidem.2023.100708) and [Novozhilov(2008)](./refs/Novozhilov2008.pdf), with different settings.)

This approach still focuses only on changing susceptibility (not changes in average infectiousness of infectors); this should be the next thing to think about. We're also not addressing locality, which will be more difficult, but we should not shy away from that, because there is too much shying away going on.

The result is elegant and match well with the past, but we need more clear descriptions for this susceptibility-only network framework (Conjecture: in-degree model)

We are particularly interested in mining more from the MSV configuration approach, and seeing if we can come up with simple approaches that account for changes in both infectiousness and susceptibility. We anticipate some difficulties here in terms of accounting for recovery.

Romanescu claim that their result for effective reproduction number could be applied on specific degree distributions in random-network framework, but have some logic and definition problems and provide little validation to simulation. As mentioned earlier, we find their result agree with the susceptibility only heterogeneity results (ours with network, Novozhilov without network)

Based on their result, we further applied corrections to the derivation and generate the "ideal" case reproductive number $\mathcal{R}^*_c(t)$: the (counterfactual) expected number of infections caused by a new case during its infectious duration, if it is infected at time $t$ and not consider competing infections of its susceptible neighbours (all susceptible neighbors of the focal case can only infected by the focal case). We generate an nice expression for $\mathcal{R}^*_c(t)$.

We further considered the probability of competing infection and derived an extra ODE for $p(t)$, the "real" infection probability susceptible neighbour of a case infected at time $t$. As a final value problem, $p(t)$ could be solved reversely from the MSV ODE system. We are trying to find a good estimation for initial value $p(0)$ but not yet successful.

Using $p(t)$ to modify $\mathcal{R}^*_c(t)$ expression give us the "real" case reproductive number $\mathcal{R}_c(t)$ that considering the competing infection. $\mathcal{R}_c(t)$ is no lager than $\mathcal{R}^*_c(t)$ and converge to $\mathcal{R}^*_c(t)$ eventually as outbreak ends. If the network size $N$ is large enough such that loops in the network are rare enough, we expecting $\mathcal{R}_c(0)$ is close to $\mathcal{R}^*_c(0)$. But we are not able to find a good estimation for the difference.

On the other hand, we could define the instantaneous effective reproductive number $\mathcal{R}_i$ as the average number of new infection caused by a randomly chosen infectious case at time $t$, i.e.$$\mathcal{R}_i(t)=-\frac{\dot{S}}{\gamma I}$$
Unlike in the homogeneous model where $\mathcal{R}_i(t)$ are monotonically decreasing, in MSV network with a purely naive population, $\mathcal{R}_i(t)$ might increase for a while to reach the transmission's "eigen-status" and then monotonically decrease to 0 until the outbreak ends. We have find expression for the peak of $\mathcal{R}_i(t)$ and use the eigenvector idea find out the initial-status s.t. $\mathcal{R}_i(0)$ agree with the peak.

The problem of $\mathcal{R}_i(t)$ in MSV is that we does not have explicit expression for $I(t)$, so it is hard to derive the $\mathcal{R}_i(t)$ curve without solve the whole ODE system for $I(t)$. An ODE of $\dot{\mathcal{R}}_i(t)$ that not dependent on $I(t)$ is generated, but we not yet find useful understanding yet. We are also thinking about connection between $\mathcal{R}_i(t)$ and $\mathcal{R}_c(t)$

We created vertex based and a much faster (and new?) edge-based Gillespie algorithm on network transmission. The simulation results match our expression in number of runs, but we might need to verify this in a more systematic way.


<!-- Dictation code

Ignore this.

:%s/Arnold\c/R0/g
:%s/Arthur\c/Reff/g
:%s/Sigma\c/σ/g
:%s/Nova\c/Novozhilov/g

-->

