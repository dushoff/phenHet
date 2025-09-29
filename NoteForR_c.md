# 1. Reproductive number by case $\mathcal{R}^*_c$ with no competing infection
As mentioned in [JR_Negbinom_Result](JR_Negbinom_Result.md), the Zhao 2 result leads to 
$$\begin{align}
\mathcal{R}_c=\mathbb{E}[X_t]&=\mathbb{E}_{K_I^*}[\mathbb{E}[X_t|K_I^*]]
\\
& =\mathbb{E}_{K_I^*}[\mu(K_I^*-1)]
\\
& =\mu \mathbb{E}[K_I^*-1]
\\
& = \mu (\mathbb{E}[K_I^*]-1)
\\
& = \mu \phi\frac{G''_p(\phi)}{G'_p(\phi)}
\end{align}$$

$\mu$ should be the probability that a random neighbor of the newly infected focal vertex is eventually infected by the focal vertex.
- For each stub(half-edge) connected to the newly infected focal vertex, except the one connected to its infector (thus the edge belongs to $K_I^*-1$).
	- Consider a random paring process, we know one side of the stubs connects to the newly infected focal vertex, so such edge it forms has not transmitted the infection. So the pairing stub is not chosen from all stubs in the graph, but only from those stubs belongs to $\phi=\phi_S+\phi_I+\phi_R$.
	- Since we want to count the number of infection caused by the focal vertex, the numerator of the probability is $\phi_S$, which corresponding to the proportion of susceptible neighbors.
	- The $\frac{\phi_S}{\phi}$ term is the probability that such stub connects to a susceptible vertex. 
- Denote $T$ is the exponential random variable for recover time of the focal (or a random) infected vertex with recovery rate $\gamma$, the transmission rate is $\beta$. 
	- **If we assume that every susceptible neighbor is only infected by the focal vertex**, then the probability a susceptible susceptible is infected by the focal vertex before its recovery is $$1-e^{-\beta T}$$
	- The expectation of such probability is the transmissibility $\tau$ $$\tau=\mathbb{E}_{T}[1-e^{\beta T}]=\int_{0}^{\infty}(1-e^{-\beta T})\times \gamma e^{-\gamma T}dT=\frac{\beta}{\beta+\gamma}$$
	- With such assumption, we have the probability $\mu^*$, such that $$\mu^*=\tau \times \frac{\phi_S}{\phi}=\frac{\beta}{\beta+\gamma}\times \frac{\phi_S}{\phi}=\frac{\beta}{\beta+\gamma}\times\frac{G'_p(\phi)}{\phi G'_p(1)}$$, we have:$$\mathcal{R}^*_c= \mu (\mathbb{E}[K_I^*]-1)=\frac{\beta}{\beta+\gamma}\times \frac{\phi_S}{\phi} \times (\mathbb{E}[K_I^*]-1)=\frac{\beta}{\beta+\gamma}\frac{G''_p(\phi)}{\delta}$$


# 2. Reproductive number by case $\mathcal{R}_c$ with competing infection
However, as Todd mentioned, the assumption is only true at the beginning when there is no other infected vertex to compete with the focal vertex. 

For the true value of $\mathcal{R}_c$ and $\mu$ it should be a proportion of the $\mathcal{R}^*_c$ and $\mu^*$ based on time $t$, corresponding to the proportion of susceptible neighbors of the newly infected focal vertex, but only infected by it during its period of infection $T$.

Consider the definition of $\phi_S$, which is defined as:
$$\phi_S(t)=\sum_k (\frac{k\times p_k}{\delta}\times\phi(t)^{k-1})=\frac{G'_p(\phi(t))}{\delta}$$
- The $\frac{k\times p_k}{\delta}$ term is the probability of the neighbor of focal vertex have degree $k$.
- $\phi(t)^{k-1}$ is the probability that the neighbor is not infected through other edges at time $t$, except by the one connected to the focal vertex.
Similarly, the newly infected focal vertex with recovery time $T$ at time $t$, the probability of its neighbor is not infected by any other vertices at time $t+T$ is:
$$\phi_S(t+T)=\sum_k (\frac{k\times p_k}{\delta}\times\phi(t+T)^{k-1})=\frac{G'_p(\phi(t+T))}{\delta}$$
Then we need to consider the exponential distributed recovery time $T$, then the expectation among all possible $T$ value would be:
$$\hat{\phi}_S=\int_0^{+\infty}\sum_k [\frac{k\times p_k}{\delta}\times\phi(t+T)^{k-1}]\gamma e^{-\gamma T}dT=\frac{1}{\delta}\int_0^{+\infty}G'_p(\phi(t+T)) \gamma e^{-\gamma T} dT$$
Accordingly, we have the "corrected" reproductive number by case:
$$
\mathcal{R}_c= \mathcal{R}^*_c\times\frac{\hat{\phi}_S}{\phi_S}=\frac{\beta}{\beta+\gamma}\times \frac{\hat{\phi}_S}{\phi} \times (\mathbb{E}[K_I^*]-1)
$$
The problem now comes to how to find $\hat{\phi}_S$ with such integration.

Some thoughts for simplifying the integration:
$$\begin{align}
& \int_0^{+\infty}\phi(t+T) \gamma e^{-\gamma T} dT 
\\
= &[\phi(t+T)(-e^{-\gamma T})]
\bigg| ^{+\infty}_{0}- \int_0^{+\infty}\dot{\phi}(t+T) (-e^{-\gamma T}) dT
\\
=&\phi(t)-\int_0^{+\infty}\dot{\phi}(t+T) (-e^{-\gamma T}) dT
\\
=&\phi(t)+\int_0^{+\infty}[-\beta\phi(t+T)+\beta\frac{G'_p(\phi(t+T))}{\delta}+\gamma(1-\phi(t+T))] (e^{-\gamma T}) dT
\\
=&\phi(t)-\frac{\beta+\gamma}{\gamma}\int_0^{+\infty}\phi(t+T) \gamma e^{-\gamma T} dT+\int_0^{+\infty}\gamma e^{-\gamma T} dT+\frac{1}{\gamma}\int_0^{+\infty}G'_p(\phi(t+T)) \gamma e^{-\gamma T} dT
\\
=&\phi(t)+\frac{1}{\gamma}-\frac{\beta+\gamma}{\gamma}\int_0^{+\infty}\phi(t+T) \gamma e^{-\gamma T} dT+\frac{1}{\gamma}\int_0^{+\infty}G'_p(\phi(t+T)) \gamma e^{-\gamma T} dT
\end{align}$$
Rearrange the term give us:
$$
\int_0^{+\infty}G'_p(\phi(t+T)) \gamma e^{-\gamma T}dT=-\gamma\phi(t)-1+(\beta+2 \gamma)\int_0^{+\infty}\phi(t+T) \gamma e^{-\gamma T} dT
$$

If ITB with $G'_p(\phi)$ directly, we have:
$$\begin{align}
& \int_0^{+\infty}G'_p(\phi(t+T)) \gamma e^{-\gamma T}dT
\\
=& [G'_p(\phi(t+T))(-e^{-\gamma T})]
\bigg| ^{+\infty}_{0}- \int_0^{+\infty}G''_p(\phi(t+T))\dot{\phi}(t+T) (-e^{-\gamma T}) dT
\\
=& G'_p(\phi(t))-\int_0^{+\infty}G''_p(\phi(t+T))\dot{\phi}(t+T) (-e^{-\gamma T}) dT
\end{align}$$
Might get something if we take in the NegBinom distribution????

## 2.1 ODE Correction Term

This is a translated version based on [Todd's notes](outputs/Rc.pdf).

This seems to (and should) be equivalent to previous integration results for $\hat{\phi}_S$ (To Do: verify). The improvement is unlike the troublesome integration, we are now able to describe the result with a short ODE and being able to solve it numerically together with ODE system of MSV framework.

Here is the derivation:

We would like to compute:
$$p(t)=\mathbb{P}\{\text{the focal infected node infected at time } t\text{ infect random one of its neighbour}\}$$
- As in MSV framework, we assume infection of neighbors are independent, thus the number of infected neighbors follows binomial distribution. Follow the idea of Zhao2 result/derivation:
	- Correspond to $\mu^*$ we have: $$\mu(t)=p(t)\times\frac{\phi_S(t)}{\phi(t)}$$
	- Furthermore, Correspond to $\mathcal{R}^*_{c}$, we have $$\mathcal{R}_c(t)=\mu(t) \times (\mathbb{E}[K_I^*]-1)=\frac{p(t)}{\tau}\times\tau\frac{\phi_S}{\phi}\times(\mathbb{E}[K_I^*]-1)=\frac{p(t)}{\tau}\times\mathcal{R}^*_{c}(t)$$
	- This further indicate $$p(t)=\tau\times\frac{\hat{\phi}_S}{\phi_S}$$if the two results are equivalent (to be verified).
- Initially at $t=0$, there should be no competing infection among neighbors, all neighbors of focal infected node can only be infected by the focal node.
	- This leads to $$p(0)=\tau=\frac{\beta}{\beta+\gamma}$$ as initial value of $p(t)$
	- Furthermore, this agree with the initial value for $\mathcal{R}_c$, s.t.$$\mathcal{R}_{c}(0)=1\times\mathcal{R}^*_{c}(0)=\mathcal{R}_{c,0}=\frac{\beta}{\beta+\gamma}\times\frac{G''_p(1)}{\delta}$$

To find $p(t)$ in the stochastic process, we observe that the focal infected node will infect its neighbor:
- if the focal node makes an infectious contact with the neighbor under rate $\beta$ before the focal node's recovery under rate $\gamma$
- $\textbf{AND}$ the neighbor is not infected (by the neighbor's neighbor other than the focal) at the time of contact by focal.
With the MSV framework, the probability that a randomly chosen neighbor of the focal node is not infected at time $t$ is $$\phi_S(t)=G_q(\phi(t))=\frac{G'_p(\phi(t))}{\delta}$$
In the random events, lets set the following random variables for time:
- $T_r$: the time after infection $t$ that the focal infected node recovers. Based on the recovery rate $\gamma$ and exponential distribution, we have$$\mathbb{P}(T_r>s)=e^{-\gamma s}$$
- $T_c$: the time after infection $t$ that the focal node makes it infectious contact with the neighbor through the edge connecting them. Based on the infection rate $\beta$ and exponential distribution, we have $$\mathbb{P}(T_c>s)=e^{-\beta s}$$
- $T_n$: the first time that the neighbor has an infectious contact from one of its other neighbors than the focal node. We further have $$\mathbb{P}(T_n>t+s)=\phi_S(t+s)$$





## 2.2 Counterargument
A counterargument that RZ and BB have for this:
- For large enough random network in MSV's configuration framework, it has been proved loops in network exist, but will be really rare (a.s. no loop as $N$ increase). 
- Like in the percolation theory, considering the transmissibility of edges is a probability smaller than 1, loops linked by "occupied" edges are even more rare.
As a result, competing infection can only happens if we have such "occupied" loops, so the impact of competing infection should be very limited, even if the infection proportion is large.
This counter argument only valid for configuration network. If we consider small world network, then this impact would be much significant.
Ideas and comments are welcome!

