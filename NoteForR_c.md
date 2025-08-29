As mentioned in [JR_Negbinom_Result](JR_Negbinom_Result.md), the Zhao 2 result leads to 
$$
\begin{align}
\mathcal{R}_c=\mathbb{E}[X_t]&=\mathbb{E}_{K_I^*}[\mathbb{E}[X_t|K_I^*]]
\\
& =\mathbb{E}_{K_I^*}[\mu(K_I^*-1)]
\\
& =\mu \mathbb{E}[K_I^*-1]
\\
& = \mu (\mathbb{E}[K_I^*]-1)
\\
& = \mu \phi\frac{G''_p(\phi)}{G'_p(\phi)}
\end{align}
$$

$\mu$ should be the probability that a random neighbor of the newly infected focal vertex is eventually infected by the focal vertex.
- For each stub(half-edge) connected to the newly infected focal vertex, except the one connected to its infector (thus the edge belongs to $K_I^*-1$).
	- Consider a random paring process, we know one side of the stubs connects to the newly infected focal vertex, so such edge it forms has not transmitted the infection. So the pairing stub is not chosen from all stubs in the graph, but only from those stubs belongs to $\phi=\phi_S+\phi_I+\phi_R$.
	- Since we want to count the number of infection caused by the focal vertex, the numerator of the probability is $\phi_S$, which corresponding to the proportion of susceptible neighbors.
	- The $\frac{\phi_S}{\phi}$ term is the probability that such stub connects to a susceptible vertex. 
- Denote $T$ is the exponential random variable for recover time of the focal (or a random) infected vertex with recovery rate $\gamma$, the transmission rate is $\beta$. 
	- **If we assume that every susceptible neighbor is only infected by the focal vertex**, then the probability a susceptible susceptible is infected by the focal vertex before its recovery is $$1-e^{-\beta T}$$
	- The expectation of such probability is the transmissibility $\tau$ $$\tau=\mathbb{E}_{T}[1-e^{\beta T}]=\int_{0}^{\infty}(1-e^{-\beta T})\times \gamma e^{-\gamma T}dT=\frac{\beta}{\beta+\gamma}$$
	- With such assumption, we have the probability $\mu^*$, such that $$\mu^*=\tau \times \frac{\phi_S}{\phi}=\frac{\beta}{\beta+\gamma}\times \frac{\phi_S}{\phi}=\frac{\beta}{\beta+\gamma}\times\frac{G'_p(\phi)}{\phi G'_p(1)}$$, we have:$$\mathcal{R}^*_c= \mu (\mathbb{E}[K_I^*]-1)=\frac{\beta}{\beta+\gamma}\times \frac{\phi_S}{\phi} \times (\mathbb{E}[K_I^*]-1)=\frac{\beta}{\beta+\gamma}\frac{G''_p(\phi)}{\delta}$$
However, as Todd mentioned, the assumption is only true at the beginning when there is no other infected vertex to compete with the focal vertex. 

For the true value of $\mathcal{R}_c$ and $\mu$ it should be a proportion of the $\mathcal{R}^*_c$ and $\mu^*$ based on time $t$, corresponding to the proportion of susceptible neighbors of the newly infected focal vertex, but only infected by it during its period of infection $T$.

Consider the definition of $\phi_S$, which is defined as:
$$\phi_S(t)=\sum_k (\frac{k\times p_k}{\delta}\times\phi(t)^{k-1})=\frac{G'_p(\phi(t))}{\delta}$$
- The $\frac{k\times p_k}{\delta}$ term is the probability of the neighbor of focal vertex have degree $k$.
- $\phi(t)^{k-1}$ is the probability that the neighbor is not infected through other edges at time $t$, except by the one connected to the focal vertex.
Similarly, the newly infected focal vertex with recovery time $T$ at time $t$, the probability of its neighbor is not infected by any other vertices at time $t+T$ is:
$$\phi_S(t+T)=\sum_k (\frac{k\times p_k}{\delta}\times\phi(t+T)^{k-1})=\frac{G'_p(\phi(t+T))}{\delta}$$
Then we need to consider the exponential distributed recovery time $T$, then the expectation among all possible $T$ value would be:
$$\int_0^{+\infty}\sum_k [\frac{k\times p_k}{\delta}\times\phi(t+T)^{k-1}]\gamma e^{-\gamma T}dT=\frac{1}{\delta}\int_0^{+\infty}G'_p(\phi(t+T)) \gamma e^{-\gamma T} dT$$

???? Conditional Probability at time $t$ ?????