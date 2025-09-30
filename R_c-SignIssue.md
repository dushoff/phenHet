This is a translated version based on [Todd's notes](outputs/Rc.pdf).
We are able to describe the result with a short ODE and being able to solve it numerically together with ODE system of MSV framework.

We would like to compute:
$$p(t)=\mathbb{P}\{\text{the focal infected node infected at time } t\text{ infect random one of its neighbour}\}$$
- As in MSV framework, we assume infection of neighbors are independent, thus the number of infected neighbors follows binomial distribution. Follow the idea of Zhao2 result/derivation:
	- We count the condition that the neighbor must be susceptible $\phi_S(t)$ in $p(t)$ s.t. $$\mu(t)=p(t)$$
	- Furthermore, correspond to $\mathcal{R}^*_{c}$, we have $$\mathcal{R}_c(t)=\mu(t) \times (\mathbb{E}[K_I^*]-1)=p(t)\times\phi(t)\frac{G''_p(\phi(t))}{G'_p(\phi(t))}$$
- Initially at $t=0$, there should be no competing infection among neighbors, all neighbors of focal infected node can only be infected by the focal node.
	- This leads to $$p(0)=\tau=\frac{\beta}{\beta+\gamma}$$ as initial value of $p(t)$.
	- Furthermore, this agree with the initial value for $\mathcal{R}_c$, s.t.$$\mathcal{R}_{c}(0)=\mathcal{R}^*_{c}(0)=\mathcal{R}_{c,0}=\frac{\beta}{\beta+\gamma}\times\frac{G''_p(1)}{\delta}$$

To find $p(t)$ in the stochastic process, we observe that the focal infected node will infect its neighbor:
- if the focal node makes an infectious contact with the neighbor under rate $\beta$ before the focal node's recovery under rate $\gamma$
- $\textbf{AND}$ the neighbor is not infected (by the neighbor's neighbor other than the focal) at the time of contact by focal.
With the MSV framework, the probability that a randomly chosen neighbor of the focal node is not infected at time $t$ is $$\phi_S(t)=G_q(\phi(t))=\frac{G'_p(\phi(t))}{\delta}$$
In the random events, lets set the following random variables for time:
- $T_r$: the time after infection $t$ that the focal infected node recovers. Based on the recovery rate $\gamma$ and exponential distribution, we have$$\mathbb{P}(T_r>s)=e^{-\gamma s}$$
- $T_c$: the time after infection $t$ that the focal node makes it infectious contact with the neighbor through the edge connecting them. Based on the infection rate $\beta$ and exponential distribution, we have $$\mathbb{P}(T_c>s)=e^{-\beta s}$$
- $T_n$: the first time that the neighbor has an infectious contact from one of its other neighbors than the focal node. We further have $$\mathbb{P}(T_n>t+s)=\phi_S(t+s)$$
With these RVs, we can interpret $p(t)$ in the following probability: $$p(t)=\mathbb{P}\{ t+T_c < (t+T_r) \wedge T_n \}$$ where $t+T_c < (t+T_r) \wedge T_n \ \Leftrightarrow \min((t+T_r),T_n)$.

Since $T_r$ and $T_n$ are independent, we can further derive: 
$$\begin{align}
\mathbb{P}\{ (t+T_r) \wedge T_n > t+u \} & = \mathbb{P}\{ T_r>u\} \mathbb{P}\{T_n>t+u \}
\\
& = e^{-\gamma u} \phi_S(t+u)
\end{align}$$
Then we can rewrite the expression for $p(t)$ based on law of total probability:
$$\begin{align}
p(t) &=\mathbb{P}\{ t+T_c < (t+T_r) \wedge T_n \}
\\
&=\int_0^{\infty}\mathbb{P}\{t+u<(t+T_r) \wedge T_n | t+T_c =t+u\}\times\mathbb{P}\{ t+T_c =t+u\}du
\\
&=\int_0^{\infty}\mathbb{P}\{t+u<(t+T_r) \wedge T_n | T_c =u\}\times\mathbb{P}\{ T_c =u\}du
\\
&=\int_0^{\infty}\mathbb{P}\{t+u<(t+T_r) \wedge T_n\}\times\mathbb{P}\{ T_c =u\}du \quad \text{as }T_n, T_c, T_r \text{ are independent}
\\
&=\int_0^{\infty}[e^{-\gamma u} \phi_S(t+u)]\times [\beta e^{-\beta u}]du
\\
&=\int_0^{\infty}\beta e^{-(\beta+\gamma)u}\phi_S(t+u)du
\end{align}$$
and we have
$$\begin{align}
\frac{d}{dt}p(t) &=\frac{d}{dt}\int_0^{\infty}\beta e^{-(\beta+\gamma)u}\phi_S(t+u)du
\\
&=\int_0^{\infty}\beta e^{-(\beta+\gamma)u} \times [\frac{d}{dt}\phi_S(t+u)] du
\\
&= \int_0^{\infty}\beta e^{-(\beta+\gamma)u} \times [\dot{\phi}_S(t+u)] du
\\
&= \left.[\beta e^{-(\beta+\gamma)u}\phi_S(t+u)]\right |^{\infty}_{u=0} -\int_0^{\infty}-(\beta+\gamma) \beta e^{-(\beta+\gamma)u} \phi_S(t+u) du  \quad \text{(I.B.P.)}
\\
&= -\beta \phi_S(t)+(\beta+\gamma)\int_0^{\infty}\beta e^{-(\beta+\gamma)u} \phi_S(t+u) du
\\
&=-\beta \phi_S(t)+(\beta+\gamma)p(t)
\end{align}$$
as a ODE of $p(t)$.

At $t=0$, we expect to have $\phi_S(0)=1$ and $p(0)=\frac{\beta}{\beta+\gamma}$, so take these initial value into the ODE gives us $$\frac{d}{dt}p(t)=-\beta+(\beta+\gamma)\times\frac{\beta}{\beta+\gamma}=0$$ which agree with our expectation.

### Sign Problem
This notes fixed some typo in [Todd's original notes](outputs/Rc.pdf) and the final ODE agree with the original one. I review the derivation myself twice and it seems correct.

But we expect $p(t)$ be a probability and monotonically decreasing from its initial value since the following two factor:
- $\phi_S$ is decreasing as the infection spread out.
- We have more competing infection happens, i.e. it is more likely to have $T_n<t+T_c$ 
However, in numeric solving, this ODE for $p(t)$ leads to increasing $p(t)$ and getting larger than 1. 

See [R_c-SignIssue.R](R_c-SignIssue.R) for numeric solutions `CM_P`, the ODE is implemented at line `#107`.

Just an observation: if alter the all signs in the ODE, the results of $p(t)$ seems to be much more reasonable. 
So my conjecture now would be that we mess up the sign in the probability arguments somewhere, perhaps in the integration expression of $p(t)$.
But I am having difficulty find the problem, and a second opinion would be really appreciated.