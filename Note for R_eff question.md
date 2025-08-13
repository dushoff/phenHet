For homogeneous SIR model:
$$
\begin{cases}
    \dot{S}(t)=-\beta S(t)I(t)
    \\
    \dot{I}(t)=+\beta S(t) I(t)-\gamma I(t)
    \\
    \dot{R}(t)=+\gamma I(t)
\end{cases}
$$
One can define the effective reproductive number $\mathcal{R}_{\text{eff}}(t)$ as the (expected/average) new incidence caused by any randomly chosen infected individual in $I(t)$ at the moment $t$ before they recover.

Consider the expected recovery time from the exponential distribution is given by $\frac{1}{\gamma}$, this gives us:
$$\mathcal{R}_{\text{eff}}=-\frac{\dot{S}(t)}{I(t)} \times\frac{1}{\gamma}=\frac{\beta}{\gamma}S(t)$$
for homogeneous model.
If we take the general definition without specifying the incidence term $\dot{S}(t)$: 
$$\mathcal{R}_{\text{eff}}=-\frac{\dot{S}(t)}{I(t)} \times\frac{1}{\gamma}$$it still make sense for the definition.
Also,
$$\mathcal{R}_{\text{eff}}=-\frac{\dot{S(t)}}{I(t)} \times\frac{1}{\gamma}=1 \Leftrightarrow \dot{I}(t)=0$$
which fits the intuitive understanding for the effective reproductive number.

Now consider the system from [J.C. Miller, A.C. Slim & E.M. Volz(2011)](./refs/MillerSlimVolz2011.pdf), such that the system is now a combination of compartment of nodes:
$$
    \begin{cases}
      S(t)=G_p(\theta(t))
      \\
      I(t)=1-S(t)-R(t)
      \\
      \dot{R}(t)=\gamma I(t)
    \end{cases}
$$
and ODE for edges:
$$\dot{\theta}=-\beta\phi_I=-\beta(\theta-\phi_S-\phi_R)=-\beta\theta+\beta\frac{G_p'(\theta)}{\delta}+\gamma(1-\theta)$$
We still have the 