# $I^q$ case under [Novozhilov](Novozhilov2008.pdf) Framework

Following discussion in [notes_NovoANDNetwork](notes_NovoANDNetwork.Rmd)
-  Assuming the transmission rate $\beta$ is determined by the traits $\omega_s$ and $\omega_i$ independently, such that $\beta=\beta(\omega_s,\omega_i)=\beta_s(\omega_s) \beta_i(\omega_i)$
-   For each traits of S, the density $s(t,\omega_s)$ in trait $\omega_s$ is governed by ODE:
	$$
	\begin{align}
	    \frac{\partial s(t,\omega_s)}{t} & =-s(t,\omega_s)\beta_s(\omega_s) \int_{\Omega_i} \beta_i(\omega_i) i(t,\omega_i) d \omega_i
	    \\
	    &=-s(t,\omega_s)\beta_s(\omega_s)I(t) \int_{\Omega_i} \beta_i(\omega_i) p_i(t,\omega_i) d \omega_i
	    \\
	    &=-s(t,\omega_s)\beta_s(\omega_s)I(t) \bar{\beta_i}(t)
    \end{align}
    $$

	- $\psi(\omega_i,\omega_i')$ is the probability that a newly infected individual gets trait value $\omega_i$ if infected by an individual with trait value $\omega_i'$

- More specifically, this frame work require assumption: the infectivity traits $\omega_i$ of new infected individual must be inherited from its infector, s.t. 
	$$
	\psi(\omega_i,\omega_i')=\delta(\omega_i'-\omega_i)
	$$
with Dirac delta function $\delta(x)$
-   For each traits of $I$, the density $i(t,\omega_i)$ in trait $\omega_i$ is governed by ODE:
$$
\begin{align}
	\frac{\partial i(t,\omega_i)}{t} & =\int_{\Omega_s} \int_{\Omega_i} {s(t,\omega_s) \beta_s(\omega_s) \psi(\omega_i,\omega_i') i(t,\omega_i') \beta_i(\omega_i')} d \omega_i' d \omega_s - \gamma i(t,\omega_i)
	\\
	&=\int_{\Omega_s}s(t,\omega_s)\beta_s(\omega_s) d \omega_s \times \int_{\Omega_i}\delta(\omega_i-\omega_i')i(t,\omega_i')\beta_i(\omega_i')d \omega_i'- \gamma i(t,\omega_i)
	\\
	&=S(t)\int_{\Omega_s}p_s(t,\omega_s)\beta_s(\omega_s)d \omega_s \times \beta_i(\omega_i)i(t,\omega_i)
	\\
	&=S(t) \bar{\beta_s}(t)\beta_i(\omega_i)i(t,\omega_i)- \gamma i(t,\omega_i)
\end{align}
$$

## SI case
Corollary 2 in [Novozhilov2008](Novozhilov2008.pdf)

For each trait $\omega_s$, the PDE for $s(t,\omega_s)$ can be written as
$$
\frac{\partial}{\partial t} s(t,\omega_s)=-\beta_s(\omega_s)s(t,\omega_s)\int_{\Omega_i}\beta_i(\omega_i)i(t,\omega_i)d \omega_i=-\beta_s(\omega_s)s(t,\omega_s)\bar{\beta_i}(t) I(t)
$$
Similarly, for each infected trait $\omega_i$, the PDE for $i(t,\omega_i)$ can be written as 
$$
\frac{\partial}{\partial t} i(t,\omega_i)=\beta_i(\omega_i)i(t,\omega_i)\int_{\Omega_s}\beta_s(\omega_s)s(t,\omega_s)d \omega_s=\beta_i(\omega_i)i(t,\omega_i)\bar{\beta_s}(t) S(t)
$$

where 
$$
\bar{\beta_s}(t)=\int_{\Omega_s} \beta_s(\omega_s)\frac{s(t,\omega_s)}{S(t)} d \omega_i
$$
and
$$
\bar{\beta_i}(t)=\int_{\Omega_i} \beta_i(\omega_i)\frac{i(t,\omega_i)}{I(t)} d \omega_i
$$

By Theorem 1 & Proposition 2 in [Novozhilov2008](Novozhilov2008.pdf), the model is equivalent to 
$$
\begin{align}
	\frac{d}{dt}S(t)&=-h_s(S(t)) h_i(I(t)) 
	\\
	\frac{d}{dt}I(t)&=h_s(S(t)) h_i(I(t))
\end{align}
$$
Here, $h_s$ and $h_i$ is given by 
$$
\begin{align}
	h_s(S)&=S_0[\frac{d}{d\xi}M_s^{-1}(0,\xi)|_{\xi=S(t)/S_0}]^{-1}=S_0 [\frac{d}{d\lambda}M_s(0,\lambda)|_{\lambda=M_s^{-1}(0,S(t)/S_0)}]
	\\
	h_i(S)&=I_0[\frac{d}{d\xi}M_i^{-1}(0,\xi)|_{\xi=I(t)/I_0}]^{-1}=I_0 [\frac{d}{d\lambda}M_i(0,\lambda)|_{\lambda=M_i^{-1}(0,I(t)/I_0)}]
\end{align}
$$
where
- $M_s$ is the MGF of $p_s(0,\omega_s)=\frac{s(0,\omega_s)}{S_0}$
- $M_i$ is the MGF of $p_i(0,\omega_i)=\frac{i(0,\omega_i)}{I_0}$
- Inverse function theorem: 
$$
\frac{d}{dx}(f^{-1})|_{x=b}=\frac{1}{f'(f^{-1})}|_{x=b}
$$

Now assume both susceptibility and infectivity are initially gamma-distributed with corresponding parameter $\alpha$, $\eta$:
- $p_s(0,\omega_s) = \frac{\eta_s^k}{\Gamma(\alpha_s)}\omega_s^{\alpha_s-1} e^{-\eta_s \omega_s}, \omega_s\geq0, \alpha_s>0, \eta_s>0$
- $p_i(0,\omega_i) = \frac{\eta_i^k}{\Gamma(\alpha_i)}\omega_i^{\alpha_i-1} e^{-\eta_i \omega_i}, \omega_i\geq0, \alpha_i>0, \eta_i>0$

==Also, assume that $\beta_s(\omega_s) = \omega_s$ and $\beta_s(\omega_s)=\omega_s$, i.e., the transmission coefficient
takes the values from the domain of traits with the probability corresponding to $p_s$ or $p_i$.==

MGF of gamma-distribution is given by $M(0,\lambda)=(1-\lambda/\eta)^{-\alpha}$. Then inverse MGF at $t=0$ is:
$$
\begin{align}
	&M^{-1}(0,\xi)=\eta (1-\xi^{-\frac{1}{\alpha}})
	\\
	\Rightarrow &\frac{d}{d\xi}M^{-1}(0,\xi)=\frac{\eta}{\alpha}\xi^{-\frac{1}{\alpha}-1}
\end{align}
$$

Take the expression to each case of $s$ and $i$ give us:
$$
\begin{align}
	h_s(S)&=S_0[\frac{\eta_s}{\alpha_s}(\frac{S}{S_0})^{-1/\alpha_s-1}]^{-1}=\frac{\alpha_s}{\eta_s} (S_0)^{-1/\alpha_s} S^{1+\frac{1}{\alpha_s}}=\frac{\alpha_s}{\eta_s} (S_0)^{-1/\alpha_s} S^{p}
\\
	h_i(I)&=I_0[\frac{\eta_i}{\alpha_i}(\frac{I}{I_0})^{-1/\alpha_i0-1}]^{-1}=\frac{\alpha_i}{\eta_i} (I_0)^{-1/\alpha_i} I^{1+\frac{1}{\alpha_i}}=\frac{\alpha_i}{\eta_i} (I_0)^{-1/\alpha_i} I^{q}
\end{align}
$$
where $p=1+\frac{1}{\alpha_s} >1$ and $q=1+\frac{1}{\alpha_i}>1$.

## SIR case
Now for each trait $\omega_i$, we add the recovery term based on trait.
$$
\begin{align}
	\frac{\partial}{\partial t} i(t,\omega_i)&=\beta_i(\omega_i)i(t,\omega_i)\int_{\Omega_s}\beta_s(\omega_s)s(t,\omega_s)d \omega_s-\gamma i(t,\omega_i)
	\\
	&=\beta_i(\omega_i)i(t,\omega_i)\bar{\beta_s}(t) S(t)-\gamma i(t,\omega_i)
	\\
	&=i(t,\omega_i)[-\gamma+\beta_i(\omega_i)\bar{\beta_s}(t)S(t)]
	\\
	&=i(t,\omega_i)[-\gamma+\omega_i \bar{\beta_s}(t)S(t)]
\end{align}
$$

According to Thm 1, the system can be written as:
$$
\begin{align}
	\frac{d}{dt}S(t)&=-S(t)\bar{\beta_s}(t)\bar{\beta_i}(t)I(t)
	\\
	\frac{d}{dt}I(t)&=I(t)[-\gamma+\bar{\beta_i}(t)\bar{\beta_s}(t)S(t)]
\end{align}
$$
with ODEs with auxiliary variable $q_s(t)$ and $q_i(t)$ such that 
$$
\begin{align}
	\frac{d}{dt}q_s(t)&=-\bar{\beta_i}(t)I(t)
	\\
	\frac{d}{dt}q_i(t)&=\bar{\beta_s}(t)S(t)
\end{align}
$$
and given initial condition $q_s(0)=q_i(0)=0$.

These should also satisfies that 
$$
\begin{align}
	\bar{\beta_s}(t)&=\int_{\Omega_s} \omega_s\frac{s(t,\omega_s)}{S(t)} d \omega_s=\frac{\frac{d}{d\lambda}M_s(0,\lambda)}{M_s(0,\lambda)}|_{\lambda=q_s(t)}=\frac{d}{d\lambda}[log(M_s(0,\lambda))]|_{\lambda=q_s(t)}
	\\
	\bar{\beta_i}(t)&=\int_{\Omega_i} \omega_i\frac{i(t,\omega_i)}{i(t)} d \omega_i=\frac{\frac{d}{d\lambda}M_i(0,\lambda)}{M_i(0,\lambda)}|_{\lambda=q_i(t)}=\frac{d}{d\lambda}[log(M_i(0,\lambda))]|_{\lambda=q_i(t)}
\end{align}
$$
This is because: [G.P. Karev(2005)](https://www.aimsciences.org/article/doi/10.3934/proc.2005.2005.487)
$$
M(t,\lambda)=\frac{M(0,\lambda+q(t))}{M(0,q(t))}
$$

Follow previous ODE for $S(t)$
$$
\begin{align}
	&\frac{d}{dt}S(t)=-S(t)\bar{\beta_s}(t)\bar{\beta_i}(t)I(t)
	\\
	\Rightarrow & \frac{1}{S(t)}\frac{d}{dt}S(t)=-\bar{\beta_s}(t)\bar{\beta_i}(t)I(t)=\bar{\beta_s}(t)\frac{d}{dt}q_s(t)
	\\
	\Rightarrow & \frac{d}{dt}log(S(t))=\frac{d}{dt}log(M_s(0,q_s(t)))
\end{align}
$$
Take into the initial condition $S(0)=S_0$, $q_s(0)=0$ gives
$$
S(t)/S_0=M_s(0,q_s(t))
$$

Similarly, for $I(t)$ we have
$$
\begin{align}
	&\frac{d}{dt}I(t)=I(t)[-\gamma+\bar{\beta_i}(t)\bar{\beta_s}(t)S(t)]
	\\
	\Rightarrow & \frac{1}{I(t)}\frac{d}{dt}I(t)=-\gamma+\bar{\beta_i}(t)\frac{d}{dt}q_i(t)
	\\
	\Rightarrow & \frac{d}{dt}log(I(t))=-\gamma+\frac{d}{dt}log(M_i(0,q_i(t)))
\end{align}
$$
Take into the initial condition $I(0)=I_0$, $q_i(0)=0$ gives
$$
I(t)/I_0=e^{-\gamma t}M_i(0,q(t))
$$

Since MGFs are absolutely monotonic function for non negative traits, we got
$$
\begin{align}
	q_s(t)&=M_s^{-1}(0,\frac{S(t)}{S_0})
	\\
	q_i(t)&=M_i^{-1}(0,\frac{I(t)e^{\gamma t}}{I_0})
\end{align}
$$

Now again, assume both susceptibility and infectivity are initially gamma-distributed with MGF $M(0,\lambda)=(1-\lambda/\eta)^{-\alpha}$ and $M^{-1}(0,\xi)=\eta (1-\xi^{-\frac{1}{\alpha}})$
and
$$
\frac{d}{d\lambda}log(M(0,\lambda))=\frac{d}{d\lambda}(-\alpha log(1-\lambda/\eta))=-\alpha \times \frac{1}{1-\lambda/\eta} \times(-\frac{1}{\eta})=\frac{\alpha}{\eta-\lambda}
$$

Then
$$
\begin{align}
	q_s(t)&=M_s^{-1}(0,\frac{S(t)}{S_0})=\eta_s(1-(\frac{S(t)}{S_0})^{-\frac{1}{\alpha_s}})
	\\
	q_i(t)&=M_i^{-1}(0,\frac{I(t)e^{\gamma t}}{I_0})=\eta_i(1-(\frac{I(t)}{I_0})^{-\frac{1}{\alpha_s}}e^{-\frac{\gamma t}{\alpha_s}})
\end{align}
$$
Take these into $\bar\beta_s(t)$ and $\bar\beta_s(t)$ gives us 
$$
\begin{align}
	\bar{\beta_s}(t)&=\frac{d}{d\lambda}[log(M_s(0,\lambda))]|_{\lambda=q_s(t)}=\frac{\alpha_s}{\eta_s-q_s(t)}=\frac{\alpha_s}{\eta_s}(\frac{S(t)}{S_0})^{1/\alpha_s}
	\\
	\bar{\beta_i}(t)&=\frac{d}{d\lambda}[log(M_i(0,\lambda))]|_{\lambda=q_i(t)}=\frac{\alpha_i}{\eta_i-q_i(t)}=\frac{\alpha_i}{\eta_i}(\frac{I(t)e^{\gamma t}}{I_0})^{1/\alpha_i}
\end{align}
$$
Eventually, take this back to ODEs of $S(t)$ and $I(t)$ gives us:
$$
\begin{align}
	\frac{d}{dt}S(t)&=-S(t)\bar{\beta_s}(t)\bar{\beta_i}(t)I(t)=-\frac{\alpha_s\alpha_i}{\eta_s\eta_i}
	(S_0)^{-1/\alpha_s}(I_0)^{-1/\alpha_i}S(t)^{1+1/\alpha_s}I(t)^{1+1/\alpha_i}e^{\frac{\gamma t}{\alpha_i}}=-A S^p I^q e^{\gamma t(q-1)}
	\\
	\frac{d}{dt}I(t)&=I(t)[-\gamma+\bar{\beta_i}(t)\bar{\beta_s}(t)S(t)]=-\gamma I(t)+A S^p I^q e^{\gamma t(q-1)}
\end{align}
$$
where $p=1+\frac{1}{\alpha_s} >1$ and $q=1+\frac{1}{\alpha_i}>1$.

# Exponential Incident

Can we get some specific type of nonlinear incident form from [Novozhilov2008](refs/Novozhilov2008.pdf) framework?

Consider first the "exponential" incident like mentioned in [GranichEtAl2009](https://doi.org/10.1016/s0140-6736(08)61697-9), such the incident term like
$$
\lambda e^{-\alpha I}SI
$$
which let transmission rate decrease with the prevalence exponentially.