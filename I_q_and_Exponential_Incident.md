# I^q case under [Novozhilov](https://arxiv.org/abs/0802.2059) Framework

Following discussion in [notes_NovoANDNetwork](notes_NovoANDNetwork.Rmd) - Assuming the transmission rate $\beta$ is determined by the traits $\omega_s$ and $\omega_i$ independently, such that $\beta=\beta(\omega_s,\omega_i)=\beta_s(\omega_s) \beta_i(\omega_i)$ - For each traits of S, the density $s(t,\omega_s)$ in trait $\omega_s$ is governed by ODE: 
$$
    \begin{align}
        \frac{\partial s(t,\omega_s)}{t} & =-s(t,\omega_s)\beta_s(\omega_s) \int_{\Omega_i} \beta_i(\omega_i) i(t,\omega_i) d \omega_i
        \\
        &=-s(t,\omega_s)\beta_s(\omega_s)I(t) \int_{\Omega_i} \beta_i(\omega_i) p_i(t,\omega_i) d \omega_i
        \\
        &=-s(t,\omega_s)\beta_s(\omega_s)I(t) \bar{\beta_i}(t)
    \end{align}
$$

-   $\psi(\omega_i,\omega_i')$ is the probability that a newly infected individual gets trait value $\omega_i$ if infected by an individual with trait value $\omega_i'$

-   More specifically, this frame work require assumption: the infectivity traits $\omega_i$ of new infected individual must be inherited from its infector, s.t. 
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

Corollary 2 in [Novozhilov2008](https://arxiv.org/abs/0802.2059)

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

By Theorem 1 & Proposition 2 in [Novozhilov2008](https://arxiv.org/abs/0802.2059), the model is equivalent to 
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
    h_i(I)&=I_0[\frac{d}{d\xi}M_i^{-1}(0,\xi)|_{\xi=I(t)/I_0}]^{-1}=I_0 [\frac{d}{d\lambda}M_i(0,\lambda)|_{\lambda=M_i^{-1}(0,I(t)/I_0)}]
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
- $p_s(0,\omega_s) = \frac{1}{\Gamma(\alpha_s)}\eta_s^{\alpha_s} \omega_s^{\alpha_s-1} e^{-\eta_s \omega_s}, \omega_s\geq0, \alpha_s>0, \eta_s>0$ 
- $p_i(0,\omega_i) = \frac{1}{\Gamma(\alpha_i)}\eta_i^{\alpha_i}\omega_i^{\alpha_i-1} e^{-\eta_i \omega_i}, \omega_i\geq0, \alpha_i>0, \eta_i>0$

==Also, assume that $\beta_s(\omega_s) = \omega_s$ and $\beta_i(\omega_i)=\omega_i$, i.e., the transmission coefficient takes the values from the domain of traits with the probability corresponding to $p_s$ or $p_i$.==

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

Note that based on property of gamma distribution, mean of the initial susceptibility distribution $\bar{\beta_s}(0)$ and mean of the initial infectivity distribution $\bar{\beta_i}(0)$ at $t=0$ is just: 
$$
\begin{align}
\bar{\beta_s}(0)&=\int_{\Omega_s}\omega_s p(0,\omega_s)d\omega_s=\frac{\alpha_s}{\eta_s}
\\
\bar{\beta_i}(0)&=\int_{\Omega_i}\omega_i p(0,\omega_i)d\omega_i=\frac{\alpha_i}{\eta_i}
\end{align}
$$ 
Therefore, we further have: 
$$
\begin{align}
    h_s(S)&=\frac{\alpha_s}{\eta_s} (S_0)^{-1/\alpha_s} S^{p}=\bar{\beta_s}(0) (S_0)^{1-p} S^{p}
\\
    h_i(I)&=\frac{\alpha_i}{\eta_i} (I_0)^{-1/\alpha_i} I^{q}=\bar{\beta_i}(0) (I_0)^{1-q} I^{q}
\end{align}
$$

Eventually, the homogeneous version of the SI model is given by: 
$$
\begin{align}
    \frac{d}{dt}S(t)&=-h_s(S(t)) h_i(I(t))= -\bar{\beta_s}(0) \bar{\beta_i}(0) (S_0)^{1-p} (I_0)^{1-q} S^{p} I^{q}
    \\
    \frac{d}{dt}I(t)&=h_s(S(t)) h_i(I(t))=\bar{\beta_s}(0) \bar{\beta_i}(0) (S_0)^{1-p} (I_0)^{1-q} S^{p} I^{q}
\end{align}
$$

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

Now again, assume both susceptibility and infectivity are initially gamma-distributed with MGF $M(0,\lambda)=(1-\lambda/\eta)^{-\alpha}$ and $M^{-1}(0,\xi)=\eta (1-\xi^{-\frac{1}{\alpha}})$ and
$$
\frac{d}{d\lambda}log(M(0,\lambda))=\frac{d}{d\lambda}(-\alpha log(1-\lambda/\eta))=-\alpha \times \frac{1}{1-\lambda/\eta} \times(-\frac{1}{\eta})=\frac{\alpha}{\eta-\lambda}
$$

Then 
$$
\begin{align}
    q_s(t)&=M_s^{-1}(0,\frac{S(t)}{S_0})=\eta_s(1-(\frac{S(t)}{S_0})^{-\frac{1}{\alpha_s}})
    \\
    q_i(t)&=M_i^{-1}(0,\frac{I(t)e^{\gamma t}}{I_0})=\eta_i(1-(\frac{I(t)e^{\gamma t}}{I_0})^{-\frac{1}{\alpha_i}})
\end{align}
$$ 
Take these into $\bar\beta_s(t)$ and $\bar\beta_s(t)$ gives us 
$$
\begin{align}
    \bar{\beta_s}(t)&=\frac{d}{d\lambda}[log(M_s(0,\lambda))]|_{\lambda=q_s(t)}=\frac{\alpha_s}{\eta_s-q_s(t)}=\frac{\alpha_s}{\eta_s}(\frac{S(t)}{S_0})^{1/\alpha_s}=\bar{\beta_s}(0)(\frac{S(t)}{S_0})^{1/\alpha_s}
    \\
    \bar{\beta_i}(t)&=\frac{d}{d\lambda}[log(M_i(0,\lambda))]|_{\lambda=q_i(t)}=\frac{\alpha_i}{\eta_i-q_i(t)}=\frac{\alpha_i}{\eta_i}(\frac{I(t)e^{\gamma t}}{I_0})^{1/\alpha_i}=\bar{\beta_i}(0)(\frac{I(t) e^{\gamma t}}{I_0})^{1/\alpha_i}
\end{align}
$$

Eventually, take this back to ODEs of $S(t)$ and $I(t)$ gives us: 
$$
\begin{align}
    \frac{d}{dt}S(t)&=-S(t)\bar{\beta_s}(t)\bar{\beta_i}(t)I(t)
    \\
    &=-\frac{\alpha_s\alpha_i}{\eta_s\eta_i}
    (S_0)^{-1/\alpha_s}(I_0)^{-1/\alpha_i}S(t)^{1+1/\alpha_s}I(t)^{1+1/\alpha_i}e^{\frac{\gamma t}{\alpha_i}}
    \\
    &=-\bar{\beta_s}(0)\bar{\beta_i}(0)(S_0)^{p-1}(I_0)^{q-1}S(t)^p I(t)^q e^{\gamma t (q-1)}
    \\
    &= -B(\alpha_s, \eta_s,\alpha_i,\eta_i;S_0,I_0) \times S^p I^q e^{\gamma t (q-1)}
    \\
    \frac{d}{dt}I(t)&=I(t)[-\gamma+\bar{\beta_i}(t)\bar{\beta_s}(t)S(t)]
    \\
    &=-\gamma I(t)+\bar{\beta_s}(0)\bar{\beta_i}(0)(S_0)^{p-1}(I_0)^{q-1}S(t)^p I(t)^q e^{\gamma t (q-1)}
    \\
    &=-\gamma I(t)+ B(\alpha_s, \eta_s,\alpha_i, \eta_i;S_0,I_0) \times S^p I^q e^{\gamma t(q-1)}
\end{align}
$$ 
where $p=1+\frac{1}{\alpha_s} >1$ and $q=1+\frac{1}{\alpha_i}>1$.

==Unlike the SI model, the SIR model generates the time dependent exponential term $e^{\gamma t(q-1)}$ which indicate that transmission will increase to infinity as time goes by. This comes from the assumptions:== 
- ==Unbounded infectiousness (Gamma Distribution)==
- ==Perfect fidelity/inheritance to infectiousness during transmission==
==Higher (to infinity) infectiousness traits will be dominant as time goes by.==

==This makes this framework much less useful for $I^q$ and perhaps other nonlinear $I$ case.==

## Basic Reproduction Ratio/Number $R_0$

Note for the $I^q$ case, initial condition $S_0$ and $I_0$ are part of the constant $B$, such that 
$$ 
B=\frac{\alpha_s\alpha_i}{\eta_s\eta_i}(S_0)^{1-p}(I_0)^{1-q}=\bar{\beta_s}(0)\bar{\beta_i}(0)(S_0)^{1-p}(I_0)^{1-q}
$$

As $t \rightarrow 0$, $I(t) \rightarrow I(0)=I_0 \rightarrow 0^{+}$ and $S(t)\rightarrow S(0)=S_0 \rightarrow 1^{-}$. 
When calculating 
$$
\begin{align}
R_{\text{eff}} & \approx \lim_{t \rightarrow 0}\frac{BS(t)^p I(t)^qe^{\gamma t (q-1)}}{I(t)}=B \lim_{t \rightarrow 0}S(t)^p I(t)^{q-1}e^{\gamma t (q-1)}
\\
&=\frac{\alpha_s\alpha_i}{\eta_s\eta_i}(S_0)^{1-p}(I_0)^{1-q}\lim_{t \rightarrow 0} S(t)^p I(t)^{q-1}e^{\gamma t (q-1)}
\\
&=\frac{\alpha_s\alpha_i}{\eta_s\eta_i}(S_0)^{1-p}(I_0)^{1-q} (S_0)^p (I_0)^{q-1}\lim_{t\rightarrow0}e^{\gamma t (q-1)}
\\
&=\frac{\alpha_s\alpha_i}{\eta_s\eta_i}S_0=\bar{\beta_s}(0)\bar{\beta_i}(0)S_0
\end{align} $$ 
Therefore, the $p,q$ value won't actually change the calculation of reproduction number. 
$R_{\text{eff}}$ is determined by the mean of initial gamma distribution and initial condition for $S_0$.

==This shed some lights on how we should construct nonlinearity for $I$. It should contain some form with $f(I/I_0)$ to avoid the absurd $R_0$ issue we previously have.==

To more rigorously derive $R_0$, I take the interpretation and derivation from [DiekmannHeesterbeekMetz(1990)](https://doi-org.libaccess.lib.mcmaster.ca/10.1007/BF00178324).

*The basic reproduction ratio $R_0$ is the expected number of secondary cases produced in a completely susceptible population, by a typical infected individual during its entire period of infectiousness.*

Define $A(\tau, \omega_s, \omega_i)$ be the expected infectivity of an individual, which was infected $\tau$ units of time ago, while having infectivity trait $\omega_i$ towards a susceptible individual with susceptible trait $\omega_s$. 
Based on previous assumptions (include exponential recovery time distribution with rate $\gamma$) and derivation, we have: 
$$
A(\tau, \omega_s, \omega_i)=\beta(\omega_s,\omega_i)e^{-\gamma \tau}=\beta_s(\omega_s)\beta_i(\omega_i)e^{-\gamma \tau}=\omega_s \omega_i e^{-\gamma \tau}
$$

Consider the independent gamma distribution with pdf $p_s(0,\omega_s), p_i(0, \omega_i)$ for the initial time $t=0$ should be the same as the distribution of completely susceptible population, we have: 
$$
\begin{align}
R_0 & =\int_{\Omega_s}p_s(0,\omega_s)\int_{\Omega_i} p_i(0,\omega_i)\int_0^{\infty}A(\tau,\omega_s, \omega_i)d\tau d\omega_i d\omega_s
\\
& = \int_{\Omega_s}p_s(0,\omega_s) \int_{\Omega_i} p_i(0,\omega_i)\int_0^{\infty}\omega_s \omega_i e^{-\gamma \tau}d\tau d\omega_i d\omega_s
\\
& =\int_{\Omega_s} \omega_s p_s(0,\omega_s) \int_{\Omega_i} \omega_i p_i(0,\omega_i)\int_0^{\infty} e^{-\gamma \tau}d\tau d\omega_i d\omega_s
\\
&=\frac{\bar{\beta_s}(0) \bar{\beta_i}(0)}{\gamma}
\end{align}
$$

This agree with the homogeneous SIR case and Novozhilov's special result without heterogeneity in infectivity. 
The invasion of the disease will then only determined by the mean of the heterogeneity distribution, but independent with variance.

The variance of the heterogeneity will affect outbreak dynamics as pointed out by Proposition 3 in [Novozhilov2008](https://arxiv.org/abs/0802.2059), for general case without specifying distribution. 
Higher variance in susceptibility lower the severity of early outbreak progression: $S(t)$ gets larger for $t \in(0,\epsilon)$ as variance of susceptibility increase.

Without explicit proof, [Novozhilov2008](https://arxiv.org/abs/0802.2059) also claim that larger variance of infectivity increase the severity of early outbreak.

==(TO DO?) Final infection size under this framework==

# Exponential Incident

Can we get some specific type of nonlinear incident form from [Novozhilov2008](https://arxiv.org/abs/0802.2059) framework?

Consider first the "exponential" incident like mentioned in [GranichEtAl2009](https://doi.org/10.1016/s0140-6736(08)61697-9), such the incident term like $$\lambda e^{-\alpha I}SI$$ which let transmission rate decrease with the prevalence exponentially.

An ansatz for now is we would like to have homogeneous in susceptibility and heterogenous infectivity with $p_i(0,\omega_i)$ such that corresponding MGF $M(0,\lambda)$ and corresponding $M^{-1}(0,\xi)$ satisfy:

$$
\begin{align}
h(I) =I_0 \frac{d}{d\lambda}M_i(0,\lambda)\bigg|_{\lambda=M^{-1}(0,I/I_0)}&=C I e^{-\alpha I}
\\
I_0 (\frac{d}{d\xi}M_i^{-1}(0,\xi)\bigg|_{\xi=I/I_0})^{-1}&=C I e^{-\alpha I}
\\
(\frac{d}{d\xi}M_i^{-1}(0,\xi)\bigg|_{\xi=I/I_0})^{-1}& =C\times(\frac{I}{I_0})e^{-\alpha(\frac{I}{I_0})} \quad \leftarrow \alpha \text{ absorb } I_0
\\
\frac{d}{d\xi}M_i^{-1}(0,\xi)\bigg|_{\xi=I/I_0}&=C^{-1} (\frac{I}{I_0})^{-1}e^{\alpha(\frac{I}{I_0})}
\\
\frac{d}{d\xi} M_i^{-1}(0,\xi)&= \frac{e^{\alpha \xi}}{C\xi} 
\\
M_i^{-1}(0,\xi)&=C^{-1}\int \frac{e^{\alpha \xi}}{\xi}  d\xi + K= C^{-1} \int \frac{e^{\alpha \xi}}{ \alpha\xi}  d(\alpha\xi)+K
\end{align}
$$

==(TO DO?) Difficulty here: integration and convert back to probability distribution from mgf $M$.==

[Stackexchange](https://math.stackexchange.com/questions/251795/problem-when-integrating-ex-x)

[Wikipedia](https://en.wikipedia.org/wiki/Exponential_integral)


# Other type of non-linear incident

[HethcoteDriessche(1991)](https://link.springer.com/article/10.1007/BF00160539): $$ SI \frac{I^p}{1+\alpha I^q}$$
