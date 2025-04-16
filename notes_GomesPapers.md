---
title: Notes for Gomes EtAl 2022 papers
author: Richard Zhao
date: 2025-04-14
knit: (function(inputFile, encoding) { rmarkdown::render(inputFile, encoding = encoding, output_dir = "docs") })
bibliography: phenHet.bib
---
### Related papers

- [MGM. Gomes,EtAl.(2022)](https://doi.org/10.1016/j.jtbi.2022.111063) 
- [A. Montalban, RM. Corder and MGM. Gomes(2022)](https://doi.org/10.1007/s00285-022-01771-x)
- [Novozhilov(2008)](./refs/Novozhilov2008.pdf)

### Heterogeneous/Variable Susceptibility Model
- More details in  [A. Montalban, RM. Corder and MGM. Gomes(2022)](https://doi.org/10.1007/s00285-022-01771-x)

- $x$ is the individual susceptibility to infection in relation to the mean
- $q(x)$ is the continuous distribution with mean $\int x q(x) dx=1$.
- $q(x)$ is parameterized by coefficient of variation(CV) $v=\sqrt{\int (x-1)^2 q(x) dx}$ =stdv of $x$.
- $\lambda$ is the average force of infection, such that $$\lambda=\frac{\beta}{N}\int (\rho E(x,t)+I(x,t))dx=\lambda=\frac{\beta}{N}(\rho E(t)+I(t))$$
- $\rho$ is the infectiousness ratio factor of $E$ comparing to $I$
- $\delta$ is the processing rate from $E$ to $I$
- $\phi$ is the proportion of death due to COVID/the disease

They created an SEIR model such that for each "trait" $x$, an equivalent PDE system is given by
$$
\begin{align}
\frac{\partial S(x,t)}{\partial t} & = -\lambda x S(x,t)
\\
\frac{\partial E(x,t)}{\partial t} & = +\lambda x S(x,t)-\delta E(x,t)
\\
\frac{\partial I(x,t)}{\partial t} & = +\delta E(x,t)-\gamma I(x,t)
\\
\frac{\partial R(x,t)}{\partial t} & = +(1-\phi)\gamma I(x,t) 
\end{align}
$$
**In infection term, both connectivity and infectiousness is homogeneous.** And all other quantities is homogeneous.

While the basic reproductive number is given by $$\mathcal{R}_0=\beta(\frac{\rho}{\delta}+\frac{1}{\gamma})$$
Additionally, they defined $\mathcal{R}_c$ to reflect impact of time varying factors (e.g. NPIs, behavior changes, seasonality and pathogen evolution) to reproductive number, such that $$\mathcal{R}_c=c(t)\times\mathcal{R}_0$$
They follow [Novozhilov(2008)](./refs/Novozhilov2008.pdf) and derive the SEIR version of Novozhilov, also assuming that $q(x)$ is Gamma distribution. Then the ODE system for 
$$
\begin{align}
\frac{d S}{dt} & = -\beta (\rho E+I)(\frac{S}{N})^{1+v^2}
\\
\frac{d E}{dt} & = +\beta (\rho E+I)(\frac{S}{N})^{1+v^2}-\delta E
\\
\frac{d I}{dt} & = +\delta E -\gamma I
\\
\frac{d R}{dt} & = +(1-\phi)\gamma I 
\end{align}
$$

