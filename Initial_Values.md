We would like to derive the "eigenvector" initial condition for $R(0)$ and a small enough $\phi(0)$, such that $$\mathcal{R}_{i,0}=\mathcal{R}_i(0)=\max(\mathcal{R}_i)$$ as in [NoteForR_i](NoteForR_i.md).

We start with MSV ODE system:
$$\begin{cases}
    S(t)=G_p(\phi(t))
    \\
    I(t)=1-S(t)-R(t)
    \\
    \dot{R}(t)=\gamma I(t)
\end{cases}$$
$$\dot{\phi}=-\beta\phi_I=-\beta(\phi-\phi_S-\phi_R)=-\beta\phi+\beta\frac{G_p'(\phi)}{\delta}+\gamma(1-\phi)$$where $\phi(t)$ is the probability that a randomly chosen edge has not yet transmitted the disease at time $t$.
At $t=0$, $\phi(0)\rightarrow 1$ so we take 
- $\omega=1-\phi \Leftrightarrow \phi=1-\omega$ 
- $V(t)=1-S(t)=1-G_p(\phi(t))=1-G_p(1-\omega(t))$

Now consider the ODE for $\omega(t)$ and $R(t)$ based on previous system, we have 
$$\begin{cases}
	\dot{\omega}=-\dot{\phi}=\beta(1-\omega)-\beta\frac{G_p'(1-\omega)}{\delta}-\gamma\omega
    \\
    \dot{R}(t)=\gamma (1-G_p(\phi(t))-R(t))
\end{cases}$$

Similar with the derivation of $\mathcal{R}_{c,0}$, using first order approximation, we have:
$$\begin{align}
G_p(\phi(t))&= G_p(1-\omega(t)) = \sum_kp_k(1-\omega)^k
\\
& =\sum_k p_k [1-k\omega+o(\omega^2)]
\\
& \approx \sum_k p_k- \omega\sum_k k p_k
\\
&= 1-\delta\omega
\end{align}$$
and
$$\begin{align}
G'_p(\phi(t))&= G'_p(1-\omega(t)) = \sum_k k p_k(1-\omega)^{k-1}
\\
& =\sum_k k p_k [1-(k-1)\omega+o(\omega^2)]
\\
& \approx \sum_k k p_k- \omega\sum_k k(k-1) p_k
\\
&= \delta(1-\frac{G''_p(1)}{\delta}\times\omega)
\end{align}$$

Using these approximation, the linearized ODE near $t \rightarrow 0$ is
$$\begin{cases}
	\dot{\omega}\approx \beta(1-\omega)-\beta(1-\frac{G''_p(1)}{\delta}\times\omega)-\gamma\omega=[\beta\times\frac{G''_p(1)}{\delta}-(\beta+\gamma)]\omega=\eta\omega
    \\
    \dot{R}(t)\approx\gamma (\delta\omega-R(t))
\end{cases}$$
- Note: $$\eta=\beta\times\frac{G''_p(1)}{\delta}-(\beta+\gamma)=(\beta+\gamma)(\mathcal{R}_{c,0}-1)$$
Then this linearized system has the matrix form:
$$
\begin{bmatrix}
\dot{\omega}
\\
\dot{R}
\end{bmatrix}
=\begin{bmatrix}
\eta & 0
\\
\delta\gamma & -\gamma
\end{bmatrix} \times
\begin{bmatrix}
\omega\\
R
\end{bmatrix} 
$$

The eigenvalues are just $\eta$ and $-\gamma$.

For the dominant eigenvalue$\eta$, the eigenvector satisfy:$$
\begin{bmatrix}
0
\\
0
\end{bmatrix}
=\begin{bmatrix}
0 & 0
\\
\delta\gamma & -\gamma-\eta
\end{bmatrix} \times
\begin{bmatrix}
\omega(0)\\
R(0)
\end{bmatrix} 
\Leftrightarrow
R(0)=\frac{\delta\gamma}{\gamma+\eta}\times\omega(0)
$$
This should give us a initial condition on the eigen-direction.