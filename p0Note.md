Following Todd's [[outputs/Rc.pdf]]/[[Rc.tex]] note:

Linearizing $\ell=-log(\phi) \Leftrightarrow \phi=e^{-\ell}$ near the disease-free equilibrium, s.t. $\ell(0)=\ell_{0}$ and $\ell_0 \rightarrow 0^+$.
Then we have:
$$\begin{split}
-\dot{\ell} &= \frac{d}{dt}[-log(\phi(t))]=\frac{\dot \phi}{\phi}=-\beta+\beta\frac{G'_p(\phi)}{\delta \phi}+\gamma(\frac{1}{\phi}-1)
\\
& = -\beta+\frac{\beta}{\delta}\bigg[\sum_{d=1}^{\infty}d p_d e^{-\ell(d-1)} \bigg]\times e^{\ell}+\gamma(e^{\ell}-1)
\end{split}$$
Take into $e^{\ell}=1+\ell+o(\ell^2)$ and $e^{-\ell}=1-\ell+o(\ell^2)$ near $\ell \rightarrow 0^+$, we further have:
$$\begin{split}
-\dot{\ell} &=-\beta+\frac{\beta}{\delta}\bigg[\sum_{d=1}^{\infty}d p_d (1-\ell+o(\ell^2))^{d-1} \bigg]\times (1+\ell+o(\ell^2))+\gamma\ell
\\
& =-\beta+\frac{\beta}{\delta}\bigg[\sum_{d=1}^{\infty}d p_d - \sum_{d=1}^{\infty}d (d-1) p_d \ell+o(\ell^2) \bigg]\times (1+\ell+o(\ell^2))+\gamma\ell
\\
& =-\beta+\frac{\beta}{\delta}\bigg[\delta - \ell G''_p(1) +o(\ell^2)\bigg]\times (1+\ell+o(\ell^2))+\gamma\ell
\\
& = -\beta+\frac{\beta}{\delta}\bigg[\delta - \ell G''_p(1)+\ell\delta +o(\ell^2) \bigg]+\gamma\ell
\\
& = \ell(- \frac{\beta G''_p(1)}{\delta}+\beta +\gamma )+ o(\ell^2)
\\
& = \ell(\beta+\gamma)(1-\mathcal{R}_{c,0})+ o(\ell^2)
\end{split}$$
We take only the first order term of $\dot{\ell}$ to estimate $\ell(t)$ near $\ell_0$, which provides:
$$\ell(t) \approx \ell_0 e^{(\beta+\gamma)(1-\mathcal{R}_{c,0})\times t}$$
Thus further $$\phi(t)=e^{-\ell(t)}\approx e^{-\ell_0 e^{(\beta+\gamma)(1-\mathcal{R}_{c,0})t}}$$
Since $$p(t)=\int_{0}^{\infty} \beta e^{-(\beta+\gamma)u}\phi_S(t+u)du=\int_{0}^{\infty} \beta e^{-(\beta+\gamma)u}\frac{G'_p(\phi(t+u))}{\delta}du$$, take into the approximation gives us:
$$p(0)= \int_{0}^{\infty} \beta e^{-(\beta+\gamma)u}\frac{G'_p(\phi(u))}{\delta}du \approx \int_{0}^{\infty} \beta e^{-(\beta+\gamma)u}\frac{G'_p(e^{-\ell_0 e^{(\beta+\gamma)(1-\mathcal{R}_{c,0})u}})}{\delta}du$$
Now change variable with 
$$v=\ell_0 e^{(\beta+\gamma)(1-\mathcal{R}_{c,0})u}\Leftrightarrow u=\frac{\log(v)-\log(\ell_0)}{(\beta+\gamma)(1-\mathcal{R}_{c,0})} \Leftrightarrow du=\frac{v^{-1}}{(\beta+\gamma)(1-\mathcal{R}_{c,0})} \times dv$$
The approximation then comes to
$$p(0)\approx \frac{1}{(\beta+\gamma)(1-\mathcal{R}_{c,0}) \delta} \ell_0^{-\frac{1}{\mathcal{R}_{c,0}-1}}\int_{\ell_0}^{\infty}v^{\frac{1}{\mathcal{R}_{c,0}-1}-1}G'_p(e^{-v})dv$$