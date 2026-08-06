This is the skeleton file for the manuscript with the structure and hyperlinks to other documents/codes/literatures.

The .tex file for the manuscript is [/Manuscript/main.tex](Manuscript/main.tex)

Note to RZ: Update functions of **align** to **aligned** for formulas in .md file .

1. Introduction:
	- What is PhenHet
	- Why PhenHet is important
		- Discussion of $\mathcal{R}_\text{eff}$ and $\mathcal{R}_0$ for PhenHet: [approaches](approaches.md) line 1-16
	- Literature of PhenHet
		- Old Papers: [literature.md](literature.md) & [resources.md](resources.md)
		- Todd-Dwyer: (is there a paper that can be refereed?)
		- [Novozhilov2008](./refs/Novozhilov2008.pdf)
		- [RomanescuEtAL(2023)](https://doi.org/10.1016/j.epidem.2023.100708)
		- Their limitation/flaws
2. Random Network Approaches
	- Introduction to MSV type of random network framework
		- [MillerSlimVolz2011](./refs/MillerSlimVolz2011.pdf)
		- [ZhaoMagpantay2025](https://doi.org/10.1002/mma.10963](https://doi.org/10.1002/mma.10963)
		- Section 1 of [JR_NegBinom_Result.md](JR_NegBinom_Result.md)
3. Zhao1 Result: no locality, heterogeneity of in-degree
	- Assumptions (**Need some help from @JD here to formalize the network model behind Zhao1 result**)
	- Section 2 of [JR_NegBinom_Result.md](JR_NegBinom_Result.md)
	- Connection with known results of Novozhilov and Dwyer & Parsons
	- [JD_RZ_curves.R](JD_RZ_curves.R)
4. Mean Field Social Heterogeneity (**refer to literature: MSV + Novozhilov**): no locality, heterogeneity of general, undirected degree
	- A short discussion for out-degree model
	- A short discussion for in-degree = out-degree, independent model
	- Assumptions: neighbors keep changing instantaneously (i.e. random paring for edges happens all the time)
	- **Optional: Try to derive $\mathcal{R}_\text{eff}$ and $\rho(t)$ result of this model for completeness**
5. MSV configuration model: Heterogeneity of general degree with locality
	- Assumptions: same as MSV config model
	- Section 3 of [JR_NegBinom_Result.md](JR_NegBinom_Result.md)
	- Discussion of $\mathcal{R}_c$ vs $\mathcal{R}_i$, with locality this is different
	- $\mathcal{R}_c$ Section: [NoteForR_c.md](NoteForR_c.md)
		- Zhao2 $\mathcal{R}^*_c$ Result
		- Discussion of $\mathcal{R}_c$ vs $\mathcal{R}^*_c$: competing infection
		- **Need to better understand [Todd's note](outputs/Rc.pdf)** for $p(0)$ working file: [p0Note.md](p0Note.md)
		- Comparison to simulation [Simulation.md](Simulation.md) 
	- $\mathcal{R}_i$ Section: [NoteForR_i.md](NoteForR_i.md)
		- [Initial_Values.md](Initial_Values.md): Discussion about the initial value for $R(0)$ (works well with linear approx.) and $p(0)$ (not working with linear approx.)
		- RZ's $\mathcal{R}_i$ ODE (**Need further investigation**)
		- $\mathcal{R}_{i,0}>max(d)$
	- Relation between $\mathcal{R}_i$ and $\mathcal{R}_c$ (**Need further investigation**)
		- `simMSV.R` compares the $\mathcal{R}_i$ and $\mathcal{R}_c$ result with simulations. 
6. Summary