This is the skeleton file for the manuscript with the structure and hyperlinks to other documents/codes/literatures.
The .tex file for the manuscript is [/Manuscript/main.tex](Manuscript/main.tex)

1. Introduction:
	- What is PhenHet
	- Why PhenHet is important
	- Literature of PhenHet
		- Old Papers: [literature.md](literature.md) & [resources.md](resources.md)
		- Todd-Dwyer: (is there a paper that can be refereed?)
		- [Novozhilov2008](./refs/Novozhilov2008.pdf)
		- [RomanescuEtAL(2023)](https://doi.org/10.1016/j.epidem.2023.100708)
		- Their limitation/flaws
2. Random Network Approaches
	- Introduction to MSV type of random network framework
		- [MillerSlimVolz2011](./refs/MillerSlimVolz2011.pdf)
		- [ZhaoMagpantay2025]([https://doi.org/10.1002/mma.10963](https://doi.org/10.1002/mma.10963))
		- Section 1 of [JR_NegBinom_Result.md](JR_NegBinom_Result.md)
	- Zhao1 Result: no locality, heterogeneity of in-degree
		- Assumptions (**Need some Help from @JD here**)
		- Section 2 of [JR_NegBinom_Result.md](JR_NegBinom_Result.md)
		- Connection with known results of Novozhilov and Dwyer & Parsons
		- [JD_RZ_curves.R](JD_RZ_curves.R)
	- MFSH (**refer to literature**): no locality, heterogeneity of general degree(equal in and out degree)
		- Assumptions
	- MSV configuration model: Heterogeneity of general degree with locality
		- Assumptions: same as MSV config model
		- Section 3 of [JR_NegBinom_Result.md](JR_NegBinom_Result.md)
		- [Initial_Values.md](Initial_Values.md): Discussion about the initial value for $R(0)$ (works well) and $p(0)$ (not working)
		- $\mathcal{R}^*_c$ Result: [NoteForR_c.md](NoteForR_c.md)
			- Comparison to simulation [Simulation.md](Simulation.md) 
		- Discussion of $\mathcal{R}_c$ vs $\mathcal{R}^*_c$: competing infection
		- $\mathcal{R}_i$ Result: [NoteForR_i.md](NoteForR_i.md)
		- $\mathcal{R}_{i,0}>max(d)$
		- Relation between $\mathcal{R}_i$ and $\mathcal{R}_c$ (**Need further investigation**)
		- `simMSV.R` compares the $\mathcal{R}_i$ and $\mathcal{R}_c$ result with simulations. 
3. Summary