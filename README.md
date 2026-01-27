[approaches.md](approaches.md): A quick summary document of progress

[JR_NegBinom_Result.md](JR_NegBinom_Result.md) ([.pdf](outputs/JR_NegBinom_Result.pdf)): The main background and proceeding document for current MSV random network. RZ is working on summarize recent progress to update this document.   
- Compiled as html here: https://dushoff.github.io/phenHet/JR_NegBinom_Result.html (Hyperlink in .html version not working due to folder issue)

[NoteForR_i.md](NoteForR_i.md)([.pdf](outputs/NoteForR_i.pdf)): Details of $\mathcal{R}_i$ results

[NoteForR_c.md](NoteForR_c.md)([.pdf](outputs/NoteForR_c.pdf)): Details of $\mathcal{R}_c$ and $\mathcal{R}^*_c$ results

[NoteForMu.md](NoteForMu.md): Discussion of expression of $\mu$ in $\mathcal{R}_c$ (need to update JD&TP's discussion email)


[Rc.tex](Rc.tex)([.pdf](Rc.pdf)): Todd version note for ODE of $p(t)$ and $\mathcal{R}_c$ and IV problem for $p(0)$
- [Rc_SIR.tex](Rc_SIR.tex)([.pdf](outputs/Rc_SIR.pdf)) simpler idea for shooting methods
- [p0Note.md](p0Note.md)([.pdf](outputs/p0Note.pdf)): Richard's attempt on $p(0)$ problem

[Initial_Values.md](Initial_Values.md)([.pdf](outputs/Initial_Values.pdf)): Eigen-direction idea for initial values of $R(0)$ and $\mathcal{R}_c(0)$

[Simulation.md](Simulation.md): Summary notes of simulation results.

[/Manuscript/main.tex](/Manuscript/main.tex): Main .tex file for the manuscript. Just templates now.
- [Paper.md](Paper.md): Structure skeleton of the paper
- [resources.md](./resources.md): Summary notes for related literatures and previous investigations:
- [literature.md](literature.md): Summary of literature reviews

[todo](todo.md): Todo list, not being maintained for now
## simulator info

* `NetSimulator.cpp` is a vertex-based simulator 
* `edgelist.cpp` is an edge-based simulator
* `scaleExamples.R` is set up to use *either* machinery, based on what's in the Makefile: compare `scaleEdges.Rout` with `scaleExamples.Rout`. `compareSims.R` (not yet plumbed) is the beginning of machinery to compare the two approaches.



