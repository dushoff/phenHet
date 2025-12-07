## This is phenomenological heterogeneity
## https://dushoff.github.io/phenHet/JR_NegBinom_Result.html
## https://github.com/dushoff/phenHet/blob/master/outputs/Rnotes.pdf

current: target
-include target.mk
Ignore = target.mk

-include makestuff/perl.def

vim_session:
	bash -cl "vmt todo.md resources.md README.md"

######################################################################

mirrors += resources

Sources += $(wildcard *.md)
Sources += $(wildcard *.tex)
Sources += $(wildcard *.pl)
Sources += $(wildcard *.cpp)

gillespie_tests.Rout: gillespie_tests.R

Ignore += *.html *.pdf
## https://dushoff.github.io/phenHet/approaches.html
approaches.html: approaches.md
	$(rmdh_r)

notes_NovoANDNetwork.html: notes_NovoANDNetwork.md
	$(rmdh_r)

scaleExamples.Rout: scaleExamples.R scaleFuns.R NetSimulator.cpp

scaleFancy.Rout: scaleExamples.R scaleFuns.R fenwick.cpp
	$(pipeR)

scaleEdges.Rout: scaleExamples.R scaleFuns.R edgelist.cpp
	$(pipeR)

NetworkExamples.Rout: NetworkExamples.R NetworkSimulator.R

## Current rcpp implementation
NetworkSimulator.Rout: NetworkSimulator.R

## MRE of sample inconsistency
## mini-EG.Rout: mini-EG.R

######################################################################

## Does this not work? Or only not work off line??
JR_NegBinom_Result.html: JR_NegBinom_Result.md
	$(rmdh_r)

JR_NegBinom_Result.pdf: JR_NegBinom_Result.md
	$(rmdp_r)

JD_RZ_curves.Rout: JD_RZ_curves.R

zhaoFuns.Rout: zhaoFuns.R
zhaoPlot.Rout: zhaoPlot.R zhaoFuns.rda

Ignore += *.MD
## NoteForMu.fix.pdf: NoteForMu.md noobsid.pl
%.fix.MD: %.md noobsid.pl
	$(PUSH)

%.fix.pdf: %.fix.MD
	$(rmdp_r)

NoteForMu.pdf: NoteForMu.md
	$(rmdp_r)

NoteForR_c.pdf: NoteForR_c.md
	$(rmdp_r)

NoteForR_i.pdf:  NoteForR_i.md
	$(rmdp_r)

## NoteForR_i.md.tex:  NoteForR_i.md

## Rnotes.pdf: Rnotes.tex
## Rc.pdf: Rc.tex

Paper.html: Paper.md
	$(rmdh_r)

######################################################################

autopipeR = defined 

## Pubmed stuff for scoping review

Sources += $(wildcard *.R)

Sources += pubmed_notes.md

Ignore += secrets.mk
-include secrets.mk
Rentrez.Rout: Rentrez.R

my_fetch.Rout: my_fetch.R

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/03.stamp
makestuff/%.stamp: | makestuff
	- $(RM) makestuff/*.stamp
	cd makestuff && $(MAKE) pull
	touch $@
makestuff:
	git clone --depth 1 $(msrepo)/makestuff

-include makestuff/os.mk

-include makestuff/pipeR.mk
-include makestuff/mirror.mk
## -include makestuff/rmdweb.mk
-include makestuff/texj.mk
-include makestuff/rmd.mk
-include makestuff/pandoc.mk

-include makestuff/git.mk
-include makestuff/visual.mk
