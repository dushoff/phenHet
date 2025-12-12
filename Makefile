## This is phenomenological heterogeneity
## https://dushoff.github.io/phenHet/JR_NegBinom_Result.html
## https://github.com/dushoff/phenHet/blob/master/outputs/Rnotes.pdf
## https://github.com/dushoff/phenHet/blob/master/outputs/postEdges.Rout.pdf

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

######################################################################

## scaleExamples has mysterious changes that seem to slow it down
scale.Rout: scaleExamples.R scaleFuns.R edgelist.cpp
	$(pipeR)

######################################################################

## Original AI-assisted edgelist pipeline

scaleEdges.Rout: scaleEdges.R scaleFuns.R edgelist.cpp
	$(pipeRcall)

## postEdges.Rout: postEdges.R scaleExamples.R
postEdges.Rout: postEdges.R scaleEdges.rda
	$(pipeRcall)

######################################################################

## Flex pipeline
Sources += $(wildcard slow/*)

impmakeR += params
%.params.Rout: params.R %.params.R
	$(pipeR)

## slow/seed.post.rds: post.R base.params.R
slowtarget/%.post.Rout: post.R %.netsim.rds %.params.rda
	$(pipeR)
impmakeR += post

## seed.plots.Rout: plots.R
## giant.plots.Rout: plots.R
%.plots.Rout: plots.R slow/%.post.rds
	$(pipeR)
impmakeR += plots

## netsim is still outputting a bit data frame
## We're making the netsim post-processing as the slow step
## base.netsim.Rout: netsim.R edgelist.cpp
%.netsim.Rout: netsim.R %.params.rda  scaleFuns.rda edgelist.cpp
	$(pipeR)
impmakeR += netsim

scaleFuns.Rout: scaleFuns.R
	$(wrapR)

######################################################################

## It is much more expensive to save an adjacency list than to make a network, so we're not doing that for now. Not sure if there are work-arounds.

impmakeR += net
## base.net.Rout: net.R base.params.R
%.net.Rout: net.R %.params.rda scaleFuns.rda
	$(pipeR)

impmakeR += sim
## base.sim.Rout: sim.R
%.sim.Rout: sim.R %.net.rds %.params.rda edgelist.cpp
	$(pipeR)

## slowtarget/big.sim.Rout: sim.R big.params.R
slowtarget/%.sim.Rout: sim.R %.net.rds %.params.rda edgelist.cpp
	$(pipeR)

######################################################################


## The separate compilation part has been so painful!
## STOP!!!!!!
%.cpp.Rout: simFun.R %.cpp
	$(pipeR)

## Why does this not work?
compiled.Rout: compiled.R simFun.cpp.rda scaleFuns.R
	$(pipeR)

######################################################################

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

## autopipeR = defined 

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
-include makestuff/slowtarget.mk

-include makestuff/git.mk
-include makestuff/visual.mk
