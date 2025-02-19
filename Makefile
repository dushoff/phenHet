## This is phenomenological heterogeneity

current: target
-include target.mk
Ignore = target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

Sources += $(wildcard *.R)

Sources += pubmed_notes.md

Ignore += secrets.mk
-include secrets.mk
Rentrez.Rout: Rentrez.R
	$(pipeR)

my_fetch.Rout: my_fetch.R
	$(pipeR)

######################################################################

## notes_NovoANDNetwork.html: notes_NovoANDNetwork.Rmd

Rmd = $(wildcard *.Rmd)
Rmdmd = $(Rmd:.Rmd=.md)
Rmdhtml = $(Rmd:.Rmd=.html)

Sources += $(Rmd)

## notes.html: notes.Rmd
%.html: %.Rmd
	$(render_rmd)

Ignore += $(Rmdmd) $(Rmdhtml)

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/00.stamp
makestuff/%.stamp: | makestuff
	- $(RM) makestuff/*.stamp
	cd makestuff && $(MAKE) pull
	touch $@
makestuff:
	git clone --depth 1 $(msrepo)/makestuff

-include makestuff/os.mk

-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
