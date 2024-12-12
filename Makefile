## This is phenomenological heterogeneity
## https://dushoff.github.io/phenHet/Note_Nov14_2024.html

current: target
-include target.mk
Ignore = target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

Rmd = $(wildcard *.Rmd)
Rmdmd = $(Rmd:.Rmd=.md)
Rmdhtml = $(Rmd:.Rmd=.html)

Sources += $(Rmd)

## note.html: note.Rmd
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
