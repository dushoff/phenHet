## This is phenomenological heterogeneity

current: target
-include target.mk
Ignore = target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt todo.md"

######################################################################

mirrors += resources

Sources += $(wildcard *.md)

## Why I still can't visualize???

## approaches.html: approaches.md

######################################################################

## Pubmed stuff for scoping review

Sources += $(wildcard *.R)

Sources += pubmed_notes.md

Ignore += secrets.mk
-include secrets.mk
Rentrez.Rout: Rentrez.R
	$(pipeR)

my_fetch.Rout: my_fetch.R
	$(pipeR)

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
-include makestuff/mirror.mk

-include makestuff/git.mk
-include makestuff/visual.mk
