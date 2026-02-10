## This is the MS subdirectory of phenHet

current: target
-include target.mk
Ignore = target.mk

vim_session:
	bash -cl "vmt"

## -include makestuff/perl.def

######################################################################

Sources += $(wildcard *.tex)

## main.pdf: main.tex

######################################################################

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

## ln -s ../makestuff . ## Do this first if you want a linked makestuff
Makefile: makestuff/00.stamp
makestuff/%.stamp: | makestuff
	- $(RM) makestuff/*.stamp
	cd makestuff && $(MAKE) pull
	touch $@
makestuff:
	git clone --depth 1 $(msrepo)/makestuff

-include makestuff/os.mk

-include makestuff/texj.mk

-include makestuff/git.mk
-include makestuff/visual.mk
