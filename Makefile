# Makefile for SCTK
sinclude ../make.inc

default: all

all: sctk

sctk:
	( cd src ; $(MAKE) all || exit 1 )

clean: sctk_clean #examples_clean

sctk_clean:
	( cd src ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi

distclean: clean doc_clean
