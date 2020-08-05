#
#       Documentation Makefile for DGELDA
#       Version 1.1
#       July 20, 1995
#
include ../make.inc

all:
	$(LATEX) dgelda; $(LATEX) dgelda; 
	$(LATEX) dnfind; $(LATEX) dnfind; 
	$(LATEX) dnfix1; $(LATEX) dnfix1; 
	$(LATEX) dnfotr; $(LATEX) dnfotr; 
	$(LATEX) dnfred; $(LATEX) dnfred; 
	$(LATEX) dnfval; $(LATEX) dnfval; 
	$(MAKE) clean
	
clean:
	rm -f  *.aux *.bbl *.blg *.log

cleanup: clean
	rm -f  *.dvi dgelda.ps dnfind.ps dnfix1.ps dnfotr.ps \
	       dnfred.ps dnfval.ps	

archiv: cleanup
	-(cd ..; $(MAKE) archiv)

