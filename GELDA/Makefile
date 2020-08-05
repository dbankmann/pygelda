#
#	Top Level Makefile for DGELDA
#	Version 1.1
#       July 19, 1995
# 
include make.inc

OBJS = dbdstp.o dbdtrp.o dbinom.o dgeerm.o dgeini.o dgelda.o dgenrm.o \
       dnfind.o dnfini.o dnfix1.o dnfotr.o dnfpen.o dnfred.o \
       dnfval.o drkstp.o drkrls.o drkcxs.o 

TARGETS = demo

all: $(TARGETS)

lib: blaslib lapacklib 

docu:
	(cd TeX; $(MAKE))

demo: dgelda.a demo.o
	$(FORTRAN) $(FCOPTS) -o demo demo.o \
	$(LOADOPTS) dgelda.a $(LAPACKLIB) $(BLASLIB) $(TIMELIB)

dgelda.a: 
	(cd SRC; $(MAKE))	

blaslib:	
	(cd BLAS; $(MAKE))

lapacklib:
	(cd LAPACK; $(MAKE))

.f.o:
	$(FORTRAN) $(FCOPTS) -c $<
	
clean:
	rm -f *.o core 
	
cleanup: clean
	rm -f  $(TARGETS) *.a
	(cd SRC; $(MAKE) clean)
	(cd BLAS; $(MAKE) clean)
	(cd INSTALL; $(MAKE) cleanup)
	(cd LAPACK; $(MAKE) clean)
	(cd TeX; $(MAKE) cleanup)

archiv:
	$(MAKE) cleanup
	(cd ..; tar cfv dgelda-1.1.tar dgelda-1.1/*; \
	 compress dgelda-1.1.tar)






